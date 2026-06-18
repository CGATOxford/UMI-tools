'''
===========================================================
dedup - Deduplicate reads using UMI and mapping coordinates
===========================================================

*Deduplicate reads based on the mapping co-ordinate and the UMI attached to the read*

The identification of duplicate reads is performed in an error-aware
manner by building networks of related UMIs (see
``--method``). ``dedup`` can also handle cell barcoded input (see
``--per-cell``).

Usage::

    umi_tools dedup --stdin=INFILE --log=LOGFILE [OPTIONS] > OUTFILE

Selecting the representative read
---------------------------------
For every group of duplicate reads, a single representative read is
retained.The following criteria are applied to select the read that
will be retained from a group of duplicated reads:

1. The read with the lowest number of mapping coordinates (see
``--multimapping-detection-method`` option)

2. The read with the highest mapping quality. Note that this is not
the read sequencing quality and that if two reads have the same
mapping quality then one will be picked at random regardless of the
read quality.

Otherwise a read is chosen at random.


Dedup-specific options
----------------------
"""""""""""""""""""""""""""
``--output-stats=[PREFIX]``
"""""""""""""""""""""""""""
One can use the edit distance between UMIs at the same position as an
quality control for the deduplication process by comparing with
a null expectation of random sampling. For the random sampling, the
observed frequency of UMIs is used to more reasonably model the null
expectation.

Use this option to generate a stats outfile called:

[PREFIX]_edit_distance.tsv
  Reports the (binned) average edit distance between the UMIs at each
  position. Positions with a single UMI are reported seperately.  The
  edit distances are reported pre- and post-deduplication alongside
  the null expectation from random sampling of UMIs from the UMIs
  observed across all positions. Note that separate null
  distributions are reported since the null depends on the observed
  frequency of each UMI which is different pre- and
  post-deduplication. The post-duplication values should be closer to
  their respective null than the pre-deduplication vs null comparison

In addition, this option will trigger reporting of further summary
statistics for the UMIs which may be informative for selecting the
optimal deduplication method or debugging.

Each unique UMI sequence may be observed [0-many] times at multiple
positions in the BAM. The following files report the distribution for
the frequencies of each UMI.

[PREFIX]_per_umi_per_position.tsv
  The `_per_umi_per_position.tsv` file simply tabulates the
  counts for unique combinations of UMI and position. E.g if prior to
  deduplication, we have two positions in the BAM (POSa, POSb), at
  POSa we have observed 2*UMIa, 1*UMIb and at POSb: 1*UMIc, 3*UMId,
  then the stats file is populated thus:

  ====== =============
  counts instances_pre
  ------ -------------
  1      2
  2      1
  3      1
  ====== =============


  If post deduplication, UMIb is grouped with UMIa such that POSa:
  3*UMIa, then the `instances_post` column is populated thus:

  ====== ============= ==============
  counts instances_pre instances_post
  ------ ------------- --------------
  1      2             1
  2      1             0
  3      1             2
  ====== ============= ==============

[PREFIX]_per_umi.tsv
  The `_per_umi.tsv` table provides UMI-level summary
  statistics. Keeping in mind that each unique UMI sequence can be
  observed at [0-many] times across multiple positions in the BAM,

  :times_observed: How many positions the UMI was observed at
  :total_counts: The total number of times the UMI was observed across all positions
  :median_counts: The median for the distribution of how often the UMI was observed at                  each position (excluding zeros)

  Hence, whenever times_observed=1, total_counts==median_counts.

'''

import sys
import collections
import re
import os

# required to make iteritems python2 and python3 compatible
from builtins import dict

import pysam

import pandas as pd
import numpy as np

import umi_tools
import umi_tools.Utilities as U
import umi_tools.Documentation as Documentation
import umi_tools.network as network
import umi_tools.umi_methods as umi_methods
import umi_tools.sam_methods as sam_methods

# add the generic docstring text
__doc__ = __doc__ + Documentation.GENERIC_DOCSTRING_GDC
__doc__ = __doc__ + Documentation.GROUP_DEDUP_GENERIC_OPTIONS
__doc__ = __doc__ + Documentation.GENERIC_DOCSTRING_SBCRAM_INPUT + Documentation.GENERIC_DOCSTRING_SBCRAM_OUTPUT

usage = '''
dedup - Deduplicate reads using UMI and mapping coordinates

Usage: umi_tools dedup [OPTIONS] [--stdin=IN_BAM] [--stdout=OUT_BAM]

       note: If --stdout is ommited, standard out is output. To
             generate a valid BAM file on standard out, please
             redirect log with --log=LOGFILE or --log2stderr '''


def detect_bam_features(bamfile, n_entries=1000):
    ''' read the first n entries in the bam file and identify the tags
    available detecting multimapping '''

    inbam = pysam.Samfile(bamfile)
    inbam = inbam.fetch(until_eof=True)

    tags = ["NH", "X0", "XT"]
    available_tags = {x: 1 for x in tags}

    for n, read in enumerate(inbam):
        if n > n_entries:
            break

        if read.is_unmapped:
            continue

        else:
            for tag in tags:
                if not read.has_tag(tag):
                    available_tags[tag] = 0

    return available_tags


def aggregateStatsDF(stats_df):
    ''' return a dataframe with aggregated counts per UMI'''

    grouped = stats_df.groupby("UMI")

    agg_dict = {'counts': ['median', len, 'sum']}
    agg_df = grouped.agg(agg_dict)

    agg_df.columns = ['median_counts', 'times_observed', 'total_counts']
    return agg_df


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = U.OptionParser(version="%prog version: $Id$",
                            usage=usage,
                            description=globals()["__doc__"])
    if len(argv) == 1:
        parser.print_usage()
        print ("Required options missing, see --help for more details")
        return 1

    group = U.OptionGroup(parser, "dedup-specific options")

    group.add_option("--output-stats", dest="stats", type="string",
                     default=False,
                     help="Specify location to output stats")

    group.add_option("--stats-sample-fraction", dest="stats_sample_fraction",
                     type="float", default=1.0, metavar="FRACTION",
                     help="Fraction of positions to use when computing edit "
                     "distance statistics (0.0-1.0). Reduces run time for "
                     "large BAMs. Per-UMI count statistics are unaffected "
                     "and always computed from all positions. Default: 1.0")

    group.add_option("--no-per-umi-stats", dest="per_umi_stats",
                     action="store_false", default=True,
                     help="Suppress the per-UMI summary table (*_per_umi.tsv) "
                     "that is written by default when --output-stats is used. "
                     "Useful when many distinct UMI sequences are observed "
                     "(e.g. long UMIs or high sequencing depth).")

    parser.add_option_group(group)

    # add common options (-h/--help, ...) and parse command line
    (options, args) = U.Start(parser, argv=argv, add_dedup_count_sam_options=True,
                              add_s_b_cram_format_options=True)

    U.validateSamOptions(options, group=False)

    if options.random_seed:
        np.random.seed(options.random_seed)

    if options.stdin != sys.stdin:
        in_name = options.stdin.name
        options.stdin.close()
    else:
        raise ValueError("Input on standard in not currently supported")

    in_format = U.determine_format(in_name, options.in_sam, options.in_format)
    out_name, out_format, sorted_out_name, sorted_out_format = U.output_names_and_formats(
        options.stdout,
        options.out_sam,
        options.out_format,
        options.no_sort_output,
        options.tmpdir
    )

    if options.stats and options.ignore_umi:
        raise ValueError("'--output-stats' and '--ignore-umi' options"
                         " cannot be used together")


    infile = U.open_input_alignments(in_name, in_format, options)
    outfile = U.open_output_alignments(out_name, out_format, infile, options)

    if options.paired:
        outfile = sam_methods.TwoPassPairWriter(infile, outfile)

    nInput, nOutput, input_reads, output_reads = 0, 0, 0, 0

    if options.detection_method:
        bam_features = detect_bam_features(infile.filename)

        if not bam_features[options.detection_method]:
            if sum(bam_features.values()) == 0:
                raise ValueError(
                    "There are no bam tags available to detect multimapping. "
                    "Do not set --multimapping-detection-method")
            else:
                raise ValueError(
                    "The chosen method of detection for multimapping (%s) "
                    "will not work with this bam. Multimapping can be detected"
                    " for this bam using any of the following: %s" % (
                        options.detection_method, ",".join(
                            [x for x in bam_features if bam_features[x]])))

    gene_tag = options.gene_tag
    metacontig2contig = None

    if options.chrom:
        inreads = infile.fetch(reference=options.chrom)

    else:
        if options.per_contig and options.gene_transcript_map:
            metacontig2contig = sam_methods.getMetaContig2contig(
                infile, options.gene_transcript_map)
            metatag = "MC"
            inreads = sam_methods.metafetcher(infile, metacontig2contig, metatag)
            gene_tag = metatag

        else:
            inreads = infile.fetch()

    # set up ReadCluster functor with methods specific to
    # specified options.method
    processor = network.ReadDeduplicator(options)

    bundle_iterator = sam_methods.get_bundles(
        options,
        metacontig_contig=metacontig2contig)

    if options.stats:
        # set up arrays to hold stats data
        stats_pre_df_dict = {"UMI": [], "counts": []}
        stats_post_df_dict = {"UMI": [], "counts": []}
        pre_cluster_stats = []
        post_cluster_stats = []
        pre_cluster_stats_null = []
        post_cluster_stats_null = []
        pre_cluster_sizes = []   # sizes of multi-UMI bundles; used post-loop to draw null distributions
        post_cluster_sizes = []
        topology_counts = collections.Counter()
        node_counts = collections.Counter()

    for bundle, key, status in bundle_iterator(inreads):

        nInput += sum([bundle[umi]["count"] for umi in bundle])

        while nOutput >= output_reads + 100000:
            output_reads += 100000
            U.info("Written out %i reads" % output_reads)

        while nInput >= input_reads + 1000000:
            input_reads += 1000000
            U.info("Parsed %i input reads" % input_reads)

        if options.ignore_umi:
            for umi in bundle:
                nOutput += 1
                outfile.write(bundle[umi]["read"])

        else:

            # dedup using umis and write out deduped bam
            reads, umis, umi_counts = processor(
                bundle=bundle,
                threshold=options.threshold)

            # with UMI filtering, it's possible reads
            if len(reads) == 0:
                continue

            for read in reads:
                outfile.write(read)
                nOutput += 1

            if options.stats:

                stats_pre_df_dict['UMI'].extend(bundle)
                stats_pre_df_dict['counts'].extend(
                    [bundle[UMI]['count'] for UMI in bundle])

                stats_post_df_dict['UMI'].extend(umis)
                stats_post_df_dict['counts'].extend(umi_counts)

                # Edit distance computation is O(N^2) per position and
                # dominates run time on large BAMs. Sample a fraction of
                # positions to keep this tractable; per-UMI count stats
                # above are cheap and always collected from all positions.
                if np.random.random() < options.stats_sample_fraction:

                    # Single-UMI positions have no pairs; result and null
                    # are both -1 by convention, with no null draw needed.
                    pre_size = len(bundle)
                    if pre_size == 1:
                        pre_cluster_stats.append(-1)
                        pre_cluster_stats_null.append(-1)
                    else:
                        pre_cluster_stats.append(
                            umi_methods.get_average_umi_distance(bundle.keys()))
                        pre_cluster_sizes.append(pre_size)

                    post_cluster_umis = [bundle_iterator.barcode_getter(x)[0] for x in reads]
                    post_size = len(post_cluster_umis)
                    if post_size <= 1:
                        post_cluster_stats.append(-1)
                        post_cluster_stats_null.append(-1)
                    else:
                        post_cluster_stats.append(
                            umi_methods.get_average_umi_distance(post_cluster_umis))
                        post_cluster_sizes.append(post_size)

    outfile.close()

    if not options.no_sort_output:
        # sort the output
        U.sort_output(sorted_out_name, out_name, sorted_out_format, 
                      options.output_options,
                      options.reference_filename)

    if options.stats:
        # Build the UMI frequency distribution from data already collected
        # during the main loop, then draw null edit distances in one pass.
        # Drawing nulls here rather than per-position avoids repeated
        # construction of the frequency table.
        umi_total_counts = collections.Counter()
        for umi, count in zip(stats_pre_df_dict['UMI'], stats_pre_df_dict['counts']):
            umi_total_counts[umi] += count

        all_umis = list(umi_total_counts.keys())
        total = sum(umi_total_counts.values())
        probs = np.array([umi_total_counts[u] / total for u in all_umis])

        # Re-seed so null draws start from a known state independent of
        # how many random numbers the main loop consumed.
        if options.random_seed:
            np.random.seed(options.random_seed)

        # Draw null UMIs in large batches (one np.random.choice call per
        # buffer-fill) rather than once per position. Each position still gets
        # its own independent slice, preserving the null distribution.
        _NULL_BUFFER_SIZE = 100000
        _random_buf = np.random.choice(all_umis, size=_NULL_BUFFER_SIZE, p=probs)
        _buf_ix = 0

        def _get_null_umis(n):
            nonlocal _random_buf, _buf_ix
            if _buf_ix + n > len(_random_buf):
                _random_buf = np.random.choice(
                    all_umis, size=_NULL_BUFFER_SIZE, p=probs)
                _buf_ix = 0
            umis = _random_buf[_buf_ix:_buf_ix + n].tolist()
            _buf_ix += n
            return umis

        for size in pre_cluster_sizes:
            pre_cluster_stats_null.append(
                umi_methods.get_average_umi_distance(_get_null_umis(size)))

        for size in post_cluster_sizes:
            post_cluster_stats_null.append(
                umi_methods.get_average_umi_distance(_get_null_umis(size)))

    if options.stats:
        # generate the stats dataframe
        stats_pre_df = pd.DataFrame(stats_pre_df_dict)
        stats_post_df = pd.DataFrame(stats_post_df_dict)

        # tally the counts per umi per position
        pre_counts = collections.Counter(stats_pre_df["counts"])
        post_counts = collections.Counter(stats_post_df["counts"])
        counts_index = list(set(pre_counts.keys()).union(set(post_counts.keys())))
        counts_index.sort()
        with U.openFile(options.stats + "_per_umi_per_position.tsv", "w") as outf:
            outf.write("counts\tinstances_pre\tinstances_post\n")
            for count in counts_index:
                values = (count, pre_counts[count], post_counts[count])
                outf.write("\t".join(map(str, values)) + "\n")

        if options.per_umi_stats:
            # aggregate stats pre/post per UMI
            agg_pre_df = aggregateStatsDF(stats_pre_df)
            agg_post_df = aggregateStatsDF(stats_post_df)

            agg_df = pd.merge(agg_pre_df, agg_post_df, how='left',
                              left_index=True, right_index=True,
                              sort=True, suffixes=["_pre", "_post"])

            # TS - if count value not observed either pre/post-dedup,
            # merge will leave an empty cell and the column will be cast as a float
            # see http://pandas.pydata.org/pandas-docs/dev/missing_data.html
            # --> Missing data casting rules and indexing
            # so, back fill with zeros and convert back to int
            agg_df = agg_df.fillna(0).astype(int)

            agg_df.index = [x.decode() for x in agg_df.index]
            agg_df.index.name = 'UMI'
            agg_df.to_csv(options.stats + "_per_umi.tsv", sep="\t")

        # bin distances into integer bins
        max_ed = int(max(map(max, [pre_cluster_stats,
                                   post_cluster_stats,
                                   pre_cluster_stats_null,
                                   post_cluster_stats_null])))

        cluster_bins = range(-1, int(max_ed) + 2)

        def bin_clusters(cluster_list, bins=cluster_bins):
            ''' take list of floats and return bins'''
            return np.digitize(cluster_list, bins, right=True)

        def tallyCounts(binned_cluster, max_edit_distance):
            ''' tally counts per bin '''
            return np.bincount(binned_cluster,
                               minlength=max_edit_distance + 3)

        pre_cluster_binned = bin_clusters(pre_cluster_stats)
        post_cluster_binned = bin_clusters(post_cluster_stats)
        pre_cluster_null_binned = bin_clusters(pre_cluster_stats_null)
        post_cluster_null_binned = bin_clusters(post_cluster_stats_null)

        edit_distance_df = pd.DataFrame(
            {"unique": tallyCounts(pre_cluster_binned, max_ed),
             "unique_null": tallyCounts(pre_cluster_null_binned, max_ed),
             options.method: tallyCounts(post_cluster_binned, max_ed),
             "%s_null" % options.method: tallyCounts(post_cluster_null_binned, max_ed),
             "edit_distance": cluster_bins},
            columns=["unique", "unique_null", options.method,
                     "%s_null" % options.method, "edit_distance"])

        edit_distance_df['edit_distance'] = edit_distance_df['edit_distance'].astype(str)

        # TS - set lowest bin (-1) to "Single_UMI"
        edit_distance_df.loc[0, 'edit_distance'] = "Single_UMI"

        edit_distance_df.to_csv(options.stats + "_edit_distance.tsv",
                                index=False, sep="\t")

    # write footer and output benchmark information.
    U.info(
        "Reads: %s" % ", ".join(["%s: %s" % (x[0], x[1]) for x in
                                 bundle_iterator.read_events.most_common()]))

    U.info("Number of reads out: %i" % nOutput)

    if not options.ignore_umi:  # otherwise processor has not been used
        U.info("Total number of positions deduplicated: %i" %
               processor.UMIClusterer.positions)
        if processor.UMIClusterer.positions > 0:
            U.info("Mean number of unique UMIs per position: %.2f" %
                   (float(processor.UMIClusterer.total_umis_per_position) /
                    processor.UMIClusterer.positions))
            U.info("Max. number of unique UMIs per position: %i" %
                   processor.UMIClusterer.max_umis_per_position)
        else:
            U.warn("The BAM did not contain any valid "
                   "reads/read pairs for deduplication")

    if options.filter_umi:
        U.info("%i UMIs were in a group where the top UMI was not a "
               "whitelist UMI and were therefore discarded"
               % processor.umi_whitelist_counts["Non-whitelist UMI"])

    U.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
