'''
dedup.py - Deduplicate reads that are coded with a UMI
=========================================================

:Author: Ian Sudbery
:Release: $Id$
:Date: |today|
:Tags: Python UMI

Purpose
-------

The purpose of this command is to deduplicate BAM files based
on the first mapping co-ordinate and the UMI attached to the read.
It is assumed that the FASTQ files were processed with extract_umi.py
before mapping and thus the UMI is the last word of the read name. e.g:

@HISEQ:87:00000000_AATT

where AATT is the UMI sequeuence.

If you have used an alternative method which does not separate the
read id and UMI with a "_", such as bcl2fastq which uses ":", you can
specify the separator with the option "--umi-separator=<sep>",
replacing <sep> with e.g ":".

Alternatively, if your UMIs are encoded in a tag, you can specify this
by setting the option --extract-umi-method=tag and set the tag name
with the --umi-tag option. For example, if your UMIs are encoded in
the 'UM' tag, provide the following options:
"--extract-umi-method=tag --umi-tag=UM"

By default, reads are considered identical if they have the same start
coordinate, are on the same strand, and have the same UMI. Optionally,
splicing status can be considered (see below).

The start postion of a read is considered to be the start of its alignment
minus any soft clipped bases. A read aligned at position 500 with
cigar 2S98M will be assumed to start at postion 498.

Methods
-------

dedup can be run with multiple methods to identify group of reads with
the same (or similar) UMI(s), from which a single read is
returned. All methods start by identifying the reads with the same
mapping position.

The simpliest methods, unique and percentile, group reads with
the exact same UMI. The network-based methods, cluster, adjacency and
directional, build networks where nodes are UMIs and edges connect UMIs
with an edit distance <= threshold (usually 1). The groups of reads
are then defined from the network in a method-specific manner. For all
the network-based methods, a single read is returned with the most
abundant UMI in each read group. For details about how the read is
selected, see 'Selecting the representative read' below.

  "unique"
      Reads group share the exact same UMI

  "percentile"
      Reads group share the exact same UMI. UMIs with counts < 1% of the
      median counts for UMIs at the same position are ignored.

  "cluster"
      Identify clusters of connected UMIs (based on hamming distance
      threshold). Each network is a read group

  "adjacency"
      Cluster UMIs as above. For each cluster, select the node(UMI)
      with the highest counts. Visit all nodes one edge away. If all
      nodes have been visted, stop. Otherise, repeat with remaining
      nodes until all nodes have been visted. Each step
      defines a read group.

  "directional"
      Identify clusters of connected UMIs (based on hamming distance
      threshold) and umi A counts >= (2* umi B counts) - 1. Each
      network is a read group.

Selecting the representative read
---------------------------------

The following criteria are applied to select the read that will be retained
from a group of duplicated reads:

1. The read with the lowest number of mapping coordinates (see
--multimapping-detection-method option)
2. The read with the highest mapping quality

Otherwise a read is chosen at random.


Options
-------
--extract-umi-method (choice)
      How are the UMIs encoded in the read?

      Options are:

      - "read_id" (default)
            UMIs contained at the end of the read separated as
            specified with --umi-separator option

      - "tag"
            UMIs contained in a tag, see --umi-tag option

--umi-separator (string)
      Separator between read id and UMI. See --extract-umi-method above

--umi-tag (string)
      Tag which contains UMI. See --extract-umi-method above

--method (string, choice)
      Method used to identify PCR duplicates within reads. All methods
      start by identifying the reads with the same mapping position

      Options are:

      - "unique"

      - "percentile"

      - "cluster"

      - "adjacency"

      - "directional" (default)

--edit-distance-threshold (int)
       For the adjacency and cluster methods the threshold for the
       edit distance to connect two UMIs in the network can be
       increased. The default value of 1 works best unless the UMI is
very long (>14bp)

--paired
       BAM is paired end - output both read pairs. This will also
       force the Use of the template length to determine reads with
       the same mapping coordinates.

--spliced-is-unique
       Causes two reads that start in the same position on the same
       strand and having the same UMI to be considered unique if one is spliced
       and the other is not. (Uses the 'N' cigar operation to test for
       splicing)

--soft-clip-threshold (int)
       Mappers that soft clip, will sometimes do so rather than mapping a
       spliced read if there is only a small overhang over the exon
       junction. By setting this option, you can treat reads with at least
       this many bases soft-clipped at the 3' end as spliced.

--multimapping-detection-method (string, choice)
       If the sam/bam contains tags to identify multimapping reads, you can
       specify for use when selecting the best read at a given loci.
       Supported tags are "NH", "X0" and "XT". If not specified, the read
       with the highest mapping quality will be selected

--read-length
      Use the read length as as a criteria when deduping, for e.g sRNA-Seq

--whole-contig (string)
      Consider all alignments to a single contig together. This is useful if
      you have aligned to a transcriptome multi-fasta

--output-stats (string, filename_prefix)
       Output edit distance statistics and UMI usage statistics
       using this prefix.

       Output files are:

       "[prefix]_stats_per_umi_per_position.tsv"
           Histogram of counts per position per UMI pre- and post-deduplication

       "[prefix]_stats_per_umi_per.tsv"
           Table of stats per umi. Number of times UMI was observed,
           total counts and median counts, pre- and post-deduplication

       "[prefix]_stats_edit_distance.tsv"
           Edit distance between UMIs at each position. Positions with a
           single UMI are reported seperately. Pre- and post-deduplication and
           inluding null expectations from random sampling of UMIs from the
           UMIs observed across all positions.

--further-stats (string, filename_prefix)
       Output additional statistics on the toplogies of the UMI
       clusters. Note: each position may contain multiple clusters of
       connected UMIs(nodes). Further stats are only possible when
       using "adjacency" or "cluster" methods

       Output files use the same prefix as above:

       "[prefix]_stats_topologies.tsv"
           The number of clusters with a single node, a single "hub"
           node connected to all other nodes or a more "complex"
           topology

       "[prefix]_stats_nodes.tsv"
           Hisogram of the number of nodes per cluster

--subset (float, [0-1])
      Only consider a fraction of the reads, chosen at random. This is useful
      for doing saturation analyses.

--chrom (string)
      Only consider a single chromosome. This is useful for debugging purposes

--per-contig (string)
      Deduplicate per contig. All reads with the same contig will be
      considered to have the same alignment position. This is useful
      if your library prep generates PCR duplicates with non identical
      alignment positions such as CEL-Seq. In this case, you would
      align to a reference transcriptome with one transcript per gene

--per-gene (string)
      Deduplicate per gene. As above except with this option you can
      align to a reference transcriptome with more than one transcript
      per gene. You need to also provide --gene-transcript-map option

--gene-transcript-map (string)
      File mapping genes to transripts (tab separated), e.g:

      gene1   transcript1
      gene1   transcript2
      gene2   transcript3

--gene-tag (string)
      Deduplicate per gene. As per --per-gene except here the gene
      information is encoded in the bam read tag specified so you do
      not need to supply --gene-transcript-map


-i, --in-sam/-o, --out-sam
      By default, inputs are assumed to be in BAM format and output are output
      in BAM format. Use these options to specify the use of SAM format for
      inputs or outputs.

-I    (string, filename) input file name
      The input file must be sorted and indexed.

-S    (string, filename) output file name

-L    (string, filename) log file name


Usage
-----

    python dedup -I infile.bam -S deduped.bam -L dedup.log


.. note::
   In order to get a valid sam/bam file you need to redirect logging
   information or turn it off logging via -v 0. You can redirect the
   logging to a file with -L <logfile> or use the --log2stderr option
   to send the logging to stderr.

'''
import sys
import collections


# required to make iteritems python2 and python3 compatible
from builtins import dict

from functools import partial

import pysam

import pandas as pd
import numpy as np

try:
    import umi_tools.Utilities as U
except ImportError:
    import Utilities as U

try:
    import umi_tools.network as network
except ImportError:
    import network

try:
    import umi_tools.umi_methods as umi_methods
except ImportError:
    import umi_methods


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
    ''' return a data from with aggregated counts per UMI'''

    agg_df_dict = {}

    agg_df_dict['total_counts'] = stats_df.pivot_table(
        columns="UMI", values="counts", aggfunc=np.sum)

    agg_df_dict['median_counts'] = stats_df.pivot_table(
        columns="UMI", values="counts", aggfunc=np.median)

    agg_df_dict['times_observed'] = stats_df.pivot_table(
        columns="UMI", values="counts", aggfunc=len)

    return pd.DataFrame(agg_df_dict)


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = U.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-i", "--in-sam", dest="in_sam", action="store_true",
                      help="Input file is in sam format [default=%default]",
                      default=False)
    parser.add_option("-o", "--out-sam", dest="out_sam", action="store_true",
                      help="Output alignments in sam format [default=%default]",
                      default=False)
    parser.add_option("--ignore-umi", dest="ignore_umi",
                      action="store_true", help="Ignore UMI and dedup"
                      " only on position", default=False)
    parser.add_option("--umi-separator", dest="umi_sep",
                      type="string", help="separator between read id and UMI",
                      default="_")
    parser.add_option("--umi-tag", dest="umi_tag",
                      type="string", help="tag containing umi",
                      default='RX')
    parser.add_option("--extract-umi-method", dest="get_umi_method", type="choice",
                      choices=("read_id", "tag"), default="read_id",
                      help="where is the read UMI encoded? [default=%default]")
    parser.add_option("--subset", dest="subset", type="float",
                      help="Use only a fraction of reads, specified by subset",
                      default=None)
    parser.add_option("--spliced-is-unique", dest="spliced",
                      action="store_true",
                      help="Treat a spliced read as different to an unspliced"
                           " one [default=%default]",
                      default=False)
    parser.add_option("--soft-clip-threshold", dest="soft",
                      type="float",
                      help="number of bases clipped from 5' end before"
                           "read is counted as spliced [default=%default]",
                      default=4)
    parser.add_option("--edit-distance-threshold", dest="threshold",
                      type="int",
                      default=1,
                      help="Edit distance theshold at which to join two UMIs"
                           "when clustering. [default=%default]")
    parser.add_option("--chrom", dest="chrom", type="string",
                      help="Restrict to one chromosome",
                      default=None)
    parser.add_option("--paired", dest="paired", action="store_true",
                      default=False,
                      help="paired BAM. [default=%default]")
    parser.add_option("--method", dest="method", type="choice",
                      choices=("adjacency", "directional",
                               "percentile", "unique", "cluster"),
                      default="directional",
                      help="method to use for umi deduping [default=%default]")
    parser.add_option("--output-stats", dest="stats", type="string",
                      default=False,
                      help="Specify location to output stats")
    parser.add_option("--further-stats", dest="further_stats",
                      action="store_true", default=False,
                      help="Output further stats")
    parser.add_option("--whole-contig", dest="whole_contig", action="store_true",
                      default=False,
                      help="Read whole contig before outputting bundles: guarantees that no reads"
                           "are missed, but increases memory usage")
    parser.add_option("--multimapping-detection-method",
                      dest="detection_method", type="choice",
                      choices=("NH", "X0", "XT"),
                      default=None,
                      help=("Some aligners identify multimapping using bam "
                            "tags. Setting this option to NH, X0 or XT will "
                            "use these tags when selecting the best read "
                            "amongst reads with the same position and umi "
                            "[default=%default]"))
    parser.add_option("--mapping-quality", dest="mapping_quality",
                      type="int",
                      help="Minimum mapping quality for a read to be retained"
                      " [default=%default]",
                      default=0)
    parser.add_option("--read-length", dest="read_length", action="store_true",
                      default=False,
                      help=("use read length in addition to position and UMI"
                            "to identify possible duplicates [default=%default]"))
    parser.add_option("--per-contig", dest="per_contig", action="store_true",
                      default=False,
                      help=("dedup per contig,"
                            " e.g for transcriptome where contig = gene"))
    parser.add_option("--per-gene", dest="per_gene", action="store_true",
                      default=False,
                      help=("Deduplicate per gene,"
                            "e.g for transcriptome where contig = transcript"
                            "must also provide a transript to gene map with"
                            "--gene-transcript-map [default=%default]"))
    parser.add_option("--gene-transcript-map", dest="gene_transcript_map",
                      type="string",
                      help="file mapping transcripts to genes (tab separated)",
                      default=None)
    parser.add_option("--gene-tag", dest="gene_tag",
                      type="string",
                      help=("Deduplicate per gene where gene is"
                            "defined by this bam tag [default=%default]"),
                      default=None)

    # add common options (-h/--help, ...) and parse command line
    (options, args) = U.Start(parser, argv=argv)

    if options.random_seed:
        np.random.seed(options.random_seed)

    if options.stdin != sys.stdin:
        in_name = options.stdin.name
        options.stdin.close()
    else:
        raise ValueError("Input on standard in not currently supported")

    if options.stdout != sys.stdout:
        out_name = options.stdout.name
        options.stdout.close()
    else:
        out_name = "-"

    if options.in_sam:
        in_mode = "r"
    else:
        in_mode = "rb"

    if options.out_sam:
        out_mode = "w"
    else:
        out_mode = "wb"

    if options.stats:
        if options.ignore_umi:
            raise ValueError("'--output-stats' and '--ignore-umi' options"
                             " cannot be used together")

    if options.further_stats:
        if not options.stats:
            raise ValueError("'--further-stats' options requires "
                             "'--output-stats' option")
        if options.method not in ["cluster", "adjacency"]:
            raise ValueError("'--further-stats' only enabled with 'cluster' "
                             "and 'adjacency' methods")

    if options.per_gene:
        if not options.gene_transcript_map:
            raise ValueError("--per-gene option requires --gene-transcript-map")

    infile = pysam.Samfile(in_name, in_mode)
    outfile = pysam.Samfile(out_name, out_mode, template=infile)

    if options.paired:
        outfile = umi_methods.TwoPassPairWriter(infile, outfile)

    nInput, nOutput = 0, 0

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

    # set the method with which to extract umis from reads
    if options.get_umi_method == "read_id":
        umi_getter = partial(
            umi_methods.get_umi_read_id, sep=options.umi_sep)
    elif options.get_umi_method == "tag":
        umi_getter = partial(
            umi_methods.get_umi_tag, tag=options.umi_tag)
    else:
        raise ValueError("Unknown umi extraction method")

    if options.stats:
        # set up arrays to hold stats data
        stats_pre_df_dict = {"UMI": [], "counts": []}
        stats_post_df_dict = {"UMI": [], "counts": []}
        pre_cluster_stats = []
        post_cluster_stats = []
        pre_cluster_stats_null = []
        post_cluster_stats_null = []
        topology_counts = collections.Counter()
        node_counts = collections.Counter()
        read_gn = umi_methods.random_read_generator(
            infile.filename, chrom=options.chrom, umi_getter=umi_getter)

    if options.chrom:
        inreads = infile.fetch(reference=options.chrom)
    else:
        if options.per_gene:
            metacontig2contig = umi_methods.getMetaContig2contig(
                options.gene_transcript_map)
            inreads = umi_methods.metafetcher(infile, metacontig2contig)
        else:
            inreads = infile.fetch()

    for bundle, read_events, status in umi_methods.get_bundles(
            inreads,
            ignore_umi=options.ignore_umi,
            subset=options.subset,
            quality_threshold=options.mapping_quality,
            paired=options.paired,
            spliced=options.spliced,
            soft_clip_threshold=options.soft,
            per_contig=options.per_contig,
            per_gene=options.per_gene,
            gene_tag=options.gene_tag,
            whole_contig=options.whole_contig,
            read_length=options.read_length,
            detection_method=options.detection_method,
            umi_getter=umi_getter,
            all_reads=False,
            return_unmapped=False):

        nInput += sum([bundle[umi]["count"] for umi in bundle])

        if nOutput % 10000 == 0:
            U.debug("Outputted %i" % nOutput)

        if nInput % 1000000 == 0:
            U.debug("Read %i input reads" % nInput)

        if options.stats:
            # generate pre-dudep stats
            average_distance = umi_methods.get_average_umi_distance(bundle.keys())
            pre_cluster_stats.append(average_distance)
            cluster_size = len(bundle)
            random_umis = read_gn.getUmis(cluster_size)
            average_distance_null = umi_methods.get_average_umi_distance(random_umis)
            pre_cluster_stats_null.append(average_distance_null)

        if options.ignore_umi:
            for umi in bundle:
                nOutput += 1
                outfile.write(bundle[umi]["read"])

        else:

            # set up ReadCluster functor with methods specific to
            # specified options.method
            processor = network.ReadClusterer(options.method)

            # dedup using umis and write out deduped bam
            reads, umis, umi_counts, topologies, nodes = processor(
                bundle=bundle,
                threshold=options.threshold,
                stats=options.stats,
                further_stats=options.further_stats)

            for read in reads:
                outfile.write(read)
                nOutput += 1

            if options.stats:

                # collect pre-dudupe stats
                stats_pre_df_dict['UMI'].extend(bundle)
                stats_pre_df_dict['counts'].extend(
                    [bundle[UMI]['count'] for UMI in bundle])

                # collect post-dudupe stats
                post_cluster_umis = [umi_getter(x) for x in reads]
                stats_post_df_dict['UMI'].extend(umis)
                stats_post_df_dict['counts'].extend(umi_counts)

                average_distance = umi_methods.get_average_umi_distance(post_cluster_umis)
                post_cluster_stats.append(average_distance)

                cluster_size = len(post_cluster_umis)
                random_umis = read_gn.getUmis(cluster_size)
                average_distance_null = umi_methods.get_average_umi_distance(random_umis)
                post_cluster_stats_null.append(average_distance_null)

                if options.further_stats:
                    for c_type, count in topologies.most_common():
                        topology_counts[c_type] += count
                    for c_type, count in nodes.most_common():
                        node_counts[c_type] += count

    outfile.close()

    if options.stats:

        stats_pre_df = pd.DataFrame(stats_pre_df_dict)
        stats_post_df = pd.DataFrame(stats_post_df_dict)

        # generate histograms of counts per UMI at each position
        UMI_counts_df_pre = pd.DataFrame(stats_pre_df.pivot_table(
            columns=stats_pre_df["counts"], values="counts", aggfunc=len))
        UMI_counts_df_post = pd.DataFrame(stats_post_df.pivot_table(
            columns=stats_post_df["counts"], values="counts", aggfunc=len))

        UMI_counts_df_pre.columns = ["instances"]
        UMI_counts_df_post.columns = ["instances"]

        UMI_counts_df = pd.merge(UMI_counts_df_pre, UMI_counts_df_post,
                                 how='left', left_index=True, right_index=True,
                                 sort=True, suffixes=["_pre", "_post"])

        # TS - if count value not observed either pre/post-dedup,
        # merge will leave an empty cell and the column will be cast as a float
        # see http://pandas.pydata.org/pandas-docs/dev/missing_data.html
        # --> Missing data casting rules and indexing
        # so, back fill with zeros and convert back to int
        UMI_counts_df = UMI_counts_df.fillna(0).astype(int)

        UMI_counts_df.to_csv(
            options.stats + "_per_umi_per_position.tsv", sep="\t")

        # aggregate stats pre/post per UMI
        agg_pre_df = aggregateStatsDF(stats_pre_df)
        agg_post_df = aggregateStatsDF(stats_post_df)

        agg_df = pd.merge(agg_pre_df, agg_post_df, how='left',
                          left_index=True, right_index=True,
                          sort=True, suffixes=["_pre", "_post"])

        # TS - see comment above regarding missing values
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

        edit_distance_df = pd.DataFrame({
            "unique": tallyCounts(pre_cluster_binned, max_ed),
            "unique_null": tallyCounts(pre_cluster_null_binned, max_ed),
            options.method: tallyCounts(post_cluster_binned, max_ed),
            "%s_null" % options.method: tallyCounts(post_cluster_null_binned, max_ed),
            "edit_distance": cluster_bins})

        # TS - set lowest bin (-1) to "Single_UMI"
        edit_distance_df['edit_distance'][0] = "Single_UMI"

        edit_distance_df.to_csv(options.stats + "_edit_distance.tsv",
                                index=False, sep="\t")

        if options.further_stats:
            with U.openFile(options.stats + "_topologies.tsv", "w") as outf:
                outf.write(
                    "\n".join(["\t".join((x, str(y)))
                               for x, y in topology_counts.most_common()]) + "\n")

            with U.openFile(options.stats + "_nodes.tsv", "w") as outf:
                outf.write(
                    "\n".join(["\t".join(map(str, (x, y))) for
                               x, y in node_counts.most_common()]) + "\n")

    # write footer and output benchmark information.

    U.info("%s" % ", ".join(
        ["%s: %s" % (x[0], x[1]) for x in read_events.most_common()]))
    U.info("Number of reads out: %i" % nOutput)

    U.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
