'''
dedup.py - Deduplicate reads that are coded with a UMI
=========================================================

:Author: Ian Sudbery
:Release: $Id$
:Date: |today|
:Tags: Python UMI 

Purpose
-------

The perpose of this script is to deduplicate BAM files based
on the first mapping co-ordinate and the UMI attached to the read.
It is assumed that the FASTQ files were processed with extract_umi.py
before mapping and thus the UMI is the last word of the read name.

By default, reads are considered identical if they have the same start
co-ordinate are on the same strand and have the same UMI. Optionally,
splicing status can be considered and reads with similar UMIs can be
removed to account for errors in sequencing (see below).

The start postion of a read is considered to be the start of its alignment
minus any soft clipped bases. A read aligned at position 500 with
cigar 2S98M will be assumed to start at postion 498.

The following criteria are applied to select the read that will be retained
form a group of duplicated reads:

1. The read with the lowest nubmer of hits
2. The read with the highest mapping quality

Otherwise a read is chosen at random.

The input file must be sorted.

.. note::
   In order to get a valid sam/bam file you need to redirect logging
   information or turn it off logging via -v 0. You can redirect the
   logging to a file with -L <logfile> or use the --log2stderr option
   to send the logging to stderr.

Options
-------

--method (choice, string)
      Method used to identify PCR duplicates within reads. All methods
      start by identifying the reads with the same mapping position

      Options are:

      "unique"
          Return all unique UMIs

      "percentile"
          Return all unique UMIs excluding UMIs with counts
          < 1% of the median counts per mapping position

      "cluster"
          Identify clusters of connected UMIs (based on hamming distance
          threshold). Return the UMIs with the highest counts per cluster

      "adjacency"
          Cluster UMIs as above. For each cluster, select and remove
          the node(UMI) with the highest counts and remove all
          connected nodes. Repeat with remaining nodes until no nodes
          remain. Return all selected nodes(UMIs)

      "directional-adjacency"
          Identify clusters of connected UMIs (based on hamming distance
          threshold) and umi A counts > (2* umi B counts) - 1. Return the
          UMIs with the highest counts per cluster

--edit-distance-threshold (int)
       For the adjacency and cluster methods the threshold for the
       edit distance to connect two UMIs in the network can be
       increased. The default value of 1 works best in all situations
       tested to date.

--paired
       Use the template length as a criteria when deduping and output both read
       pairs

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

--multimapping-detection-method (choice, string)
       If the sam/bam contains tags to identify multimapping reads, you can
       specify for use when selecting the best read at a given loci.
       Supported tags are "NH", "X0" and "XT". If not specified, the read
       with the highest mapping quality will be selected

--read-length
      Use the read length as as a criteria when deduping, for e.g sRNA-Seq

--whole-contig
      Consider all alignments to a single contig together. This is useful if
      you have aligned to a transcriptome multi-fasta

--output-stats (filename prefix, string)
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

--further-stats
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

--subset
      Only consider a fraction of the reads, chosen at random. This is useful
      for doing saturation analyses.

--chrom
      Only consider a single chromosome. This is useful for debugging purposes

-i, --in-sam/-o, --out-sam
      By default, inputs are assumed to be in BAM format and output are output
      in BAM format. Use these options to specify the use of SAM format for
      inputs or outputs. SAM input assumed if input from stdin.

Usage
-----

    python dedup_umi.py -I infile.bam -S deduped.bam

Command line options
--------------------

'''
import sys
import random
import collections
import itertools

# required to make iteritems python2 and python3 compatible
from future.utils import iteritems
from builtins import dict

import pysam

import pandas as pd
import numpy as np

import umi_tools.Utilities as U

import pyximport
pyximport.install(build_in_temp=False)
from umi_tools._dedup_umi import edit_distance


def breadth_first_search(node, adj_list):
    searched = set()
    found = set()
    queue = set()
    queue.update((node,))
    found.update((node,))

    while len(queue) > 0:
        node = (list(queue))[0]
        found.update(adj_list[node])
        queue.update(adj_list[node])
        searched.update((node,))
        queue.difference_update(searched)

    return found


def edit_dist(first, second):
    ''' returns the edit distance/hamming distances between
    its two arguments '''

    dist = sum([not a == b for a, b in zip(first, second)])
    return dist


def get_umi(read):
    return read.qname.split("_")[-1]


def get_average_umi_distance(umis):

    if len(umis) == 1:
        return -1

    dists = [edit_distance(x.encode('utf-8'), y.encode('utf-8')) for
             x, y in itertools.combinations(umis, 2)]
    return float(sum(dists))/(len(dists))


def remove_umis(adj_list, cluster, nodes):
    '''removes the specified nodes from the cluster and returns
    the remaining nodes '''

    # list incomprehension: for x in nodes: for node in adj_list[x]: yield node
    nodes_to_remove = set([node
                           for x in nodes
                           for node in adj_list[x]] + nodes)

    return cluster - nodes_to_remove


class ClusterAndReducer:
    '''A functor that clusters a bundle of reads,
    indentifies the parent UMIs and returns the selected reads, umis and counts

    The initiation of the functor defines the methods:

      ** get_adj_list ** - returns the edges connecting the UMIs

      ** connected_components ** - returns clusters of connected components
                                   using the edges in the adjacency list

      ** get_best ** - returns the parent UMI(s) in the connected_components

      ** reduce_clusters ** - loops through the connected components in a
                              cluster and returns the unique reads. Optionally
                              returns lists of umis and counts per umi also

    Note: The get_adj_list and connected_components methods are not required by
    all custering methods. Where there are not required, the methods return
    None or the input parameters.

    '''

    ######## "get_best" methods ##########

    def _get_best_min_account(self, cluster, adj_list, counts):
        ''' return the min UMI(s) need to account for cluster'''
        if len(cluster) == 1:
            return list(cluster)

        sorted_nodes = sorted(cluster, key=lambda x: counts[x],
                              reverse=True)

        for i in range(len(sorted_nodes) - 1):
            if len(remove_umis(adj_list, cluster, sorted_nodes[:i+1])) == 0:
                return sorted_nodes[:i+1]

    def _get_best_higher_counts(self, cluster, counts):
        ''' return the UMI with the highest counts'''
        if len(cluster) == 1:
            return list(cluster)[0]
        else:
            sorted_nodes = sorted(cluster, key=lambda x: counts[x],
                                  reverse=True)
            return sorted_nodes[0]

    def _get_best_percentile(self, cluster, counts):
        ''' return all UMIs with counts >1% of the
        median counts in the cluster '''

        if len(cluster) == 1:
            return list(cluster)
        else:
            threshold = np.median(counts.values())/100
            return [read for read in cluster if counts[read] > threshold]

    def _get_best_null(self, cluster, counts):
        ''' return all UMIs in the cluster'''

        return list(cluster)

    ######## "get_adj_list" methods ##########

    def _get_adj_list_adjacency(self, umis, counts, threshold):
        ''' identify all umis within hamming distance threshold'''

        return {umi: [umi2 for umi2 in umis if
                      edit_distance(umi.encode('utf-8'),
                                    umi2.encode('utf-8')) <= threshold]
                for umi in umis}

    def _get_adj_list_directional_adjacency(self, umis, counts, threshold):
        ''' identify all umis within the hamming distance threshold
        and where the counts of the first umi is > (2 * second umi counts)-1'''

        return {umi: [umi2 for umi2 in umis if
                      edit_distance(umi.encode('utf-8'),
                                    umi2.encode('utf-8')) == 1 and
                      counts[umi] >= (counts[umi2]*2)-1] for umi in umis}

    def _get_adj_list_null(self, umis, counts, threshold):
        ''' for methods which don't use a adjacency dictionary'''
        return None

    ######## "get_connected_components" methods ##########

    def _get_connected_components_adjacency(self, umis, graph, counts):
        ''' find the connected UMIs within an adjacency dictionary'''

        found = list()
        components = list()

        for node in sorted(graph, key=lambda x: counts[x], reverse=True):
            if node not in found:
                component = breadth_first_search(node, graph)
                found.extend(component)
                components.append(component)

        return components

    def _get_connected_components_null(self, umis, adj_list, counts):
        ''' for methods which don't use a adjacency dictionary'''
        return umis

    ######## "reduce_clusters" methods ##########

    def _reduce_clusters_multiple(self, bundle, clusters,
                                  adj_list, counts, stats=False):
        ''' collapse clusters down to the UMI(s) which account for the cluster
        using the adjacency dictionary and return the list of final UMIs'''

        # TS - the "adjacency" variant of this function requires an adjacency
        # list to identify the best umi, whereas the other variants don't
        # As temporary solution, pass adj_list to all variants

        reads = []
        final_umis = []
        umi_counts = []

        for cluster in clusters:
            parent_umis = self.get_best(cluster, adj_list, counts)
            reads.extend([bundle[umi]["read"] for umi in parent_umis])

            if stats:
                final_umis.extend(parent_umis)
                umi_counts.extend([counts[umi] for umi in parent_umis])

        return reads, final_umis, umi_counts

    def _reduce_clusters_single(self, bundle, clusters,
                                adj_list, counts, stats=False):
        ''' collapse clusters down to the UMI which accounts for the cluster
        using the adjacency dictionary and return the list of final UMIs'''

        reads = []
        final_umis = []
        umi_counts = []

        for cluster in clusters:
            parent_umi = self.get_best(cluster, counts)
            reads.append(bundle[parent_umi]["read"])

            if stats:
                final_umis.append(parent_umi)
                umi_counts.append(sum([counts[x] for x in cluster]))

        return reads, final_umis, umi_counts

    def _reduce_clusters_no_network(self, bundle, clusters,
                                    adj_list, counts, stats=False):
        ''' collapse down to the UMIs which accounts for the cluster
        and return the list of final UMIs'''

        reads = []
        final_umis = []
        umi_counts = []

        parent_umis = self.get_best(clusters, counts)
        reads.extend([bundle[umi]["read"] for umi in parent_umis])

        if stats:
            final_umis.extend(parent_umis)
            umi_counts.extend([counts[umi] for umi in parent_umis])

        return reads, final_umis, umi_counts

    def __init__(self, cluster_method="directional-adjacency"):
        ''' select the required class methods for the cluster_method'''

        if cluster_method == "adjacency":
            self.get_adj_list = self._get_adj_list_adjacency
            self.get_connected_components = self._get_connected_components_adjacency
            self.get_best = self._get_best_min_account
            self.reduce_clusters = self._reduce_clusters_multiple

        elif cluster_method == "directional-adjacency":
            self.get_adj_list = self._get_adj_list_directional_adjacency
            self.get_connected_components = self._get_connected_components_adjacency
            self.get_best = self._get_best_higher_counts
            self.reduce_clusters = self._reduce_clusters_single

        elif cluster_method == "cluster":
            self.get_adj_list = self._get_adj_list_adjacency
            self.get_connected_components = self._get_connected_components_adjacency
            self.get_best = self._get_best_higher_counts
            self.reduce_clusters = self._reduce_clusters_single

        elif cluster_method == "percentile":
            self.get_adj_list = self._get_adj_list_null
            self.get_connected_components = self._get_connected_components_null
            self.get_best = self._get_best_percentile
            self.reduce_clusters = self._reduce_clusters_no_network

        if cluster_method == "unique":
            self.get_adj_list = self._get_adj_list_null
            self.get_connected_components = self._get_connected_components_null
            self.get_best = self._get_best_null
            self.reduce_clusters = self._reduce_clusters_no_network

    def __call__(self, bundle, threshold, stats=False, further_stats=False):

        umis = bundle.keys()

        len_umis = [len(x) for x in umis]
        assert max(len_umis) == min(len_umis), (
            "not all umis are the same length(!):  %d - %d" % (
                min(len_umis), max(len(umis))))

        counts = {umi: bundle[umi]["count"] for umi in umis}

        adj_list = self.get_adj_list(umis, counts, threshold)

        clusters = self.get_connected_components(umis, adj_list, counts)

        reads, final_umis, umi_counts = self.reduce_clusters(
            bundle, clusters, adj_list, counts, stats)

        if further_stats:
            topologies = collections.Counter()
            nodes = collections.Counter()

            if len(clusters) == len(umis):
                topologies["single node"] = len(umis)
                nodes[1] = len(umis)
            else:
                for cluster in clusters:
                    if len(cluster) == 1:
                        topologies["single node"] += 1
                        nodes[1] += 1
                    else:
                        most_con = max([len(adj_list[umi]) for umi in cluster])
                        
                        if most_con == len(cluster):
                            topologies["single hub"] += 1
                            nodes[len(cluster)] += 1
                        else:
                            topologies["complex"] += 1
                            nodes[len(cluster)] += 1

        else:
            topologies = None
            nodes = None

        return reads, final_umis, umi_counts, topologies, nodes


def get_bundles(insam, ignore_umi=False, subset=None, quality_threshold=0,
                paired=False, chrom=None, spliced=False, soft_clip_threshold=0,
                per_contig=False, whole_contig=False, read_length=False,
                detection_method="MAPQ"):
    ''' Returns a dictionary of dictionaries, representing the unique reads at
    a position/spliced/strand combination. The key to the top level dictionary
    is a umi. Each dictionary contains a "read" entry with the best read, and a
    count entry with the number of reads with that position/spliced/strand/umi
    combination'''

    last_pos = 0
    last_chr = ""
    reads_dict = collections.defaultdict(
        lambda: collections.defaultdict(
            lambda: collections.defaultdict(dict)))
    read_counts = collections.defaultdict(
        lambda: collections.defaultdict(dict))

    if chrom:
        inreads = insam.fetch(reference=chrom)
    else:
        inreads = insam.fetch()

    for read in inreads:

        if subset:
            if random.random() >= subset:
                continue

        if quality_threshold:
            if read.mapq < quality_threshold:
                continue

        if read.is_unmapped:
            continue

        if read.mate_is_unmapped and paired:
            continue

        if read.is_read2:
            continue

        # TS - some methods require deduping on a per contig
        # (gene for transcriptome) basis, e.g Soumillon et al 2014
        # to fit in with current workflow, simply assign pos and key as contig
        if per_contig:

            pos = read.tid
            key = read.tid
            if not read.tid == last_chr:

                out_keys = reads_dict.keys()

                for p in out_keys:
                    for bundle in reads_dict[p].values():
                        yield bundle
                    del reads_dict[p]
                    del read_counts[p]

                last_chr = read.tid

        else:

            is_spliced = False

            if read.is_reverse:
                pos = read.aend
                if read.cigar[-1][0] == 4:
                    pos = pos + read.cigar[-1][1]
                start = read.pos

                if ('N' in read.cigarstring or
                    (read.cigar[0][0] == 4 and
                     read.cigar[0][1] > soft_clip_threshold)):
                    is_spliced = True
            else:
                pos = read.pos
                if read.cigar[0][0] == 4:
                    pos = pos - read.cigar[0][1]
                start = pos

                if ('N' in read.cigarstring or
                    (read.cigar[-1][0] == 4 and
                     read.cigar[-1][1] > soft_clip_threshold)):
                    is_spliced = True

            if whole_contig:
                do_output = not read.tid == last_chr
            else:
                do_output = start > (last_pos+1000) or not read.tid == last_chr

            if do_output:

                out_keys = [x for x in reads_dict.keys() if x <= start-1000]

                for p in out_keys:
                    for bundle in reads_dict[p].values():
                        yield bundle
                    del reads_dict[p]
                    del read_counts[p]

                last_pos = start
                last_chr = read.tid

            if read_length:
                r_length = read.query_length
            else:
                r_length = 0

            key = (read.is_reverse, spliced & is_spliced,
                   paired*read.tlen, r_length)

        if ignore_umi:
            umi = ""
        else:
            umi = read.qname.split("_")[-1]

        try:
            reads_dict[pos][key][umi]["count"] += 1
        except KeyError:
            reads_dict[pos][key][umi]["read"] = read
            reads_dict[pos][key][umi]["count"] = 1
            read_counts[pos][key][umi] = 0
        else:
            if reads_dict[pos][key][umi]["read"].mapq > read.mapq:
                continue

            if reads_dict[pos][key][umi]["read"].mapq < read.mapq:
                reads_dict[pos][key][umi]["read"] = read
                read_counts[pos][key][umi] = 0
                continue

            # TS: implemented different checks for multimapping here
            if detection_method in ["NH", "X0"]:
                tag = detection_method
                if reads_dict[pos][key][umi]["read"].opt(tag) < read.opt(tag):
                    continue
                elif reads_dict[pos][key][umi]["read"].opt(tag) > read.opt(tag):
                    reads_dict[pos][key][umi]["read"] = read
                    read_counts[pos][key][umi] = 0

            elif detection_method == "XT":
                if reads_dict[pos][key][umi]["read"].opt("XT") == "U":
                    continue
                elif read.opt("XT") == "U":
                    reads_dict[pos][key][umi]["read"] = read
                    read_counts[pos][key][umi] = 0

            read_counts[pos][key][umi] += 1
            prob = 1.0/read_counts[pos][key][umi]

            if random.random() < prob:
                reads_dict[pos][key][umi]["read"] = read

    # yield remaining bundles
    for p in reads_dict:
        for bundle in reads_dict[p].values():
            yield bundle


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


class random_read_generator:
    ''' class to generate umis at random based on the
    distributon of umis in a bamfile '''

    def __init__(self, bamfile, chrom):
        inbam = pysam.Samfile(bamfile)

        if chrom:
            self.inbam = inbam.fetch(reference=chrom)
        else:
            self.inbam = inbam.fetch()

        self.umis = collections.defaultdict(int)
        self.fill()

    def fill(self):

        for read in self.inbam:

            if read.is_unmapped:
                continue

            if read.is_read2:
                continue

            self.umis[get_umi(read)] += 1

        self.observed_umis, freq = zip(*iteritems(self.umis))
        total = sum(freq)
        self.ps = [(x+0.0)/total for x in freq]

    def getUmis(self, n):
        '''get n umis at random'''

        umi_sample = np.random.choice(self.observed_umis, n, p=self.ps)

        return list(umi_sample)


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
                      help="Input file is in sam format", default=False)
    parser.add_option("-o", "--out-sam", dest="out_sam", action="store_true",
                      help="Output alignments in sam format", default=False)
    parser.add_option("--ignore-umi", dest="ignore_umi", action="store_true",
                      help="Ignore UMI and dedup only on position",
                      default=False)
    parser.add_option("--subset", dest="subset", type="string",
                      help="Use only a fraction of reads, specified by subset",
                      default=1.1)
    parser.add_option("--spliced-is-unique", dest="spliced",
                      action="store_true",
                      help="Treat a spliced read as different to an unspliced"
                           " one",
                      default=False)
    parser.add_option("--soft-clip-threshold", dest="soft",
                      type="float",
                      help="number of bases clipped from 5' end before"
                           "read os counted as spliced",
                      default=4)
    parser.add_option("--edit-distance-threshold", dest="threshold",
                      type="int",
                      help="Edit distance theshold at which to join two UMIs"
                           "when clustering", default=1)
    parser.add_option("--chrom", dest="chrom", type="string",
                      help="Restrict to one chromosome",
                      default=None)
    parser.add_option("--paired", dest="paired", action="store_true",
                      default=False,
                      help="Use second-in-pair position when deduping")
    parser.add_option("--method", dest="method", type="choice",
                      choices=("adjacency", "directional-adjacency",
                               "percentile", "unique", "cluster"),
                      default="directional-adjacency",
                      help="method to use for umi deduping")
    parser.add_option("--output-stats", dest="stats", type="string",
                      default=False,
                      help="Specify location to output stats")
    parser.add_option("--further-stats", dest="further_stats",
                      action="store_true", default=False,
                      help="Output further stats")
    parser.add_option("--per-contig", dest="per_contig", action="store_true",
                      default=False,
                      help=("dedup per contig,"
                            " e.g for transcriptome where contig = gene"))
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
                            "amongst reads with the same position and umi"))
    parser.add_option("--mapping-quality", dest="mapping_quality",
                      type="int",
                      help="Minimum mapping quality for a read to be retained",
                      default=0)
    parser.add_option("--read-length", dest="read_length", action="store_true",
                      default=False,
                      help=("use read length in addition to position and UMI"
                            "to identify possible duplicates"))

    # add common options (-h/--help, ...) and parse command line
    (options, args) = U.Start(parser, argv=argv)

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

    infile = pysam.Samfile(in_name, in_mode)
    outfile = pysam.Samfile(out_name, out_mode,
                            template=infile)

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
        read_gn = random_read_generator(infile.filename, chrom=options.chrom)

    for bundle in get_bundles(infile,
                              ignore_umi=options.ignore_umi,
                              subset=float(options.subset),
                              quality_threshold=options.mapping_quality,
                              paired=options.paired,
                              chrom=options.chrom,
                              spliced=options.spliced,
                              soft_clip_threshold=options.soft,
                              per_contig=options.per_contig,
                              whole_contig=options.whole_contig,
                              read_length=options.read_length,
                              detection_method=options.detection_method):

        nInput += sum([bundle[umi]["count"] for umi in bundle])

        if nOutput % 10000 == 0:
            U.debug("Outputted %i" % nOutput)

        if nInput % 1000000 == 0:
            U.debug("Read %i input reads" % nInput)

        if options.stats:
            # generate pre-dudep stats
            average_distance = get_average_umi_distance(bundle.keys())
            pre_cluster_stats.append(average_distance)
            cluster_size = len(bundle)
            random_umis = read_gn.getUmis(cluster_size)
            average_distance_null = get_average_umi_distance(random_umis)
            pre_cluster_stats_null.append(average_distance_null)

        if options.ignore_umi:
            for umi in bundle:
                nOutput += 1
                outfile.write(bundle[umi]["read"])
                if options.paired:
                    try:
                        outfile.write(infile.mate(bundle[umi]["read"]))
                    except:
                        raise ValueError("Mate not found for read: %s" %
                                         bundle[umi]["read"].query_name)
        else:

            # set up ClusterAndReducer functor with methods specific to
            # specified options.method
            processor = ClusterAndReducer(options.method)

            # dedup using umis and write out deduped bam
            reads, umis, umi_counts, topologies, nodes = processor(
                bundle, options.threshold,
                options.stats, options.further_stats)

            for read in reads:
                outfile.write(read)
                nOutput += 1
                if options.paired:
                    # TS - write out paired end mate
                    try:
                        outfile.write(infile.mate(read))
                    except:
                        raise ValueError("Mate not found for read: %s" %
                                         bundle[umi]["read"].query_name)
            if options.stats:

                # collect pre-dudupe stats
                stats_pre_df_dict['UMI'].extend(bundle)
                stats_pre_df_dict['counts'].extend(
                    [bundle[UMI]['count'] for UMI in bundle])

                # collect post-dudupe stats
                post_cluster_umis = [x.qname.split("_")[-1] for x in reads]
                stats_post_df_dict['UMI'].extend(umis)
                stats_post_df_dict['counts'].extend(umi_counts)

                average_distance = get_average_umi_distance(post_cluster_umis)
                post_cluster_stats.append(average_distance)

                cluster_size = len(post_cluster_umis)
                random_umis = read_gn.getUmis(cluster_size)
                average_distance_null = get_average_umi_distance(random_umis)
                post_cluster_stats_null.append(average_distance_null)

                if options.further_stats:
                    for c_type, count in topologies.most_common():
                        topology_counts[c_type] += count
                    for c_type, count in nodes.most_common():
                        node_counts[c_type] += count

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
    U.info("Number of reads in: %i, Number of reads out: %i" %
           (nInput, nOutput))
    U.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
