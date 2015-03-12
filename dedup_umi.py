'''
DedupUMI.py - Deduplicate reads that are coded with a UMI
=========================================================

:Author: Ian Sudbery
:Release: $Id$
:Date: |today|
:Tags: Python

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
minus any soft clipped bases. Hence a read aligned at position 100 with
cigar 2S98M will be assumed to start at postion 98.

The following criteria are applied to select the read that will be retained
form a group of duplicated reads:

1. The read with the lowest nubmer of hits
2. The read with the highest mapping quality

Otherwise a read is chosen at random.

The input file must be sorted.

.. note::
   You need to redirect either logging information or output to a file or turn
   it off logging via -v 0 in order to get a valid sam/bam file.

Options
-------

--cluster-umis, --edit-distance-theshold
       Often when looking at reads mapping to a similar base, you will
       find that the umis are more similar than you would expect. This
       option causes the clustering of umis within a threshold edit
       distance of each other, and then a search for the most common
       umis that will explain the cluster.

--spliced-is-unique
       Causes two reads that start in the same position on the same
       strand and having the same UMI to be considered unique if one is spliced
       and the other is not. (Uses the 'N' cigar operation to test for
       splicing)

--soft-clip-threshold
       Mappers that soft clip, will sometimes do so rather than mapping a
       spliced read if there is only a small overhang over the exon
       junction. By setting this option, you can treat reads with at least
       this many bases soft-clipped at the 3' end as spliced.

--paired
       Use the template length as a criteria when deduping. Currently only
       the first in pair read is output, although this might change.

--subset
      Only consider a fraction of the reads, chosen at random. This is useful
      for doing saturation analyses.

--chrom
      Only consider a single chromosome.

-i, --in-sam/-o, --out-sam
      By default, inputs are assumed to be in BAM format and output are output
      in BAM format. Use these options to specify the use of SAM format for
      inputs or outputs.

Usage
-----

    python dedup_umi.py -I infile.bam -O deduped.bam

or

    cat infile.bam | python dedup_umi.py > deduped.bam 2> deduped.log

Command line options
--------------------

'''

import sys
import pysam
import CGAT.Experiment as E
import random
import collections
import itertools
import pandas as pd
import numpy as np


def breadth_first_search(node, adj_list):
    searched = set()
    components = set()
    queue = set()
    queue.update((node,))
    components.update((node,))
    
    while len(queue)>0:
        components.update(adj_list[node])
        queue.update(adj_list[node])
        searched.update((node,))
        queue.difference_update(searched)
        if len(queue)>0:
            node=(list(queue))[0]
    return components


def connected_components(graph, counts):
    found = list()
    components = list()
    for node in sorted(graph, key=lambda x: counts[x], reverse=True):
        if node not in found:
            component = breadth_first_search(node, graph)
            found.extend(component)
            components.append(component)
    return components


def edit_dist(first, second):
    ''' returns the edit distance/hamming distances between
    its two arguements '''

    dist = sum([not a == b for a, b in zip(first, second)])
    return dist


def get_adj_list(umis, counts, threshold=1):
    '''Returns dictionary where each entry represents the nodes
    adjecent to the dictionary key'''
    
    # TS - added threshold for difference in counts
    # adjacency list is therefore directional
    return {umi: [umi2 for umi2 in umis if
                  0 < edit_dist(umi, umi2) <= threshold and
                  counts[umi] >= (counts[umi2]*2)-1] for umi in umis}


def remove_umis(adj_list, cluster, nodes):
    '''removes the specified nodes from the cluster and returns
    the remaining nodes '''

    # list incomprehension: for x in nodes: for node in adj_list[x]: yield node
    nodes_to_remove = set([node
                           for x in nodes
                           for node in adj_list[x]] + nodes)

    return cluster - nodes_to_remove


def get_best(cluster, counts):
    '''Finds the nodes that best explain the cluster by successively
    removing more and more nodes until the cluster is explained.
    Nodes are removed in order of their count. '''
    if len(cluster) == 1:
        return list(cluster)[0]
    else:
        return max(cluster, key=lambda x: counts[x])


def cluster_and_reduce(bundle, threshold=1):
    ''' Recieves a bundle of reads, clusters them by edit distance and counts.
    Within each cluster the UMIs with the highest reads are selected.
    Yeilds the selected reads '''

    umis = bundle.keys()
    counts = {umi: bundle[umi]["count"] for umi in umis}

    adj_list = get_adj_list(umis, counts, threshold)
    clusters = connected_components(adj_list, counts)

    for cluster in clusters:
        parent_umi = get_best(cluster, counts)

        yield bundle[parent_umi]["read"]


def cluster_and_reduce_with_stats(bundle, threshold=1):
    ''' Recieves a bundle of reads, clusters them on edit distance
    and selects connected clusters using the edit distance threshold.
    Within each cluster the UMI with the highest reads is selected.
    Returns the selected reads, umis and counts'''

    umis = bundle.keys()
    counts = {umi: bundle[umi]["count"] for umi in umis}

    adj_list = get_adj_list(umis, counts, threshold)
    clusters = connected_components(adj_list, counts)

    reads = []
    umis = []
    umi_counts = []

    for cluster in clusters:
        parent_umi = get_best(cluster, counts)
        reads.append(bundle[parent_umi]["read"])
        umis.append(parent_umi)
        umi_counts.append(sum([counts[x] for x in cluster]))

    return reads, umis, umi_counts


def get_bundles(insam, ignore_umi=False, subset=None, paired=False,
                chrom=None, spliced=False, soft_clip_threshold=0):
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

        if read.is_unmapped:
            continue

        if read.mate_is_unmapped and paired:
            continue

        if read.is_read2:
            continue

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

        if start > (last_pos+1000) and not read.tid == last_chr:

            out_keys = [x for x in reads_dict.keys() if x <= start-1000]

            for p in out_keys:
                for bundle in reads_dict[p].itervalues():
                    yield bundle
                del reads_dict[p]
                del read_counts[p]
 
            last_pos = start
            last_chr = read.tid

        if ignore_umi:
            umi = ""
        else:
            umi = read.qname.split("_")[-1]

        key = (read.is_reverse, spliced & is_spliced, paired*read.tlen)

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

            if reads_dict[pos][key][umi]["read"].opt("NH") < read.opt("NH"):
                continue
            elif reads_dict[pos][key][umi]["read"].opt("NH") > read.opt("NH"):
                reads_dict[pos][key][umi]["read"] = read
                read_counts[pos][key][umi] = 0

            read_counts[pos][key][umi] += 1
            prob = 1.0/read_counts[pos][key][umi]

            if random.random() < prob:
                reads_dict[pos][key][umi]["read"] = read

    # yeild remaining bundles
    for p in reads_dict:
                for bundle in reads_dict[p].itervalues():
                    yield bundle

class random_read_generator:
    ''' class to generate umis at random based on the 
    distributon of umis in a bamfile '''

    def __init__(self, bamfile):
        inbam = pysam.Samfile(bamfile)
        self.umis = [None] * (inbam.mapped/2)
        self.inbam = inbam.fetch()
        self.max = (inbam.mapped /10) - int((inbam.mapped /2)*0.1)
        self.filled = 0
                
    def fillto(self, n):
        
        for i in range(self.filled, n+1):
            try:
                read = self.inbam.next()
                while read.is_unmapped or read.is_read2:
                    read = self.inbam.next()
            except:
                print i
                raise 
                
            self.umis[i] = get_umi(read)
            
        self.filled = n
        
    def getUmis(self, n):
        '''get n umis at random'''
        
        umis = []
        for umi in range(n):
            
            i = random.randrange(0,self.max)
            if i > self.filled:
                self.fillto(i)
            
            umis.append(self.umis[i])
            
        assert len(umis) == n
        return umis


def aggregateStatsDF(stats_df):
    ''' return a data from with aggregated counts per UMI'''

    agg_df_dict={}

    agg_df_dict['total_counts'] = stats_df.pivot_table(
        columns="UMI", values="counts", aggfunc=np.sum)

    agg_df_dict['median_counts'] = stats_df.pivot_table(
        columns="UMI", values="counts", aggfunc=np.median)

    agg_df_dict['times_observed'] = stats_df.pivot_table(
        columns="UMI", values="counts", aggfunc=len)

    return pd.DataFrame(agg_df_dict)


def get_umi(read):
    return read.qname.split("_")[-1]


def get_average_umi_distance(umis):
    if len(umis) == 1:
        return -1
    dists = [edit_dist(*pair) for pair in itertools.combinations(umis, 2)]
    return float(sum(dists))/(len(dists))


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
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
    parser.add_option("--cluster-umis", dest="cluster_umis",
                      action="store_true",
                      help="Select best reads by clustering umis and removing"
                           "potential sequencing errors", default=False)
    parser.add_option("--edit-distance-theshold", dest="threshold",
                      type="int",
                      help="Edit distance theshold at which to join two UMIs"
                           "when clustering", default=1)
    parser.add_option("--chrom", dest="chrom", type="string",
                      help="Restrict to one chromosome",
                      default=None)
    parser.add_option("--paired", dest="paired", action="store_true",
                      default=False,
                      help="Use second-in-pair position when deduping")
    parser.add_option("--output-stats", dest="stats", type="string",
                      default=None,
                      help="Specify location to output stats")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if options.stdin != sys.stdin:
        in_name = options.stdin.name
        options.stdin.close()
    else:
        in_name = "-"

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

    infile = pysam.Samfile(in_name, in_mode)
    outfile = pysam.Samfile(out_name, out_mode,
                            template=infile)

    nInput, nOutput = 0, 0

    if options.stats and not options.ignore_umi:
        # set up arrays to hold stats data
        stats_pre_df_dict = {"UMI":[], "counts":[]}
        stats_post_df_dict = {"UMI":[], "counts":[]}
        pre_cluster_stats = []
        post_cluster_stats = []
        pre_cluster_stats_null = []
        post_cluster_stats_null = []
        read_gn = random_read_generator(infile.filename)

    for bundle in get_bundles(infile,
                              ignore_umi=options.ignore_umi,
                              subset=float(options.subset),
                              paired=options.paired,
                              chrom=options.chrom,
                              spliced=options.spliced,
                              soft_clip_threshold=options.soft):

        nOutput += 1
        nInput += sum([bundle[umi]["count"] for umi in bundle])

        if nOutput % 10000 == 0:
            E.debug("Outputted %i" % nOutput)

        if nInput % 1000000 == 0:
            E.debug("Read %i input reads" % nInput)

        if options.stats:
            # generate pre-dudep stats
            average_distance = get_average_umi_distance(bundle.keys())
            pre_cluster_stats.append(average_distance)
            cluster_size = len(bundle)
            random_umis = read_gn.getUmis(cluster_size)
            average_distance_null = get_average_umi_distance(random_umis)
            pre_cluster_stats_null.append(average_distance_null)
            
        if options.ignore_umi or not options.cluster_umis:
            for umi in bundle:
                outfile.write(bundle[umi]["read"])

        else:

            if options.stats:
                # collect pre-dudep stats 
                stats_pre_df_dict['UMI'].extend(bundle)
                stats_pre_df_dict['counts'].extend(
                    [bundle[UMI]['count'] for UMI in bundle])

                # dedup using umis and write out deduped bam
                reads, umis, umi_counts = cluster_and_reduce_with_stats(bundle, 1)
                for read in reads:
                    outfile.write(read)
                    if options.paired:
                    # TS - write out paired end mate
                        outfile.write(infile.mate(read))

                # collect post-dudep stats        
                post_cluster_umis = [x.qname.split("_")[-1] for x in reads]
                stats_post_df_dict['UMI'].extend(umis)
                stats_post_df_dict['counts'].extend(umi_counts)
                    
                average_distance = get_average_umi_distance(post_cluster_umis)
                post_cluster_stats.append(average_distance)
                
                cluster_size = len(post_cluster_umis)
                random_umis = read_gn.getUmis(cluster_size)
                average_distance_null = get_average_umi_distance(random_umis)
                post_cluster_stats_null.append(average_distance_null)

            else:
                # dedup on umi without collecting stats
                for read in cluster_and_reduce(bundle, options.threshold):
                    outfile.write(read)
                    if options.paired:
                    # TS - write out paired end mate
                        outfile.write(infile.mate(read))

    if options.stats:
        
        stats_pre_df = pd.DataFrame(stats_pre_df_dict)
        stats_post_df = pd.DataFrame(stats_post_df_dict)

        # generate histograms of counts per UMI at each position        
        UMI_counts_df_pre = pd.DataFrame(stats_pre_df.pivot_table(
            columns=stats_pre_df["counts"], values="counts", aggfunc=len))
        UMI_counts_df_post = pd.DataFrame(stats_post_df.pivot_table(
            columns=stats_post_df["counts"], values="counts", aggfunc=len))

        UMI_counts_df_pre.columns=["instances"]
        UMI_counts_df_post.columns=["instances"]

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
        cluster_bins = range(-1, int(max(pre_cluster_stats))+1)

        def bin_clusters(cluster_list, bins=cluster_bins):
            return np.digitize(cluster_list, bins, right=True)

        pre_cluster_binned = bin_clusters(pre_cluster_stats)
        post_cluster_binned = bin_clusters(post_cluster_stats)
        pre_cluster_null_binned = bin_clusters(pre_cluster_stats_null)
        post_cluster_null_binned = bin_clusters(post_cluster_stats_null)

        # tally counts across bins
        pre_cluster_tally = np.bincount(pre_cluster_binned)
        post_cluster_tally = np.bincount(post_cluster_binned)
        pre_cluster_null_tally = np.bincount(pre_cluster_null_binned)
        post_cluster_null_tally = np.bincount(post_cluster_null_binned)


        edit_distance_df = pd.DataFrame({"pre": pre_cluster_tally,
                                         "pre_null": pre_cluster_null_tally,
                                         "post": post_cluster_tally,
                                         "post_null": post_cluster_null_tally,
                                         "edit_distance": cluster_bins})
        # TS - set lowest bin (-1) to "Single_UMI"
        edit_distance_df['edit_distance'][0] = "Single_UMI"

        edit_distance_df.to_csv(options.stats + "_edit_distance.tsv",
                                index=False, sep="\t")

    # write footer and output benchmark information.
    E.info("Number of reads in: %i, Number of reads out: %i" %
           (nInput, nOutput))
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
