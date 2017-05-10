'''
count_tab.py - Count reads per gene from flatfile using UMIs
=================================================================

:Author: Ian Sudbery, Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Python UMI

Purpose
-------

The purpose of this command is to count the number of reads per gene
based on gene assigned to each each and read UMI.  It is assumed that
the FASTQ files were processed with extract_umi.py before mapping and
thus the UMI is the last word of the read name. e.g:

@HISEQ:87:00000000_AATT

where AATT is the UMI sequeuence.


Methods
-------

count can be run with multiple methods to identify group of reads with
the same (or similar) UMI(s), from which a single read is
returned. All methods start by identifying the reads with the same
mapping position.

The simpliest methods, unique and percentile, group reads with
the exact same UMI. The network-based methods, cluster, adjacency and
directional, build networks where nodes are UMIs and edges connect UMIs
with an edit distance <= threshold (usually 1). The groups of reads
are then defined from the network in a method-specific manner. For all
the network-based methods, each read group is equivalent to one read
count for the gene.

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

-I    (string, filename) input file name
      The input file must be sorted and indexed.

-S    (string, filename) output file name

-L    (string, filename) log file name


Usage
-----

The input must be in the following format (tab separated), where the
first column is the read identifier and the second column is the
assigned gene. The input must be sorted by the gene identifier:

NS500668:144:H5FCJBGXY:2:22309:18356:15843_TCTAA    ENSG00000279457.3
NS500668:144:H5FCJBGXY:3:23405:3971:19716_CGATG     ENSG00000225972.1

You can perform any required file transformation and pipe the output
directly to count_tab. For example to pipe output from featureCounts
with the '-R' option you can do the following:

    awk '$2=="Assigned" {print $1"\t"$3}' my.bam.featureCounts| sort -k2 |
    umi_tools count_tab -S gene_counts.tsv -L count.log

The tab file is assumed to contain each read id once only. For paired
end reads with featureCounts you must include the "-p" option so each
read id is included once only.

'''

import sys
import collections
import re

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
    parser.add_option("--umi-separator", dest="umi_sep",
                      type="string", help="separator between read id and UMI",
                      default="_")
    parser.add_option("--edit-distance-threshold", dest="threshold",
                      type="int",
                      default=1,
                      help="Edit distance theshold at which to join two UMIs"
                           "when clustering. [default=%default]")
    parser.add_option("--method", dest="method", type="choice",
                      choices=("adjacency", "directional",
                               "percentile", "unique", "cluster"),
                      default="directional",
                      help="method to use for umi deduping [default=%default]")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = U.Start(parser, argv=argv)

    if options.random_seed:
        np.random.seed(options.random_seed)

    if options.stdin != sys.stdin:
        infile = U.openFile(options.stdin.name, "r")
        options.stdin.close()
    else:
        infile = sys.stdin

    nInput, nOutput = 0, 0

    # set the method with which to extract umis from reads
    umi_getter = partial(
        umi_methods.get_umi_read_string, sep=options.umi_sep)

    options.stdout.write("%s\t%s\n" % ("gene", "count"))

    for gene, counts in umi_methods.get_gene_count_tab(
            infile,
            umi_getter=umi_getter):

        umis = counts.keys()

        nInput += sum(counts.values())

        # set up UMIClusterer functor with methods specific to
        # specified options.method
        processor = network.UMIClusterer(options.method)

        # group the umis
        groups = processor(
            umis,
            counts,
            threshold=options.threshold)

        gene_count = len(groups)
        options.stdout.write("%s\t%i\n" % (gene, gene_count))
        nOutput += gene_count

    U.info("Number of reads counted: %i" % nOutput)

    U.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
