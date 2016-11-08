'''
group.py - Group reads based on their UMI
=========================================

:Author: Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Python UMI

Purpose
-------

The purpose of this command is to identify groups of reads based on
their genomic coordinate and UMI. It is assumed that the FASTQ files
were processed with extract_umi.py before mapping and thus the UMI is
the last word of the read name. e.g:

@HISEQ:87:00000000_AATT

where AATT is the UMI sequeuence.

By default, reads are considered identical if they have the same start
coordinate, are on the same strand, and have the same UMI. Optionally,
splicing status can be considered (see below).

The start postion of a read is considered to be the start of its alignment
minus any soft clipped bases. A read aligned at position 500 with
cigar 2S98M will be assumed to start at postion 498.

Methods
-------

group can be run with multiple methods to identify group of reads with
the same (or similar) UMI(s). All methods start by identifying the
reads with the same mapping position.

The simpliest methods, unique and percentile, group reads with
the exact same UMI. The network-based methods, cluster, adjacency and
directional, build networks where nodes are UMIs and edges connect UMIs
with an edit distance <= threshold (usually 1). The groups of reads
are then defined from the network in a method-specific manner.

  "unique"
      Reads group share the exact same UMI

  "percentile"
      Reads group share the exact same UMI. UMIs with counts < 1% of the
      median counts for UMIs at the same position are ignored.

  "cluster"
      Identify clusters of connected UMIs (based on hamming distance
      threshold). Each network is a read group

  "adjacency"
      Cluster UMIs as above. For each cluster, select and remove the
      node(UMI) with the highest counts and remove all connected
      nodes. Repeat with remaining nodes until no nodes remain. Each
      removal step defines a read group.

  "directional"
      Identify clusters of connected UMIs (based on hamming distance
      threshold) and umi A counts >= (2* umi B counts) - 1. Each
      network is a read group.


Options
-------

--method (choice, string)
      Method used to identify PCR duplicates within reads. All methods
      start by identifying the reads with the same mapping position

      Options are:

      - "unique"

      - "percentile"

      - "cluster"

      - "adjacency"

      - "directional-adjacency"

--edit-distance-threshold (int)
       For the adjacency and cluster methods the threshold for the
       edit distance to connect two UMIs in the network can be
       increased. The default value of 1 works best unless the UMI is
       very long (>14bp)

--paired
       BAM is paired end - output both read pairs. This will also
       force the use of the template length to determine reads with
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

--subset
      Only consider a fraction of the reads, chosen at random. This is useful
      for doing saturation analyses.

--chrom
      Only consider a single chromosome. This is useful for debugging purposes

-i, --in-sam/-o, --out-sam
      By default, inputs are assumed to be in BAM format and output are output
      in BAM format. Use these options to specify the use of SAM format for
      inputs or outputs.

-I    input file name. The input file must be sorted and indexed.

-S    output file name.

-L    log file name.


Usage
-----

    python group -I infile.bam -S grouped.bam -L group.log


.. note::
   In order to get a valid sam/bam file you need to redirect logging
   information or turn it off logging via -v 0. You can redirect the
   logging to a file with -L <logfile> or use the --log2stderr option
   to send the logging to stderr.

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

try:
    import umi_tools.Utilities as U
    import umi_tools.network as network
    import umi_tools.umi_methods as umi_methods
except:
    import Utilities as U
    import network
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
    parser.add_option("-o", "--out-sam", dest="out_sam", action="store_true",
                      help="Output alignments in sam format [default=%default]",
                      default=False)
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
    parser.add_option("--per-contig", dest="per_contig", action="store_true",
                      default=False,
                      help=("dedup per contig,"
                            " e.g for transcriptome where contig = gene"))
    parser.add_option("--whole-contig", dest="whole_contig", action="store_true",
                      default=False,
                      help="Read whole contig before outputting bundles: guarantees that no reads"
                           "are missed, but increases memory usage")
    parser.add_option("--read-length", dest="read_length", action="store_true",
                      default=False,
                      help=("use read length in addition to position and UMI"
                            "to identify possible duplicates [default=%default]"))
    parser.add_option("--group-out", dest="tsv", type="string",
                      help="Outfile name for file mapping read id to read group",
                      default=None)
    parser.add_option("--output-bam", dest="output_bam", action="store_true",
                      default=False,
                      help=("output a bam file with read groups tagged using the UG tag"
                            "[default=%default]"))

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
        assert options.output_bam, (
            "To output a bam you must include --output-bam option")
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

    if options.output_bam:
        outfile = pysam.Samfile(out_name, out_mode, template=infile)
        if options.paired:
            outfile = umi_methods.TwoPassPairWriter(infile, outfile)
    else:
        outfile = None

    if options.tsv:
        mapping_outfile = U.openFile(options.tsv, "w")
        mapping_outfile.write(
            "read_id\tcontig\tposition\tumi\tcount\tfinal_umi\tunique_id\n")

    nInput, nOutput, unique_id = 0, 0, 0

    for bundle in umi_methods.get_bundles(
            infile,
            ignore_umi=False,
            subset=False,
            quality_threshold=0,
            paired=options.paired,
            chrom=options.chrom,
            spliced=options.spliced,
            soft_clip_threshold=options.soft,
            per_contig=options.per_contig,
            whole_contig=options.whole_contig,
            read_length=options.read_length,
            all_reads=True):

        nInput += sum([bundle[umi]["count"] for umi in bundle])

        if nOutput % 10000 == 0:
            U.debug("Outputted %i" % nOutput)

        if nInput % 1000000 == 0:
            U.debug("Read %i input reads" % nInput)

        # set up ReadCluster functor with methods specific to
        # specified options.method
        processor = network.ReadClusterer(options.method)

        clusters, umis, umi_counts, topologies, nodes = processor(
            bundle=bundle,
            threshold=options.threshold,
            stats=True,
            deduplicate=False)

        for ix, umi_group in enumerate(clusters):
            for umi in umi_group:
                reads = bundle[umi]['read']
                for read in reads:
                    if outfile:
                        # Add the 'UG' tag to the read
                        read.tags += [('UG', unique_id)]
                        read.tags += [('FU', umis[ix])]

                        if options.paired:
                            # if paired, we need to supply the tags to
                            # add to the paired read
                            outfile.write(read, add_tag=True,
                                          unique_id=unique_id, umi=umis[ix])
                        else:
                            outfile.write(read)

                    if options.tsv:
                        mapping_outfile.write("%s\n" % "\t".join(map(str, (
                            read.query_name, read.reference_name,
                            umi_methods.get_read_position(read, options.soft)[1],
                            umi_methods.get_umi(read),
                            umi_counts[ix], umis[ix], unique_id))))

                    nOutput += 1

            unique_id += 1

    if outfile:
        outfile.close()

    if options.tsv:
        mapping_outfile.close()

    # write footer and output benchmark information.
    U.info("Number of reads in: %i, Number of reads out: %i, Number of groups: %i" %
           (nInput, nOutput, unique_id))
    U.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
