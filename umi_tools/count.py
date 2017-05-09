'''
count.py - Count reads per gene using UMIs
==========================================

:Author: Ian Sudbery, Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Python UMI

Purpose
-------

The purpose of this command is to counts the number of reads per gene based
on the mapping co-ordinate and the UMI attached to the read.
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

By default, reads are considered identical (and therefore only counted
once) if they have the same start coordinate, are on the same strand,
and have the same UMI. Optionally, splicing status can be considered
(see below).

The start postion of a read is considered to be the start of its alignment
minus any soft clipped bases. A read aligned at position 500 with
cigar 2S98M will be assumed to start at postion 498.

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

--read-length
      Use the read length as as a criteria when deduping, for e.g sRNA-Seq

--subset (float, [0-1])
      Only consider a fraction of the reads, chosen at random. This is useful
      for doing saturation analyses.

--chrom (string)
      Only consider a single chromosome. This is useful for debugging purposes

--per-contig (string)
      Count per contig (field 3 in BAM; RNAME).
      All reads with the same contig will be
      considered to have the same alignment position. This is useful
      if your library prep generates PCR duplicates with non identical
      alignment positions such as CEL-Seq. In this case, you could
      align to a reference transcriptome with one transcript per gene

--per-gene (string)
      Count per gene. As above except with this option you can
      align to a reference transcriptome with more than one transcript
      per gene. You need to also provide --gene-transcript-map option

--gene-transcript-map (string)
      File mapping genes to transripts (tab separated), e.g:

      gene1   transcript1
      gene1   transcript2
      gene2   transcript3

--gene-tag (string)
      Count per gene. As per --per-gene except here the gene
      information is encoded in the bam read tag specified so you do
      not need to supply the --gene-transcript-map

--skip-tags-regex (string)
      Used in conjunction with the --gene-tag option. Skip any reads
      where the gene tag matches this regex.
      Defualt matches anything which starts with "__" or "Unassigned":
      ("^[__|Unassigned]")

-i, --in-sam
      By default, inputs are assumed to be in BAM format.
      Use this option to specify the use of SAM format.

-I    (string, filename) input file name
      The input file must be sorted and indexed.

-S    (string, filename) output file name

-L    (string, filename) log file name


Usage
-----

    python count -I infile.bam -S gene_counts.tsv -L dedup.log


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
    parser.add_option("--umi-tag", dest="umi_tag",
                      type="string", help="tag containing umi",
                      default='RX')
    parser.add_option("--extract-umi-method", dest="get_umi_method", type="choice",
                      choices=("read_id", "tag"), default="read_id",
                      help="where is the read UMI encoded? [default=%default]")
    parser.add_option("--subset", dest="subset", type="float",
                      help="Use only a fraction of reads, specified by subset",
                      default=None)
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
    parser.add_option("--mapping-quality", dest="mapping_quality",
                      type="int",
                      help="Minimum mapping quality for a read to be retained"
                      " [default=%default]",
                      default=0)
    parser.add_option("--per-contig", dest="per_contig", action="store_true",
                      default=False,
                      help=("dedup per contig (field 3 in BAM; RNAME),"
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
    parser.add_option("--skip-tags-regex", dest="skip_regex",
                      type="string",
                      help=("Used with --gene-tag. "
                            "Ignore reads where the gene-tag matches this regex"),
                      default="^[__|Unassigned]")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = U.Start(parser, argv=argv)

    if options.random_seed:
        np.random.seed(options.random_seed)

    if options.stdin != sys.stdin:
        in_name = options.stdin.name
        options.stdin.close()
    else:
        raise ValueError("Input on standard in not currently supported")

    if options.in_sam:
        in_mode = "r"
    else:
        in_mode = "rb"

    if options.per_gene:
        if not options.gene_transcript_map and not options.gene_map:
            raise ValueError("--per-gene option requires --gene-transcript-map "
                             "or --gene-tag")
    try:
        re.compile(options.skip_regex)
    except re.error:
        raise ValueError("skip-regex '%s' is not a "
                         "valid regex" % options.skip_regex)

    infile = pysam.Samfile(in_name, in_mode)

    nInput, nOutput = 0, 0

    # set the method with which to extract umis from reads
    if options.get_umi_method == "read_id":
        umi_getter = partial(
            umi_methods.get_umi_read_id, sep=options.umi_sep)
    elif options.get_umi_method == "tag":
        umi_getter = partial(
            umi_methods.get_umi_tag, tag=options.umi_tag)
    else:
        raise ValueError("Unknown umi extraction method")

    if options.chrom:
        inreads = infile.fetch(reference=options.chrom)
    else:
        if options.per_gene and options.gene_transcript_map:
            metacontig2contig = umi_methods.getMetaContig2contig(
                infile, options.gene_transcript_map)
            metatag = "MC"
            inreads = umi_methods.metafetcher(infile, metacontig2contig, metatag)
            gene_tag = metatag

        else:
            inreads = infile.fetch()
            gene_tag = options.gene_tag

    options.stdout.write("%s\t%s\n" % ("gene", "count"))
    for gene, bundle, read_events in umi_methods.get_gene_count(
            inreads,
            subset=options.subset,
            quality_threshold=options.mapping_quality,
            paired=options.paired,
            per_contig=options.per_contig,
            gene_tag=options.gene_tag,
            skip_regex=options.skip_regex,
            umi_getter=umi_getter):

        umis = bundle.keys()
        counts = {umi: bundle[umi]["count"] for umi in umis}

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

    # output reads events and benchmark information.
    for event in read_events.most_common():
        U.info("%s: %s" % (event[0], event[1]))

    U.info("Number of reads counted: %i" % nOutput)

    U.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
