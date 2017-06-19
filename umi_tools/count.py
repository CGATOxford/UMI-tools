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

Finally, if you have used umis to extract the UMI +/- cell barcode,
you can specify --extract-umi-method=umis

By default, reads are considered identical (and therefore only counted
once) if they are assigned to the same gene. Additionally, if you have
multiple cells in a single input, you can specify the --per-cell
option to count per cell

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
      How are the barcodes encoded in the read?

      Options are:

      - "read_id" (default)
            Barcodes are contained at the end of the read separated as
            specified with --umi-separator option

      - "tag"
            Barcodes contained in a tag(s), see --umi-tag/--cell-tag
            options

      - "umis"
            Barcodes were extracted using umis (https://github.com/vals/umis)

--edit-distance-threshold (int)
       For the adjacency and cluster methods the threshold for the
       edit distance to connect two UMIs in the network can be
       increased. The default value of 1 works best unless the UMI is
very long (>14bp)

--per-contig
      Count per contig (field 3 in BAM; RNAME).
      All reads with the same contig will be
      considered to have the same alignment position. This is useful
      if your library prep generates PCR duplicates with non identical
      alignment positions such as CEL-Seq. In this case, you could
      align to a reference transcriptome with one transcript per gene

--gene-transcript-map (string)
      File mapping genes to transripts (tab separated), e.g:

      gene1   transcript1
      gene1   transcript2
      gene2   transcript3


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
    parser.add_option("--cell-tag", dest="cell_tag",
                      type="string", help="tag containing cell",
                      default=None)
    parser.add_option("--extract-umi-method", dest="get_umi_method", type="choice",
                      choices=("read_id", "tag", "umis"), default="read_id",
                      help=("how is the read UMI +/ cell barcode encoded? "
                            "[default=%default]"))
    parser.add_option("--subset", dest="subset", type="float",
                      help="Use only a fraction of reads, specified by subset",
                      default=None)
    parser.add_option("--edit-distance-threshold", dest="threshold",
                      type="int",
                      default=1,
                      help="Edit distance theshold at which to join two UMIs "
                           "when grouping UMIs. [default=%default]")
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
    parser.add_option("--per-cell", dest="per_cell", action="store_true",
                      default=False,
                      help=("Deduplicate per cell,"
                            "e.g for transcriptome where contig = transcript"
                            "must also provide a transript to gene map with"
                            "--gene-transcript-map [default=%default]"))
    parser.add_option("--wide-format-cell-counts", dest="wide_format_cell_counts",
                      action="store_true",
                      default=False,
                      help=("output the cell counts in a wide format "
                            "(rows=genes, columns=cells)"))
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

    if not options.gene_transcript_map and not options.gene_tag:
        raise ValueError("need to use either --gene-transcript-map "
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
        barcode_getter = partial(
            umi_methods.get_barcode_read_id,
            cell_barcode=options.per_cell,
            sep=options.umi_sep)

    elif options.get_umi_method == "tag":
        if options.per_cell and options.cell_tag is None:
            raise ValueError("Need to supply the --cell-tag option")
        barcode_getter = partial(
            umi_methods.get_barcode_tag,
            umi_tag=options.umi_tag,
            cell_barcode=options.per_cell,
            cell_tag=options.cell_tag)

    elif options.get_umi_method == "umis":
        barcode_getter = partial(
            umi_methods.get_barcode_umis,
            cell_barcode=options.per_cell)

    else:
        raise ValueError("Unknown UMI extraction method")

    if options.chrom:
        inreads = infile.fetch(reference=options.chrom)
    else:
        if options.gene_transcript_map:
            metacontig2contig = umi_methods.getMetaContig2contig(
                infile, options.gene_transcript_map)
            metatag = "MC"
            inreads = umi_methods.metafetcher(infile, metacontig2contig, metatag)
            gene_tag = metatag

        else:
            inreads = infile.fetch()
            gene_tag = options.gene_tag

    if options.wide_format_cell_counts:
        tmp_out_name = U.getTempFilename()
        outfile = U.openFile(tmp_out_name, "w")
        final_out_name = options.stdout.name
    else:
        outfile = options.stdout

    if options.per_cell:
        outfile.write("%s\t%s\t%s\n" % ("gene", "cell", "count"))
    else:
        outfile.write("%s\t%s\n" % ("gene", "count"))

    for gene, cell, bundle, read_events in umi_methods.get_gene_count(
            inreads,
            subset=options.subset,
            quality_threshold=options.mapping_quality,
            paired=options.paired,
            per_contig=options.per_contig,
            gene_tag=options.gene_tag,
            skip_regex=options.skip_regex,
            barcode_getter=barcode_getter):

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

        if options.per_cell:
            outfile.write("%s\t%s\t%i\n" % (gene, cell.decode(), gene_count))
        else:
            outfile.write("%s\t%i\n" % (gene, gene_count))
        nOutput += gene_count

    # pivot the counts table and write out
    if options.wide_format_cell_counts:
        counts_df = pd.read_table(tmp_out_name)
        print(counts_df)
        os.unlink(tmp_out_name) # delete the tempfile

    # output reads events and benchmark information.
    for event in read_events.most_common():
        U.info("%s: %s" % (event[0], event[1]))

    U.info("Number of reads counted: %i" % nOutput)

    U.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
