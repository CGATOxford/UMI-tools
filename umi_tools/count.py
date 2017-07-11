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

    group = U.OptionGroup(parser, "count-specific options")

    parser.add_option("--wide-format-cell-counts", dest="wide_format_cell_counts",
                      action="store_true",
                      default=False,
                      help=("output the cell counts in a wide format "
                            "(rows=genes, columns=cells)"))

    parser.add_option_group(group)

    # add common options (-h/--help, ...) and parse command line
    (options, args) = U.Start(parser, argv=argv, add_group_dedup_options=False)
    
    options.per_gene = True # hardcodes counting to per-gene only

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

    if options.gene_tag and options.gene_transcript_map:
        raise ValueError("need to use either --gene-transcript-map "
                         "OR --gene-tag, please do not provide both")

    if not options.gene_transcript_map and not options.gene_tag:
        raise ValueError("need to use either --gene-transcript-map "
                         "or --gene-tag")
    try:
        re.compile(options.skip_regex)
    except re.error:
        raise ValueError("skip-regex '%s' is not a "
                         "valid regex" % options.skip_regex)

    infile = pysam.Samfile(in_name, in_mode)

    # write out to tempfile and then sort to stdout
    tmpfilename = U.getTempFilename()
    tmpfile = U.openFile(tmpfilename, mode="w")

    nInput, nOutput = 0, 0

    gene_tag = options.gene_tag
    metacontig2contig = None

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

    bundle_iterator = umi_methods.get_bundles(
        options,
        all_reads=True, # need to retain all reads to tally their number
        return_read2=False,
        metacontig_contig=metacontig2contig)

    for bundle, key, status in bundle_iterator(inreads):
        if status == "single_read":
            continue

        gene, cell = key

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
            tmpfile.write("%s\n" % ":".join((gene, cell.decode(), str(gene_count))))
        else:
            tmpfile.write("%s\n" % ":".join((gene, str(gene_count))))
        nOutput += gene_count

    tmpfile.close()



    if options.per_cell:

        if options.wide_format_cell_counts:  # pivot the counts table and write out
            counts_df = pd.read_table(tmpfilename, sep=":", header=None)
            counts_df.columns = ["gene", "cell", "count"]
            counts_df = pd.pivot_table(counts_df, values='count',
                                       index='gene', columns='cell')  # pivot
            counts_df = counts_df.fillna(0).astype(int)  # replace NA with 0
            counts_df.to_csv(options.stdout, index=True, sep="\t")

        else:
            gene_counts_dict = collections.defaultdict(collections.Counter)

            options.stdout.write("%s\t%s\t%s\n" % ("gene", "cell", "count"))
            with U.openFile(tmpfilename, mode="r") as inf:
                for line in inf:
                    gene, cell, gene_count = line.strip().split(":")
                    gene_counts_dict[gene][cell] = gene_count
                for gene in sorted(list(gene_counts_dict.keys())):
                    for cell in sorted(list(gene_counts_dict[gene].keys())):
                        gene_count = gene_counts_dict[gene][cell]
                        options.stdout.write("%s\t%s\t%s\n" % (
                            gene, cell, gene_count))
    else:
        gene_counts_dict = collections.Counter()

        options.stdout.write("%s\t%s\n" % ("gene", "count"))

        with U.openFile(tmpfilename, mode="r") as inf:

            for line in inf:
                gene, gene_count = line.strip().split(":")
                gene_counts_dict[gene] = gene_count
            for gene in sorted(list(gene_counts_dict.keys())):
                gene_count = gene_counts_dict[gene]
                options.stdout.write("%s\t%s\n" % (gene, gene_count))

    os.unlink(tmpfilename)

    # output reads events and benchmark information.
    for event in bundle_iterator.read_events.most_common():
        U.info("%s: %s" % (event[0], event[1]))

    U.info("Number of reads counted: %i" % nOutput)

    U.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
