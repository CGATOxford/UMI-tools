'''
count.py - Count reads per gene from BAM using UMIs
===================================================

:Author: Ian Sudbery, Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Python UMI

Purpose
-------

The purpose of this command is to count the number of reads per gene based
on the mapping co-ordinate and the UMI attached to the read.

'''

import sys
import collections
import re
import os

# required to make iteritems python2 and python3 compatible
from builtins import dict

import pysam

import numpy as np

import umi_tools.Utilities as U
import umi_tools.network as network
import umi_tools.umi_methods as umi_methods

# add the generic docstring text
__doc__ = __doc__ + U.GENERIC_DOCSTRING_GDC


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

    options.per_gene = True  # hardcodes counting to per-gene only

    U.validateSamOptions(options)

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

    infile = pysam.Samfile(in_name, in_mode)

    # write out to tempfile and then sort to stdout
    tmpfilename = U.getTempFilename()
    tmpfile = U.openFile(tmpfilename, mode="w")

    nInput, nOutput, input_reads = 0, 0, 0

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
        only_count_reads=True,
        metacontig_contig=metacontig2contig)

    for bundle, key, status in bundle_iterator(inreads):
        if status == "single_read":
            continue

        gene, cell = key

        umis = bundle.keys()
        counts = {umi: bundle[umi]["count"] for umi in umis}

        nInput += sum(counts.values())

        while nInput >= input_reads + 1000000:
            input_reads += 1000000
            U.info("Parsed %i input reads" % input_reads)

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
            tmpfile.write("%s\n" % "\t".join((gene, cell.decode(), str(gene_count))))
        else:
            tmpfile.write("%s\n" % "\t".join((gene, str(gene_count))))

        nOutput += gene_count

    tmpfile.close()

    if options.per_cell:

        gene_counts_dict = {}

        with U.openFile(tmpfilename, mode="r") as inf:
            genes = set()
            cells = set()
            for line in inf:
                gene, cell, gene_count = line.strip().split("\t")
                genes.add(gene)
                cells.add(cell)

                if gene not in gene_counts_dict:
                    gene_counts_dict[gene] = {}

                gene_counts_dict[gene][cell] = gene_count

        if options.wide_format_cell_counts:  # write out in wide format

            options.stdout.write(
                "%s\t%s\n" % ("gene", "\t".join(sorted(cells))))

            for gene in sorted(genes):
                counts = []
                for cell in sorted(cells):
                    if cell in gene_counts_dict[gene]:
                        counts.append(gene_counts_dict[gene][cell])
                    else:
                        counts.append(0)
                options.stdout.write(
                    "%s\t%s\n" % (gene, "\t".join(map(str, counts))))

        else:  # write out in long format
            options.stdout.write("%s\t%s\t%s\n" % ("gene", "cell", "count"))
            for gene in sorted(genes):
                for cell in sorted(list(gene_counts_dict[gene].keys())):
                    options.stdout.write("%s\t%s\t%s\n" % (
                        gene, cell, gene_counts_dict[gene][cell]))
    else:
        options.stdout.write("%s\t%s\n" % ("gene", "count"))

        with U.openFile(tmpfilename, mode="r") as inf:
            for line in inf:
                options.stdout.write(line)

    os.unlink(tmpfilename)

    # output reads events and benchmark information.
    for event in bundle_iterator.read_events.most_common():
        U.info("%s: %s" % (event[0], event[1]))

    U.info("Number of reads counted: %i" % nOutput)

    U.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
