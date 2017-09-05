'''
tabulate.py - Output grouped UMIs and the genes they belong to in tabular form
==============================================================================

:Author: Ian Sudbery, Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Python UMI

Purpose
-------

WRITE ME
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

    group = U.OptionGroup(parser, "tabulate-specific options")

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

        # for each UMI group, output group size and total read count
        for umi_group in groups:
            cell_code = cell.decode() if options.per_cell else ''
            top_umi = umi_group[0]
            umi_count = len(umi_group)
            read_count = sum(counts[umi] for umi in umi_group)
            tmpfile.write("%s\n" % "\t".join((
                gene, cell_code, '-1', '-1', top_umi, str(umi_count), str(read_count))))

    tmpfile.close()

    umi_counts_dict = collections.defaultdict(collections.Counter)
    with U.openFile(tmpfilename, mode="r") as inf:
        for line in inf:
            gene, cell, pos, end, umi, rawumis, reads = line.strip().split("\t")
            pos, end, rawumis, reads = int(pos), int(end), int(rawumis), int(reads)
            umi_counts_dict[(gene,cell,pos,end,umi)] = (rawumis, reads)

    os.unlink(tmpfilename)

    options.stdout.write("\t".join(("gene", "cell", "pos", "end", "umi", "rawumis", "reads") if  options.per_cell else
                                   ("gene", "pos", "end", "umi", "rawumis", "reads")) + "\n")
    for gene, cell, pos, end, umi in sorted(list(umi_counts_dict.keys())):
        rawumis, reads = umi_counts_dict[(gene,cell,pos,end,umi)]

        if options.per_cell:
            options.stdout.write("%s\t%s\t%d\t%d\t%s\t%d\t%d\n" % (
                gene, cell, pos, end, umi, rawumis, reads))
        else:
            options.stdout.write("%s\t%d\t%d\t%s\t%d\t%d\n" % (
                gene, pos, end, umi, rawumis, reads))

    # output reads events and benchmark information.
    for event in bundle_iterator.read_events.most_common():
        U.info("%s: %s" % (event[0], event[1]))

    U.info("Number of reads counted: %i" % nOutput)

    U.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
