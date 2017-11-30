'''
extract.py - Extract UMI from fastq
====================================================

:Author: Ian Sudbery, Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Python UMI

Purpose
-------

Extract UMI barcode from a read and add it to the read name, leaving
any sample barcode in place. Can deal with paired end reads and UMIs
split across the paired ends. Can also optionally extract cell
barcodes and append these to the read name also. See the section below
for an explanation for how to encode the barcode pattern(s) to
specficy the position of the UMI +/- cell barcode.


Filtering and correcting cell barcodes
--------------------------------------

umi_tools extract can optionally filter cell barcodes
(--filter-cell-barcode) against a user-supplied whitelist
(--whitelist). If a whitelist is not available for your data, e.g
if you have performed droplet-based scRNA-Seq, you can use the
whitelist tool.

Cell barcodes which do not match the whitelist (user-generated or
automatically generated) can also be optionally corrected using the
--error-correct-cell option.

The whitelist should be in the following format (tab-separated):

    AAAAAA	AGAAAA
    AAAATC
    AAACAT
    AAACTA	AAACTN,GAACTA
    AAATAC
    AAATCA	GAATCA
    AAATGT	AAAGGT,CAATGT

Where column 1 is the whitelisted cell barcodes and column 2 is
the list (comma-separated) of other cell barcodes which should be
corrected to the barcode in column 1. If the --error-correct-cell
option is not used, this column will be ignored. Any additional columns
in the whitelist input, such as the counts columns from the output of
umi_tools whitelist, will be ignored.


Usage:
------

For single ended reads:
        umi_tools extract --extract-method=string
        --bc-pattern=[PATTERN] -L extract.log [OPTIONS]

reads from stdin and outputs to stdout.

For paired end reads:
        umi_tools extract --extract-method=string
        --bc-pattern=[PATTERN] --bc-pattern2=[PATTERN]
        --read2-in=[FASTQIN] --read2-out=[FASTQOUT] -L extract.log [OPTIONS]

reads end one from stdin and end two from FASTQIN and outputs end one to stdin
and end two to FASTQOUT.


Using regex and filtering against a whitelist of cell barcodes:
        umi_tools extract --extract-method=regex --filter-cell-barcode
        --bc-pattern=[REGEX] --whitlist=[WHITELIST_TSV]
        -L extract.log [OPTIONS]

'''
from __future__ import absolute_import
import sys
import regex
import collections
import optparse

# python 3 doesn't require izip
try:
    # Python 2
    from itertools import izip
except ImportError:
    # Python 3
    izip = zip

import umi_tools.Utilities as U
import umi_tools.umi_methods as umi_methods

# add the generic docstring text
__doc__ = __doc__ + U.GENERIC_DOCSTRING_WE


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = U.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    # (Experimental option) Retain the UMI in the sequence read"
    parser.add_option("--retain-umi", dest="retain_umi", action="store_true",
                      help=optparse.SUPPRESS_HELP)
    parser.add_option("-p", "--bc-pattern", dest="pattern", type="string",
                      help="Barcode pattern")
    parser.add_option("--bc-pattern2", dest="pattern2", type="string",
                      help="Barcode pattern for paired reads")
    parser.add_option("--3prime", dest="prime3", action="store_true",
                      help="barcode is on 3' end of read.")
    parser.add_option("--read2-in", dest="read2_in", type="string",
                      help="file name for read pairs")
    parser.add_option("--read2-out", dest="read2_out", type="string",
                      help="file to output processed paired read to")
    parser.add_option("--read2-stdout", dest="read2_stdout",
                      action="store_true",
                      help="Paired reads, send read2 to stdout, discarding read1")
    parser.add_option("--quality-filter-threshold",
                      dest="quality_filter_threshold", type="int",
                      help=("Remove reads where any UMI base quality score "
                            "falls below this threshold"))
    parser.add_option("--quality-filter-mask",
                      dest="quality_filter_mask", type="int",
                      help=("If a UMI base has a quality below this threshold, "
                            "replace the base with 'N'"))
    parser.add_option("--quality-encoding",
                      dest="quality_encoding", type="choice",
                      choices=["phred33", "phred64", "solexa"],
                      help=("Quality score encoding. Choose from 'phred33'"
                            "[33-77] 'phred64' [64-106] or 'solexa' [59-106]"))
    parser.add_option("--extract-method",
                      dest="extract_method", type="choice",
                      choices=["string", "regex"],
                      help=("How to extract the umi +/- cell barcodes, Choose "
                            "from 'string' or 'regex'"))
    parser.add_option("--filter-cell-barcode",
                      dest="filter_cell_barcode",
                      action="store_true",
                      help="Filter the cell barcodes")
    parser.add_option("--error-correct-cell",
                      dest="error_correct_cell",
                      action="store_true",
                      help=("Correct errors in the cell barcode"))
    parser.add_option("--whitelist",
                      dest="whitelist", type="string",
                      help=("A whitelist of accepted cell barcodes"))
    parser.add_option("--blacklist",
                      dest="blacklist", type="string",
                      help=("A blacklist of accepted cell barcodes"))
    parser.add_option("--reads-subset",
                      dest="reads_subset", type="int",
                      help=("Only extract from the first N reads. If N is "
                            "greater than the number of reads, all reads will "
                            "be used"))
    parser.add_option("--reconcile-pairs",
                      dest="reconcile", action="store_true",
                      help=("Allow the presences of reads in read2 input that are"
                            "not present in read1 input. This allows cell barcode"
                            "filtering of read1s without considering read2s"))
    parser.set_defaults(extract_method="string",
                        filter_cell_barcodes=False,
                        whitelist=None,
                        blacklist=None,
                        error_correct_cell=False,
                        pattern=None,
                        pattern2=None,
                        read2_in=None,
                        read2_out=False,
                        read2_stdout=False,
                        quality_filter_threshold=None,
                        quality_encoding=None,
                        reconcile=False)

    # add common options (-h/--help, ...) and parse command line

    (options, args) = U.Start(parser, argv=argv,
                              add_group_dedup_options=False,
                              add_sam_options=False)

    if options.retain_umi and not options.extract_method == "regex":
        raise ValueError("option --retain-umi only works with --extract-method=regex")

    if options.quality_filter_threshold or options.quality_filter_mask:
        if not options.quality_encoding:
            U.error("must provide a quality encoding (--quality-"
                    "encoding) to filter UMIs by quality (--quality"
                    "-filter-threshold) or mask low quality bases "
                    "with (--quality-filter-mask)")

    if not options.pattern and not options.pattern2:
        if not options.read2_in:
            U.error("Must supply --bc-pattern for single-end")
        else:
            U.error("Must supply --bc-pattern and/or --bc-pattern "
                    "if paired-end ")

    if options.pattern2:
        if not options.read2_in:
            U.error("must specify a paired fastq ``--read2-in``")

        if not options.pattern2:
            options.pattern2 = options.pattern

    extract_cell = False
    extract_umi = False

    # If the pattern is a regex we can compile the regex(es) prior to
    # ExtractFilterAndUpdate instantiation
    if options.extract_method == "regex":
        if options.pattern:
            try:
                options.pattern = regex.compile(options.pattern)
            except regex.error:
                U.error("barcode_regex '%s' is not a "
                        "valid regex" % options.pattern)

        if options.pattern2:
            try:
                options.pattern2 = regex.compile(options.pattern2)
            except regex.Error:
                U.error("barcode_regex2 '%s' is not a "
                        "valid regex" % options.pattern2)

    # check whether the regex contains a umi group(s) and cell groups(s)
    if options.extract_method == "regex":
        if options.pattern:
            for group in options.pattern.groupindex:
                if group.startswith("cell_"):
                    extract_cell = True
                elif group.startswith("umi_"):
                    extract_umi = True
        if options.pattern2:
            for group in options.pattern2.groupindex:
                if group.startswith("cell_"):
                    extract_cell = True
                elif group.startswith("umi_"):
                    extract_umi = True

    # check whether the pattern string contains umi/cell bases
    elif options.extract_method == "string":
        if options.pattern:
            if "C" in options.pattern:
                extract_cell = True
            if "N" in options.pattern:
                extract_umi = True
        if options.pattern2:
            if "C" in options.pattern2:
                extract_cell = True
            if "N" in options.pattern2:
                extract_umi = True

    if not extract_umi:
        if options.extract_method == "string":
            U.error("barcode pattern(s) do not include any umi bases "
                    "(marked with 'Ns') %s, %s" % (
                        options.pattern, options.pattern2))
        elif options.extract_method == "regex":
            U.error("barcode regex(es) do not include any umi groups "
                    "(starting with 'umi_') %s, %s" (
                        options.pattern, options.pattern2))

    if options.filter_cell_barcodes:

        if not options.whitelist:
                U.error("must provide a whitelist (--whitelist) if using "
                        "--filter-cell-barcode option")

        if not extract_cell:
            if options.extract_method == "string":
                U.error("barcode pattern(s) do not include any cell bases "
                        "(marked with 'Cs') %s, %s" % (
                            options.pattern, options.pattern2))
            elif options.extract_method == "regex":
                U.error("barcode regex(es) do not include any cell groups "
                        "(starting with 'cell_') %s, %s" (
                            options.pattern, options.pattern2))

    read1s = umi_methods.fastqIterate(options.stdin)

    # set up read extractor
    ReadExtractor = umi_methods.ExtractFilterAndUpdate(
        options.extract_method,
        options.pattern,
        options.pattern2,
        options.prime3,
        extract_cell,
        options.quality_encoding,
        options.quality_filter_threshold,
        options.quality_filter_mask,
        options.filter_cell_barcode,
        options.retain_umi)

    if options.filter_cell_barcode:
        cell_whitelist, false_to_true_map = umi_methods.getUserDefinedBarcodes(
            options.whitelist, options.error_correct_cell)

        ReadExtractor.cell_whitelist = cell_whitelist
        ReadExtractor.false_to_true_map = false_to_true_map

    if options.blacklist:
        blacklist = set()
        with U.openFile(options.blacklist, "r") as inf:
            for line in inf:
                blacklist.add(line.strip().split("\t")[0])
        ReadExtractor.cell_blacklist = blacklist

    # variables for progress monitor
    progCount = 0
    displayMax = 100000
    U.info("Starting barcode extraction")

    if options.read2_in is None:
        for read in read1s:

            # incrementing count for monitoring progress
            progCount += 1

            # Update display in every 100kth iteration
            if progCount % displayMax == 0:
                U.info("Parsed {} reads".format(progCount))

            new_read = ReadExtractor(read)

            if options.reads_subset:
                if (ReadExtractor.read_counts['Input Reads'] >
                    options.reads_subset):
                    break

            if not new_read:
                continue

            options.stdout.write(str(new_read) + "\n")

    else:
        read2s = umi_methods.fastqIterate(U.openFile(options.read2_in))

        if options.read2_out:
            read2_out = U.openFile(options.read2_out, "w")

        if options.reconcile:
            strict = False
        else:
            strict = True

        for read1, read2 in umi_methods.joinedFastqIterate(
                read1s, read2s, strict):

            # incrementing count for monitoring progress
            progCount += 1

            # Update display in every 100kth iteration
            if progCount % displayMax == 0:
                U.info("Parsed {} reads".format(progCount))
                sys.stdout.flush()

            reads = ReadExtractor(read1, read2)

            if options.reads_subset:
                if (ReadExtractor.read_counts['Input Reads'] >
                    options.reads_subset):
                    break

            if not reads:
                continue
            else:
                new_read1, new_read2 = reads

            if options.read2_stdout:
                options.stdout.write(str(new_read2) + "\n")
            else:
                options.stdout.write(str(new_read1) + "\n")

                if options.read2_out:
                    read2_out.write(str(new_read2) + "\n")

    if options.read2_out:
        read2_out.close()

    for k, v in ReadExtractor.getReadCounts().most_common():
        U.info("%s: %s" % (k, v))

    U.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
