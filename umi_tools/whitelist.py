'''
whitelist.py - Identify the true cell barcodes
====================================================

:Author: Tom Smith, Ian Sudbery
:Release: $Id$
:Date: |today|
:Tags: Python UMI

Purpose
-------

Extract cell barcodes and identify the most likely true barcodes using
the 'knee' method


Barcode extraction
------------------

There are two methods enabled to extract the umi barocode (+/- cell
barcode). For both methods, the patterns should be provided using the
--bc-pattern and --bc-pattern options. The method is specified using
the --extract-method option

-'string':
       This should be used where the barcodes are always in the same
       place in the read.

       - N = UMI position (required)
       - C = cell barcode position (optional)
       - X = sample position (optional)

       Bases with Ns and Cs will be extracted and added to the read
       name. The corresponding sequence qualities will be removed from
       the read. Bases with an X will be reattached to the read.

       E.g. If the pattern is NNNNCC,
       Then the read:
       @HISEQ:87:00000000 read1
       AAGGTTGCTGATTGGATGGGCTAG
       DA1AEBFGGCG01DFH00B1FF0B
       +
       will become:
       @HISEQ:87:00000000_TT_AAGG read1
       GCTGATTGGATGGGCTAG
       1AFGGCG01DFH00B1FF0B
       +

       where 'TT' is the cell barcode and 'AAGG' is the UMI.

-'regex'
       This method allows for more flexible barcode extraction and
       should be used where the cell barcodes are variable in
       length. Alternatively, the regex option can also be used to
       filter out reads which do not contain an expected adapter
       sequence.

       The expected groups in the regex are:

       umi_n = UMI positions, where n can be any value (required)
       cell_n = cell barcode positions, where n can be any value (optional)
       discard_n = positions to discard, where n can be any value (optional)

       UMI positions and cell barcode positions will be extrated and
       added to the read name. The corresponding sequence qualities
       will be removed from the read. Discard bases and the
       corresponding quality scores will be removed from the read. All
       bases matched by other groups or componentts of the regex will
       reattached to the read sequence

       For example, the following regex can be used to extract reads
       from the Klein et al inDrop data:

       (?P<cell_1>.{8,12})(?P<discard_1>GAGTGATTGCTTGTGACGCCTT)(?P<cell_2>.{8})(?P<umi_1>.{6})T{3}.*

       Where only reads with a 3' T-tail and GAGTGATTGCTTGTGACGCCTT in
       the correct position to yield two cell barcodes of 8-12 and 8bp
       respectively, and a 6bp UMI will be retained.

       You can also specify fuzzy matching to allow errors. For example if
       the discard group above was specified as below this would enable
       matches with up to 2 errors in the discard_1 group.

       (?P<discard_1>GAGTGATTGCTTGTGACGCCTT{s<=2})

       Note that all UMIs must be the same length for downstream
       processing with dedup, group or count commands


Identifying the true cell barcodes
----------------------------------

We use the distribution of counts per UMI to identify the
cut-off for 'true' UMIs (the 'knee'). See this blog post for a more
detailed explanation:

https://cgatoxford.wordpress.com/2017/05/18/estimating-the-number-of-true-cell-barcodes-in-single-cell-rna-seq/

You can supply the --plot-prefix option to visualise the threshold set
for true cell barcodes.

Options
-------

--3prime
       By default the barcode is assumed to be on the 5' end of the
       read, but use this option to sepecify that it is on the 3' end
       instead. This option only works with --extact-method=string
       since 3' encoding can be specified explicitly with a regex, e.g
       ".*(?P<umi_1>.{5})$"

-L (string, filename)
       Specify a log file to retain logging information and final statistics

Usage:
------

For single ended reads:
        umi_tools whitelist --bc-pattern=[PATTERN] -L extract.log
        [OPTIONS]

reads from stdin and outputs to stdout.

For paired end reads where the cell barcodes is split across the read pairs:
        umi_tools whitelist --bc-pattern=[PATTERN]
        --bc-pattern2=[PATTERN] --read2-in=[FASTQIN] -L extract.log
        [OPTIONS]

reads end one from stdin and end two from FASTQIN and outputs to stdin


Output:
-------

The whitlist is outputted as 4 tab-separated columns:

    1. whitelisted cell barcode
    2. Other cell barcode(s) (comma-separated) to correct to the
       whitelisted barcode
    3. Count for whitelisted cell barcodes
    4. Count(s) for the other cell barcode(s) (comma-separated)

example output:

    AAAAAA      AGAAAA          146	1
    AAAATC		        22
    AAACAT		        21
    AAACTA	AAACTN,GAACTA	27	1,1
    AAATAC		        72
    AAATCA	GAATCA	        37	3
    AAATGT	AAAGGT,CAATGT	41	1,1
    AAATTG	CAATTG	        36	1
    AACAAT		        18
    AACATA		        24

If --error-correct-threshold is set to 0, columns 2 and 4 will be empty.

Command line options
--------------------

'''
import sys
import regex
import collections
import pyximport
pyximport.install(build_in_temp=False)

# python 3 doesn't require izip
try:
    # Python 2
    from itertools import izip
except ImportError:
    # Python 3
    izip = zip

try:
    import umi_tools.Utilities as U
except ImportError:
    import Utilities as U

try:
    import umi_tools.network as network
except ImportError:
    import netwxoork

try:
    from umi_tools._dedup_umi import edit_distance
except:
    from _dedup_umi import edit_distance

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

    parser.add_option("-p", "--bc-pattern", dest="pattern", type="string",
                      help="Barcode pattern")
    parser.add_option("--bc-pattern2", dest="pattern2", type="string",
                      help="Barcode pattern for paired reads")
    parser.add_option("--3prime", dest="prime3", action="store_true",
                      help="barcode is on 3' end of read.")
    parser.add_option("--read2-in", dest="read2_in", type="string",
                      help="file name for read pairs")
    parser.add_option("--extract-method",
                      dest="extract_method", type="choice",
                      choices=["string", "regex"],
                      help=("How to extract the umi +/- cell barcodes, Choose "
                            "from 'string' or 'regex'"))
    parser.add_option("--plot-prefix",
                      dest="plot_prefix", type="string",
                      help=("Prefix for plots to visualise the automated "
                            "detection of the number of 'true' cell barcodes"))
    parser.add_option("--subset-reads",
                      dest="subset_reads", type="int",
                      help=("Use only the first N reads to automatically "
                            "identify the true cell barcodes. If N is greater "
                            "than the number of reads, all reads will be used"))
    parser.add_option("--error-correct-threshold",
                      dest="error_correct_threshold",
                      type="int",
                      help=("Hamming distance for correction of "
                            "barcodes to whitelist barcodes"))
    parser.set_defaults(extract_method="string",
                        filter_cell_barcodes=False,
                        whitelist_tsv=None,
                        blacklist_tsv=None,
                        error_correct_threshold=1,
                        pattern=None,
                        pattern2=None,
                        read2_in=None,
                        plot_prefix=None,
                        subset_reads=100000000)

    # add common options (-h/--help, ...) and parse command line

    (options, args) = U.Start(parser, argv=argv,
                              add_group_dedup_options=False,
                              add_sam_options=False)

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
                options.pattern2 = regex.compile(options.barcode_regex2)
            except regex.Error:
                U.error("barcode_regex2 '%s' is not a "
                        "valid regex" % options.barcode_regex2)

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
    if not extract_cell:
        if options.extract_method == "string":
            U.error("barcode pattern(s) do not include any cell bases "
                    "(marked with 'Cs') %s, %s" % (
                        options.pattern, options.pattern2))
        elif options.extract_method == "regex":
            U.error("barcode regex(es) do not include any cell groups "
                    "(starting with 'cell_') %s, %s" (
                        options.pattern, options.pattern2))

    if options.stdin == sys.stdin:
        read1s = umi_methods.fastqIterate(U.openFile(options.stdin))
    else:
        read1s = umi_methods.fastqIterate(U.openFile(options.stdin.name))

    # set up read extractor
    ReadExtractor = umi_methods.ExtractFilterAndUpdate(
        method=options.extract_method,
        pattern=options.pattern,
        pattern2=options.pattern2,
        prime3=options.prime3,
        extract_cell=extract_cell)

    cell_barcode_counts = collections.Counter()

    n_reads = 0
    n_cell_barcodes = 0

    if not options.read2_in:
        for read1 in read1s:
            n_reads += 1
            cell_barcode = ReadExtractor.getCellBarcode(read1)
            if cell_barcode:
                cell_barcode_counts[cell_barcode] += 1
                n_cell_barcodes += 1
            if options.subset_reads:
                if n_cell_barcodes > options.subset_reads:
                    break
    else:
        read2s = umi_methods.fastqIterate(U.openFile(options.read2_in))
        for read1, read2 in izip(read1s, read2s):
            n_reads += 1
            cell_barcode = ReadExtractor.getCellBarcode(read1, read2)
            if cell_barcode:
                cell_barcode_counts[cell_barcode] += 1
            if options.subset_reads:
                if n_reads > options.subset_reads:
                    break

    cell_whitelist, true_to_false_map = umi_methods.getCellWhitelist(
        cell_barcode_counts,
        options.error_correct_threshold,
        options.plot_prefix)

    # re-make the reads1s iterator
    read1s = umi_methods.fastqIterate(U.openFile(options.stdin.name))

    for barcode in sorted(list(cell_whitelist)):

        if true_to_false_map:
            corrected_barcodes = ",".join(
                sorted(true_to_false_map[barcode]))
            corrected_barcode_counts = ",".join(
                map(str, [cell_barcode_counts[x] for x
                          in sorted(true_to_false_map[barcode])]))
        else:
            corrected_barcodes, corrected_barcode_counts = "", ""

        options.stdout.write("%s\t%s\t%s\t%s\n" % (
            barcode, corrected_barcodes, cell_barcode_counts[barcode],
            corrected_barcode_counts))

    U.info("Parsed %i reads" % n_reads)
    U.info("%i reads matched the barcode pattern" % n_cell_barcodes)
    U.info("Found %i unique cell barcodes" % len(cell_barcode_counts))

    U.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
