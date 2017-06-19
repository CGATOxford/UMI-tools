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


Filtering and correcting cell barcodes
--------------------------------------

umi_tools extract can optionally filter cell barcodes
(--filter-cell-barcode). This can either be done against a
user-supplied whitelist (--whitelist-tsv). If a whitelist is not
supplied, the whitelist will be generated computationally via the
'knee' method. This method uses the distribution of counts per UMI to
identify the cut-off for 'true' UMIs (the 'knee'). See this blog post
for a more detailed explanation:

https://cgatoxford.wordpress.com/2017/05/18/estimating-the-number-of-true-cell-barcodes-in-single-cell-rna-seq/

You can supply the --plot-prefix option to visualise the threshold set
for true cell barcodes.

Cell barcodes which do not match the whitelist (user-generated or
automatically generated) can also be optionally corrected using the
--error-correct-cell option. All UMIs which do not match the whitelist
but are within --error-correct-threshold (default 1) of a single
whitelisted UMI will be "corrected" to this UMI.

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
        --bc-pattern=[REGEX] --whitlist-tsv=[WHITELIST_TSV]
        -L extract.log [OPTIONS]


TO DO:
Add option to filter out cell barcodes within the auto-generate
whitelist which are within the error threshold of another whitelisted
barcode with greater counts. This is to prevent in silico separation
of a cell into two appararently distinct cells whilst guarding against
merging two truly different cells.


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
    parser.add_option("--read2-out", dest="read2_out", type="string",
                      help="file to output processed paired read to")
    parser.add_option("--read2-out-only", dest="read2_out_only",
                      action="store_true",
                      help="Paired reads, only output the second read in the pair")
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
    parser.add_option("--error-correct-threshold",
                      dest="error_correct_threshold",
                      type="int",
                      help=("Hamming distance allowed for correction"))
    parser.add_option("--plot-prefix",
                      dest="plot_prefix", type="string",
                      help=("Prefix for plots to visualise the automated "
                            "detection of the number of 'true' cell barcodes"))
    parser.add_option("--output-whitelist",
                      dest="output_whitelist", type="string",
                      help=("Write out the automatically generated whitelist"))
    parser.add_option("--whitelist-tsv",
                      dest="whitelist_tsv", type="string",
                      help=("A whitelist of accepted cell barcodes"))
    parser.add_option("--blacklist-tsv",
                      dest="blacklist_tsv", type="string",
                      help=("A blacklist of accepted cell barcodes"))
    parser.add_option("--cell-barcode-subset",
                      dest="cell_barcode_subset", type="int",
                      help=("Use only the first N reads to automatically "
                            "identify the true cell barcodes. If N is greater "
                            "than the number of reads, all reads will be used"))
    parser.add_option("--reads-subset",
                      dest="reads_subset", type="int",
                      help=("Only extract from the first N reads. If N is "
                            "greater than the number of reads, all reads will "
                            "be used"))
    parser.set_defaults(extract_method="string",
                        filter_cell_barcodes=False,
                        whitelist_tsv=None,
                        blacklist_tsv=None,
                        error_correct_cell=False,
                        error_correct_threshold=1,
                        pattern=None,
                        pattern2=None,
                        read2_in=None,
                        read2_out=False,
                        read2_out_only=False,
                        quality_filter_threshold=None,
                        quality_encoding=None,
                        plot_prefix=None,
                        output_whitelist=None,
                        cell_barcode_subset=50000000)

    # add common options (-h/--help, ...) and parse command line

    (options, args) = U.Start(parser, argv=argv)

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

    if options.whitelist_tsv:
        if options.blacklist_tsv:
            U.error("Do not supply a blacklist and a whitelist. Just "
                    "remove the blacklist barcodes from the whitelist!")

    if not extract_umi:
        if options.extract_method == "string":
            U.error("barcode pattern(s) do not include any umi bases "
                    "(marked with 'Ns') %s, %s" % (
                        options.pattern, options.pattern2))
        elif options.extract_method == "regex":
            U.error("barcode regex(es) do not include any umi groups "
                    "(starting with 'umi_') %s, %s" (
                        options.pattern, options.pattern2))

    if options.stdin == sys.stdin:
        if not options.whitelist_tsv and options.filter_cell_barcode:
            U.error("cannot support reading from stdin if correcting cell barcode")
        read1s = umi_methods.fastqIterate(U.openFile(options.stdin))
    else:
        read1s = umi_methods.fastqIterate(U.openFile(options.stdin.name))

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
        options.filter_cell_barcode)

    if options.filter_cell_barcode:
        if (not options.whitelist_tsv) or options.error_correct_cell:
            cell_barcode_counts = collections.Counter()

            n_reads = 0
            if not options.read2_in:
                for read1 in read1s:
                    n_reads += 1
                    cell_barcode = ReadExtractor.getCellBarcode(read1)
                    if cell_barcode:
                        cell_barcode_counts[cell_barcode] += 1
                    if options.cell_barcode_subset:
                        if (n_reads > options.cell_barcode_subset):
                            break
            else:
                read2s = umi_methods.fastqIterate(U.openFile(options.read2_in))
                for read1, read2 in izip(read1s, read2s):
                    n_reads += 1
                    cell_barcode = ReadExtractor.getCellBarcode(read1, read2)
                    if cell_barcode:
                        cell_barcode_counts[cell_barcode] += 1
                    if options.cell_barcode_subset:
                        if (n_reads > options.cell_barcode_subset):
                            break

            if options.blacklist_tsv:
                cell_blacklist = umi_methods.getUserDefinedBarcodes(
                    options.blacklist_tsv)
                for cell in cell_blacklist:
                    del cell_barcode_counts[cell]

            if options.whitelist_tsv:
                cell_whitelist = umi_methods.getUserDefinedBarcodes(
                    options.whitelist_tsv)
                error_correct_mappings = umi_methods.getErrorCorrectMappings(
                    cell_barcode_counts.keys(), cell_whitelist,
                    options.error_correct_threshold)
            else:
                # getCellWhitelist has not been properly defined yet!
                cell_whitelist, error_correct_mappings = umi_methods.getCellWhitelist(
                    cell_barcode_counts,
                    options.error_correct_threshold,
                    options.plot_prefix)

            # re-make the reads1s iterator
            read1s = umi_methods.fastqIterate(U.openFile(options.stdin.name))

        else:
            cell_whitelist = umi_methods.getUserDefinedBarcodes(
                options.whitelist_tsv)
            error_correct_mappings = None, None

        false_to_true_map, true_to_false_map = error_correct_mappings

        if options.output_whitelist:

            with U.openFile(options.output_whitelist, "w") as outf:

                columns = ["barcode", "count", "corrected_barcodes",
                           "corrected_barcode_counts"]
                outf.write("\t".join(columns) + "\n")

                for barcode in sorted(list(cell_whitelist)):

                    if true_to_false_map:
                        corrected_barcodes = ",".join(true_to_false_map[barcode])
                        corrected_barcode_counts = ",".join(
                            map(str, [cell_barcode_counts[x] for x
                                      in true_to_false_map[barcode]]))

                    else:
                        corrected_barcodes, corrected_barcode_counts = "", ""

                    outf.write("%s\t%s\t%s\t%s\n" % (
                        barcode, cell_barcode_counts[barcode],
                        corrected_barcodes, corrected_barcode_counts))

        ReadExtractor.cell_whitelist = cell_whitelist
        ReadExtractor.false_to_true_map = false_to_true_map

    if options.read2_in is None:
        for read in read1s:
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

        for read1, read2 in izip(read1s, read2s):
            reads = ReadExtractor(read1, read2)

            if options.reads_subset:
                if (ReadExtractor.read_counts['Input Reads'] >
                    options.reads_subset):
                    break

            if not reads:
                continue
            else:
                new_read1, new_read2 = reads

            if not options.read2_out_only:
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
