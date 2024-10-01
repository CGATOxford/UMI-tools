'''==================================================
whitelist - Identify the likely true cell barcodes
==================================================

*Extract cell barcodes and identify the most likely true cell barcodes*

Usage:
------

For single ended reads, the following reads from stdin and outputs to stdout::

        umi_tools whitelist --bc-pattern=[PATTERN] -L extract.log
        [OPTIONS]

For paired end reads where the cell barcodes is split across the read
pairs, the following reads end one from stdin and end two from FASTQIN
and outputs to stdin::

        umi_tools whitelist --bc-pattern=[PATTERN]
        --bc-pattern2=[PATTERN] --read2-in=[FASTQIN] -L extract.log
        [OPTIONS]


Output:
-------

The whitelist is outputted as 4 tab-separated columns:

    1. whitelisted cell barcode
    2. Other cell barcode(s) (comma-separated) to correct to the
       whitelisted barcode
    3. Count for whitelisted cell barcodes
    4. Count(s) for the other cell barcode(s) (comma-separated)

example output::

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

If ``--error-correct-threshold`` is set to 0, columns 2 and 4 will be empty.

Identifying the true cell barcodes
----------------------------------

In the absence of the ``--set-cell-number`` option, ``whitelist``
finds the knee in the curve for the cumulative read counts per CB or
unique UMIs per CB (``--method=[reads|umis]``). This point is referred
to as the 'knee'. Previously this point was identified using the
distribution of read counts per CB or unique UMIs per CB. The old
behaviour can be activated using ``--knee-method=density``

See this blog post for a more detailed exploration of the previous method:

https://cgatoxford.wordpress.com/2017/05/18/estimating-the-number-of-true-cell-barcodes-in-single-cell-rna-seq/

Counts per cell barcode can be performed using either read or unique
UMI counts. Use ``--method=[read|umis]`` to set the counting method.

The process of selecting the "best" local minima with
``--knee-method=density`` is not completely foolproof. We recommend
users always run whitelist with the ``--plot-prefix`` option to
visualise the set of thresholds considered for defining cell
barcodes. This option will also generate a table containing the
thresholds which were rejected if you want to manually adjust the
threshold. In addition, if you expect that a local minima will not be
found, you can use the ``--allow-threshold-error`` option to allow
``whitelist`` to proceed proceed past this stage. In addition, if you
have some prior expectation on the maximum number of cells which may
have been sequenced, you can provide this using the option
``--expect-cells`` (see below).

If you don't mind if ``whitelist --knee-method=density`` cannot
identify a suitable threshold as you intend to inspect the plots and
identify the threshold manually, provide the following options:
``--allow-threshold-error``, ``--plot-prefix=[PLOT_PREFIX]``

We expect that the default distance-based knee method should be more
robust than the density-based method. However, we haven't extensively
tested this method. If you have a dataset where you believe the
density-based method is better, please share this information with us:
https://github.com/CGATOxford/UMI-tools/issues

Finally, in some datasets there may be a risk that CBs above the
selected threshold are actually errors from another CB. We can detect
potential instances of this by looking for CBs within one error
(substition, insertion or deletion) of another CB with higher
counts. One can then either take a conservate approach (remove CB with
lower counts), or a more relaxed approach (correct CB with lower
counts to CB with higher counts). Note that correction is only
possible for substitutions since insertions & deletions may also affect
the UMI so these are always discarded. See
``--ed-above-threshold=[discard/correct]`` below. Of course, the risk
with the relaxed approach is that this may erroneously merge two truly
different CBs together and create an in-silico "doublet". The end of
the log file (--log) will detail the number of reads from CBs above
the threshold which may be errors. In most cases, we expect the number
of reads to be a very small fraction of the total reads and therefore
recommend taking the conservative approach. See
https://cgatoxford.wordpress.com/2017/05/23/estimating-the-number-of-true-cell-barcodes-in-single-cell-rna-seq-part-2/
for an analysis of errors in barcodes above the knee threshold.


whitelist-specific options
--------------------------

""""""""""""
``--method``
""""""""""""
       "reads" or "umis". Use either reads or unique UMI counts per cell

""""""""""""
``--knee-method``
""""""""""""
       "distance" or "density". Two methods are available to detect
       the 'knee' in the cell barcode count distributions. "distance"
       identifies the maximum distance between the cumulative
       distribution curve and a straight line between the first and
       last points on the cumulative distribution curve. "density"
       transforms the counts per UMI into a gaussian density and then
       finds the local minima which separates "real" from "error" cell
       barcodes. The gaussian method was the only method available
       prior to UMI-tools v1.0.0. "distance" is now the default
       method.

"""""""""""""""""""""
``--set-cell-number``
"""""""""""""""""""""
        Use this option to explicity set the number of cell barcodes
        which should be accepted. Note that the exact number of cell
        barcodes in the outputted whitelist may be slightly less than
        this if there are multiple cells observed with the same
        frequency at the threshold between accepted and rejected cell
        barcodes.

""""""""""""""""""
``--expect-cells``
""""""""""""""""""
        An upper limit estimate for the number of inputted cells. The
        knee method will now select the first threshold (order
        ascendingly) which results in the number of cell barcodes
        accepted being <= EXPECTED_CELLS and > EXPECTED_CELLS *
        0.1. Note: This is not compatible with the default
        ``--knee-method=distance`` since there is always as single
        solution using this method.

"""""""""""""""""""""""""""
``--allow-threshold-error``
"""""""""""""""""""""""""""
        This is useful if you what the command to exit with just a
        warning if a suitable threshold cannot be selected

"""""""""""""""""""""""""""""
``--error-correct-threshold``
"""""""""""""""""""""""""""""
       Hamming distance for correction of barcodes to whitelist
       barcodes. This value will also be used for error detection
       above the knee if required (``--ed-above-threshold``)

"""""""""""""""""
``--plot-prefix``
"""""""""""""""""
        Use this option to indicate the prefix for the plots and table
        describing the set of thresholds considered for defining cell barcodes

""""""""""""""""""""""""""""""""""""""""""
``--ed-above-threshold=[discard|correct]``
""""""""""""""""""""""""""""""""""""""""""
        Detect CBs above the threshold which may be sequence
        errors:
            - "discard"
                  Discard all putative error CBs.
            - "correct"
                  Correct putative substituion errors in CBs above the
                  threshold. Discard putative insertions/deletions. Note that
                  correction is only possible when the CB contains only
                  substituions since insertions and deletions may cause errors
                  in the UMI sequence too

        Where a CB could be corrected to two other CBs, correction is
        not possible. In these cases, the CB will be discarded
        regardless of which option is used.

""""""""""""""""""
``--subset-reads``
""""""""""""""""""
        Use the first N reads to automatically identify the true cell
        barcodes. If N is greater than the number of reads, all reads
        will be used. Default is 100000000 (100 Million).



'''
import sys
import regex
import collections

import umi_tools.Utilities as U
import umi_tools.Documentation as Documentation
import umi_tools.umi_methods as umi_methods
import umi_tools.extract_methods as extract_methods
import umi_tools.whitelist_methods as whitelist_methods

# python 3 doesn't require izip
try:
    # Python 2
    from itertools import izip
except ImportError:
    # Python 3
    izip = zip

# add the generic docstring text
__doc__ = __doc__ + Documentation.GENERIC_DOCSTRING_WE

usage = '''
whitelist - Generates a whitelist of accepted cell barcodes

Usage:

   Single-end:
      umi_tools whitelist [OPTIONS] [-I IN_FASTQ[.gz]] [-S OUT_TSV[.gz]]

   Paired end:
      umi_tools whitelist [OPTIONS] [-I IN_FASTQ[.gz]] [-S OUT_TSV[.gz]] --read2-in=IN2_FASTQ[.gz]

   note: If -I/-S are ommited standard in and standard out are used
         for input and output.  Input/Output will be (de)compressed if a
         filename provided to -S/-I/--read2-in ends in .gz
'''


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = U.OptionParser(version="%prog version: $Id$",
                            usage=usage,
                            description=globals()["__doc__"])
    if len(argv) == 1:
        parser.print_usage()
        print ("Required options missing, see --help for more details")
        return 1

    group = U.OptionGroup(parser, "whitelist-specific options")

    group.add_option("--plot-prefix",
                     dest="plot_prefix", type="string",
                     help=("Prefix for plots to visualise the automated "
                           "detection of the number of 'true' cell barcodes"))
    group.add_option("--subset-reads",
                     dest="subset_reads", type="int",
                     help=("Use the first N reads to automatically identify "
                           "the true cell barcodes. If N is greater than the "
                           "number of reads, all reads will be used. "
                           "Default is 100,000,000"))
    group.add_option("--error-correct-threshold",
                     dest="error_correct_threshold",
                     type="int",
                     help=("Hamming distance for correction of barcodes to "
                           "whitelist barcodes. This value will also be used "
                           "for error detection above the knee if required "
                           "(--ed-above-threshold)"))
    group.add_option("--method",
                     dest="method",
                     choices=["reads", "umis"],
                     help=("Use reads or unique umi counts per cell"))
    group.add_option("--knee-method",
                     dest="knee_method",
                     choices=["distance", "density"],
                     help=("Use distance or density methods for detection of knee"))
    group.add_option("--expect-cells",
                     dest="expect_cells",
                     type="int",
                     help=("Prior expectation on the upper limit on the "
                           "number of cells sequenced"))
    group.add_option("--allow-threshold-error",
                     dest="allow_threshold_error", action="store_true",
                     help=("Don't select a threshold. Will still "
                           "output the plots if requested (--plot-prefix)"))
    group.add_option("--set-cell-number",
                     dest="cell_number",
                     type="int",
                     help=("Specify the number of cell barcodes to accept"))

    parser.add_option("--ed-above-threshold",
                      dest="ed_above_threshold", type="choice",
                      choices=["discard", "correct"],
                      help=("Detect CBs above the threshold which may be "
                            "sequence errors from another CB and either "
                            "'discard' or 'correct'. Default=None (No correction)"))
    parser.add_option_group(group)

    parser.set_defaults(method="reads",
                        knee_method="distance",
                        extract_method="string",
                        whitelist_tsv=None,
                        blacklist_tsv=None,
                        error_correct_threshold=1,
                        pattern=None,
                        pattern2=None,
                        read2_in=None,
                        plot_prefix=None,
                        subset_reads=100000000,
                        expect_cells=False,
                        allow_threshold_error=False,
                        cell_number=False,
                        ed_above_threshold=None,
                        ignore_suffix=False)

    # add common options (-h/--help, ...) and parse command line

    (options, args) = U.Start(parser, argv=argv,
                              add_extract_options=True,
                              add_group_dedup_options=False,
                              add_umi_grouping_options=False,
                              add_sam_options=False)

    if options.filtered_out and not options.extract_method == "regex":
        U.error("Reads will not be filtered unless extract method is"
                "set to regex (--extract-method=regex)")

    if options.expect_cells:
        if options.knee_method == "distance":
            U.error("Cannot use --expect-cells with 'distance' knee "
                    "method. Switch to --knee-method=density if you want to "
                    "provide an expectation for the number of "
                    "cells. Alternatively, if you know the number of cell "
                    "barcodes, use --cell-number")
        if options.cell_number:
            U.error("Cannot supply both --expect-cells and "
                    "--cell-number options")

    extract_cell, extract_umi = U.validateExtractOptions(options)

    if not extract_cell:
        if options.extract_method == "string":
            U.error("barcode pattern(s) do not include any cell bases "
                    "(marked with 'Cs') %s, %s" % (
                        options.pattern, options.pattern2))
        elif options.extract_method == "regex":
            U.error("barcode regex(es) do not include any cell groups "
                    "(starting with 'cell_') %s, %s" % (
                        options.pattern, options.pattern2))

    read1s = umi_methods.fastqIterate(options.stdin)

    # set up read extractor
    ReadExtractor = extract_methods.ExtractFilterAndUpdate(
        method=options.extract_method,
        pattern=options.pattern,
        pattern2=options.pattern2,
        prime3=options.prime3,
        extract_cell=extract_cell)

    cell_barcode_counts = collections.Counter()

    n_reads = 0
    n_cell_barcodes = 0

    # if using the umis method, need to keep a set of umis observed
    if options.method == "umis":
        cell_barcode_umis = collections.defaultdict(set)

    # variables for progress monitor
    displayMax = 100000
    U.info("Starting barcode extraction")

    if options.filtered_out:
        filtered_out = U.openFile(options.filtered_out, "w")

    if not options.read2_in:
        for read1 in read1s:

            # Update display in every 100kth iteration
            if n_reads % displayMax == 0:
                U.info("Parsed {} reads".format(n_reads))

            n_reads += 1
            barcode_values = ReadExtractor.getBarcodes(read1)
            if barcode_values is None:
                if options.filtered_out:
                    filtered_out.write(str(read1) + "\n")
                continue
            else:
                cell, umi, _, _, _, _, _ = barcode_values
                if options.method == "umis":
                    cell_barcode_umis[cell].add(umi)
                else:
                    cell_barcode_counts[cell] += 1
                n_cell_barcodes += 1

            if options.subset_reads:
                if n_cell_barcodes > options.subset_reads:
                    break
    else:

        if options.filtered_out2:
            filtered_out2 = U.openFile(options.filtered_out2, "w")

        read2s = umi_methods.fastqIterate(U.openFile(options.read2_in))
        for read1, read2 in izip(read1s, read2s):

            # Update display in every 100kth iteration
            if n_reads % displayMax == 0:
                U.info("Parsed {} reads".format(n_reads))

            n_reads += 1

            barcode_values = ReadExtractor.getBarcodes(read1, read2)
            if barcode_values is None:
                if options.filtered_out:
                    filtered_out.write(str(read1) + "\n")
                if options.filtered_out2:
                    filtered_out2.write(str(read2) + "\n")
                continue
            else:
                cell, umi, _, _, _, _, _ = barcode_values
                if options.method == "umis":
                    cell_barcode_umis[cell].add(umi)
                else:
                    cell_barcode_counts[cell] += 1
                n_cell_barcodes += 1

            if options.subset_reads:
                if n_reads > options.subset_reads:
                    break

    U.info("Starting - whitelist determination")

    if options.method == "umis":
        for cell in cell_barcode_umis:
            cell_barcode_counts[cell] = len(cell_barcode_umis[cell])

    if options.cell_number and options.cell_number > len(cell_barcode_counts):
        raise ValueError(
            "--set-cell-barcode option specifies more cell barcodes than the "
            "number of observed cell barcodes. This may be because "
            "--subset-reads was set to a value too low to capture reads from "
            "all cells. %s cell barcodes observed from %s parsed reads. "
            "Expected>= %s cell barcodes" % (
                len(cell_barcode_counts),
                options.subset_reads,
                options.cell_number))

    cell_whitelist, true_to_false_map = whitelist_methods.getCellWhitelist(
        cell_barcode_counts,
        options.knee_method,
        options.expect_cells,
        options.cell_number,
        options.error_correct_threshold,
        options.plot_prefix)

    if cell_whitelist:
        U.info("Top %s cell barcodes passed the selected threshold" %
               len(cell_whitelist))

    if options.ed_above_threshold:
        cell_whitelist, true_to_false_map = whitelist_methods.errorDetectAboveThreshold(
            cell_barcode_counts,
            cell_whitelist,
            true_to_false_map,
            errors=options.error_correct_threshold,
            resolution_method=options.ed_above_threshold)

    if cell_whitelist:
        U.info("Writing out whitelist")
        total_correct_barcodes = 0
        total_corrected_barcodes = 0
        for barcode in sorted(list(cell_whitelist)):

            total_correct_barcodes += cell_barcode_counts[barcode]

            if true_to_false_map:
                corrected_barcodes = ",".join(
                    sorted(true_to_false_map[barcode]))

                correct_barcode_counts = [cell_barcode_counts[x] for x in
                                          sorted(true_to_false_map[barcode])]
                total_corrected_barcodes += sum(correct_barcode_counts)

                corrected_barcode_counts = ",".join(
                    map(str, correct_barcode_counts))
            else:
                corrected_barcodes, corrected_barcode_counts = "", ""

            options.stdout.write("%s\t%s\t%s\t%s\n" % (
                barcode, corrected_barcodes, cell_barcode_counts[barcode],
                corrected_barcode_counts))
    else:
        msg = ("No local minima was accepted. Recommend checking the plot "
               "output and counts per local minima (requires `--plot-prefix`"
               "option) and then re-running with manually selected threshold "
               "(`--set-cell-number` option)")

        if options.allow_threshold_error:
            U.info(msg)
        else:
            U.error(msg)

    U.info("Parsed %i reads" % n_reads)
    U.info("%i reads matched the barcode pattern" % n_cell_barcodes)
    U.info("Found %i unique cell barcodes" % len(cell_barcode_counts))

    if cell_whitelist:
        U.info("Found %i total reads matching the selected cell barcodes" %
               total_correct_barcodes)
        U.info("Found %i total reads which can be error corrected to the "
               "selected cell barcodes" % total_corrected_barcodes)

    if options.filtered_out:
        filtered_out.close()
    if options.filtered_out2:
        filtered_out2.close()

    U.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
