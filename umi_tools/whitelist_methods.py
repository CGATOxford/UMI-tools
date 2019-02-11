'''
whitelist_methods.py - Methods for whitelisting cell barcodes
=============================================================

'''


import collections
import matplotlib
import copy
import regex

# require to run on systems with no X11
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

import numpy as np
from scipy.stats import gaussian_kde
from scipy.signal import argrelextrema

import umi_tools.Utilities as U
from umi_tools._dedup_umi import edit_distance


def getKneeEstimate(cell_barcode_counts,
                    expect_cells=False,
                    cell_number=False,
                    plotfile_prefix=None):
    ''' estimate the number of "true" cell barcodes

    input:
         cell_barcode_counts = dict(key = barcode, value = count)
         expect_cells (optional) = define the expected number of cells
         cell_number (optional) = define number of cell barcodes to accept
         plotfile_prefix = (optional) prefix for plots

    returns:
         List of true barcodes
    '''

    # very low abundance cell barcodes are filtered out (< 0.001 *
    # the most abundant)
    threshold = 0.001 * cell_barcode_counts.most_common(1)[0][1]

    counts = sorted(cell_barcode_counts.values(), reverse=True)
    counts_thresh = [x for x in counts if x > threshold]
    log_counts = np.log10(counts_thresh)

    # guassian density with hardcoded bw
    density = gaussian_kde(log_counts, bw_method=0.1)

    xx_values = 10000  # how many x values for density plot
    xx = np.linspace(log_counts.min(), log_counts.max(), xx_values)

    local_min = None

    if cell_number:  # we have a prior hard expectation on the number of cells
        threshold = counts[cell_number]

    else:
        local_mins = argrelextrema(density(xx), np.less)[0]
        local_mins_counts = []

        for poss_local_min in local_mins[::-1]:

            passing_threshold = sum([y > np.power(10, xx[poss_local_min])
                                     for x, y in cell_barcode_counts.items()])
            local_mins_counts.append(passing_threshold)

            if not local_min:   # if we have selected a local min yet
                if expect_cells:  # we have a "soft" expectation
                    if (passing_threshold > expect_cells * 0.1 and
                        passing_threshold <= expect_cells):
                        local_min = poss_local_min

                else:  # we have no prior expectation
                    # TS: In abscence of any expectation (either hard or soft),
                    # this set of heuristic thresholds are used to decide
                    # which local minimum to select.
                    # This is very unlikely to be the best way to achieve this!
                    if (poss_local_min >= 0.2 * xx_values and
                        (log_counts.max() - xx[poss_local_min] > 0.5 or
                         xx[poss_local_min] < log_counts.max()/2)):
                        local_min = poss_local_min

        if local_min is not None:
            threshold = np.power(10, xx[local_min])

    if cell_number or local_min is not None:
        final_barcodes = set([
            x for x, y in cell_barcode_counts.items() if y > threshold])
    else:
        final_barcodes = None

    if plotfile_prefix:

        # colour-blind friendly colours - https://gist.github.com/thriveth/8560036
        CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                          '#f781bf', '#a65628', '#984ea3',
                          '#999999', '#e41a1c', '#dede00']
        user_line = mlines.Line2D(
            [], [], color=CB_color_cycle[0], ls="dashed",
            markersize=15, label='User-defined')
        selected_line = mlines.Line2D(
            [], [], color=CB_color_cycle[0], ls="dashed", markersize=15, label='Selected')
        rejected_line = mlines.Line2D(
            [], [], color=CB_color_cycle[3], ls="dashed", markersize=15, label='Rejected')

        # make density plot
        fig = plt.figure()
        fig1 = fig.add_subplot(111)
        fig1.plot(xx, density(xx), 'k')
        fig1.set_xlabel("Count per cell (log10)")
        fig1.set_ylabel("Density")

        if cell_number:
            fig1.axvline(np.log10(threshold), ls="dashed", color=CB_color_cycle[0])
            lgd = fig1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,
                              handles=[user_line],
                              title="Cell threshold")

        elif local_min is None:  # no local_min was accepted
            for pos in xx[local_mins]:
                fig1.axvline(x=pos, ls="dashed", color=CB_color_cycle[3])
            lgd = fig1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,
                              handles=[selected_line, rejected_line],
                              title="Possible thresholds")
        else:
            for pos in xx[local_mins]:
                if pos == xx[local_min]:  # selected local minima
                    fig1.axvline(x=xx[local_min], ls="dashed", color=CB_color_cycle[0])
                else:
                    fig1.axvline(x=pos, ls="dashed", color=CB_color_cycle[3])

            lgd = fig1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,
                              handles=[selected_line, rejected_line],
                              title="Possible thresholds")

        fig.savefig("%s_cell_barcode_count_density.png" % plotfile_prefix,
                    bbox_extra_artists=(lgd,), bbox_inches='tight')

        # make knee plot
        fig = plt.figure()
        fig2 = fig.add_subplot(111)
        fig2.plot(range(0, len(counts)), np.cumsum(counts), c="black")

        xmax = len(counts)
        if local_min is not None:
            # reasonable maximum x-axis value
            xmax = min(len(final_barcodes) * 5, xmax)

        fig2.set_xlim((0 - (0.01 * xmax), xmax))
        fig2.set_xlabel("Rank")
        fig2.set_ylabel("Cumulative count")

        if cell_number:
            fig2.axvline(x=cell_number, ls="dashed", color=CB_color_cycle[0])
            lgd = fig2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,
                              handles=[user_line],
                              title="Cell threshold")

        elif local_min is None:  # no local_min was accepted
            for local_mins_count in local_mins_counts:
                fig2.axvline(x=local_mins_count, ls="dashed",
                             color=CB_color_cycle[3])
            lgd = fig2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,
                              handles=[selected_line, rejected_line],
                              title="Possible thresholds")

        else:
            for local_mins_count in local_mins_counts:
                if local_mins_count == len(final_barcodes):  # selected local minima
                    fig2.axvline(x=local_mins_count, ls="dashed",
                                 color=CB_color_cycle[0])
                else:
                    fig2.axvline(x=local_mins_count, ls="dashed",
                                 color=CB_color_cycle[3])

            lgd = fig2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,
                              handles=[selected_line, rejected_line],
                              title="Possible thresholds")

        fig.savefig("%s_cell_barcode_knee.png" % plotfile_prefix,
                    bbox_extra_artists=(lgd,), bbox_inches='tight')

        if local_min is not None:
            colours_selected = [CB_color_cycle[0] for x in range(0, len(final_barcodes))]
            colours_rejected = ["black" for x in range(0, len(counts)-len(final_barcodes))]
            colours = colours_selected + colours_rejected
        else:
            colours = ["black" for x in range(0, len(counts))]

        fig = plt.figure()
        fig3 = fig.add_subplot(111)
        fig3.scatter(x=range(1, len(counts)+1), y=counts,
                     c=colours, s=10, linewidths=0)
        fig3.loglog()
        fig3.set_xlim(0, len(counts)*1.25)
        fig3.set_xlabel('Barcode index')
        fig3.set_ylabel('Count')

        if cell_number:
            fig3.axvline(x=cell_number, ls="dashed", color=CB_color_cycle[0])
            lgd = fig3.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,
                              handles=[user_line],
                              title="Cell threshold")
        elif local_min is None:  # no local_min was accepted
            for local_mins_count in local_mins_counts:
                fig3.axvline(x=local_mins_count, ls="dashed",
                             color=CB_color_cycle[3])
            lgd = fig3.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,
                              handles=[selected_line, rejected_line],
                              title="Possible thresholds")
        else:
            for local_mins_count in local_mins_counts:
                if local_mins_count == len(final_barcodes):  # selected local minima
                    fig3.axvline(x=local_mins_count, ls="dashed",
                                 color=CB_color_cycle[0])
                else:
                    fig3.axvline(x=local_mins_count, ls="dashed",
                                 color=CB_color_cycle[3])

            lgd = fig3.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,
                              handles=[selected_line, rejected_line],
                              title="Possible thresholds")

        fig.savefig("%s_cell_barcode_counts.png" % plotfile_prefix,
                    bbox_extra_artists=(lgd,), bbox_inches='tight')

        if not cell_number:
            with U.openFile("%s_cell_thresholds.tsv" % plotfile_prefix, "w") as outf:
                outf.write("count\taction\n")
                for local_mins_count in local_mins_counts:
                    if local_min and local_mins_count == len(final_barcodes):
                        threshold_type = "Selected"
                    else:
                        threshold_type = "Rejected"

                    outf.write("%s\t%s\n" % (local_mins_count, threshold_type))

    return final_barcodes


def getErrorCorrectMapping(cell_barcodes, whitelist, threshold=1):
    ''' Find the mappings between true and false cell barcodes based
    on an edit distance threshold.

    Any cell barcode within the threshold to more than one whitelist
    barcode will be excluded'''

    true_to_false = collections.defaultdict(set)

    whitelist = set([str(x).encode("utf-8") for x in whitelist])

    for cell_barcode in cell_barcodes:
        match = None
        barcode_in_bytes = str(cell_barcode).encode("utf-8")
        for white_cell in whitelist:

            if barcode_in_bytes in whitelist:  # don't check if whitelisted
                continue

            if edit_distance(barcode_in_bytes, white_cell) <= threshold:
                if match is not None:  # already matched one barcode
                    match = None  # set match back to None
                    break  # break and don't add to maps
                else:
                    match = white_cell.decode("utf-8")

        if match is not None:
            true_to_false[match].add(cell_barcode)

    return true_to_false


def getCellWhitelist(cell_barcode_counts,
                     expect_cells=False,
                     cell_number=False,
                     error_correct_threshold=0,
                     plotfile_prefix=None):

    cell_whitelist = getKneeEstimate(
        cell_barcode_counts, expect_cells, cell_number, plotfile_prefix)

    U.info("Finished - whitelist determination")

    true_to_false_map = None

    if cell_whitelist and error_correct_threshold > 0:
        U.info("Starting - finding putative error cell barcodes")
        true_to_false_map = getErrorCorrectMapping(
            cell_barcode_counts.keys(), cell_whitelist,
            error_correct_threshold)
        U.info("Finished - finding putative error cell barcodes")

    return cell_whitelist, true_to_false_map


def getUserDefinedBarcodes(whitelist_tsv, getErrorCorrection=False):
    cell_whitelist = []

    if getErrorCorrection:
        false_to_true_map = {}
    else:
        false_to_true_map = None

    with U.openFile(whitelist_tsv, "r") as inf:

        for line in inf:
            if line.startswith('#'):
                continue

            line = line.strip().split("\t")
            whitelist_barcode = line[0]
            cell_whitelist.append(whitelist_barcode)

            if getErrorCorrection:
                for error_barcode in line[1].split(","):
                    false_to_true_map[error_barcode] = whitelist_barcode

    return set(cell_whitelist), false_to_true_map


def checkError(barcode, whitelist, errors=1):
    '''
    Check for errors (substitutions, insertions, deletions) between a barcode
    and a set of whitelist barcodes.

    Returns the whitelist barcodes which match the input barcode
    allowing for errors. Returns as soon as two are identified.
    '''

    near_matches = []
    comp_regex = regex.compile("(%s){e<=%i}" % (barcode, errors))
    b_length = len(barcode)
    for whitelisted_barcode in whitelist:
        w_length = len(whitelisted_barcode)

        # Don't check against itself
        if barcode == whitelisted_barcode:
            continue

        # If difference in barcode lengths > number of allowed errors, continue
        if (max(b_length, w_length) > (min(b_length, w_length) + errors)):
            continue
        if comp_regex.match(whitelisted_barcode):
            near_matches.append(whitelisted_barcode)

            # Assuming downstream processes are the same for
            # (>1 -> Inf) near_matches this is OK
            if len(near_matches) > 1:
                return near_matches

    return near_matches


def errorDetectAboveThreshold(cell_barcode_counts,
                              cell_whitelist,
                              true_to_false_map,
                              errors=1,
                              resolution_method="discard"):

    assert resolution_method in ["discard", "correct"], (
        "resolution method must be discard or correct")

    error_counter = collections.Counter()

    new_true_to_false_map = copy.deepcopy(true_to_false_map)

    discard_cbs = set()

    cell_whitelist = list(cell_whitelist)
    cell_whitelist.sort(key=lambda x: cell_barcode_counts[x])

    for ix, cb in enumerate(cell_whitelist):

        near_misses = checkError(cb, cell_whitelist[ix+1:], errors=errors)

        if len(near_misses) > 0:
            error_counter["error_discarded_mt_1"]
            discard_cbs.add(cb)  # Will always discard CB from cell_whitelist

        if resolution_method == "correct" and len(near_misses) == 1:

            # Only correct substitutions as INDELs will also mess
            # up UMI so simple correction of CB is insufficient
            if regex.match("(%s){s<=%i}" % (cb, errors), near_misses[0]):
                # add corrected barcode to T:F map
                new_true_to_false_map[near_misses[0]].add(cb)
                error_counter["substitution_corrected"] += 1
            else:
                discard_cbs.add(cb)
                error_counter["indel_discarded"] += 1
        else:
            error_counter["error_discarded"] += 1

    if resolution_method == "correct":
        U.info("CBs above the knee corrected due to possible substitutions: %i" %
               error_counter["substitution_corrected"])
        U.info("CBs above the knee discarded due to possible INDELs: %i" %
               error_counter["indel_discarded"])
        U.info("CBs above the knee discarded due to possible errors from "
               "multiple other CBs: %i" % error_counter["error_discarded_mt_1"])
    else:
        U.info("CBs above the knee discarded due to possible errors: %i" %
               len(discard_cbs))

    cell_whitelist = set(cell_whitelist).difference(discard_cbs)

    return(cell_whitelist, new_true_to_false_map)
