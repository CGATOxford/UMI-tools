'''
whitelist_methods.py - Methods for whitelisting cell barcodes
=============================================================

'''


import itertools
import collections
import matplotlib
import copy
import regex

# require to run on systems with no X11
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

import numpy as np
import numpy.matlib as npm

from scipy.stats import gaussian_kde
from scipy.signal import argrelextrema

import umi_tools.Utilities as U
from umi_tools._dedup_umi import edit_distance
import pybktree


def getKneeEstimateDensity(cell_barcode_counts,
                           expect_cells=False,
                           cell_number=False,
                           plotfile_prefix=None):
    ''' estimate the number of "true" cell barcodes using a gaussian
    density-based method

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


def getKneeEstimateDistance(cell_barcode_counts,
                            cell_number=False,
                            plotfile_prefix=None):
    ''' estimate the number of "true" cell barcodes via a knee method
    which finds the point with maximum distance

    input:
         cell_barcode_counts = dict(key = barcode, value = count)
         cell_number (optional) = define number of cell barcodes to accept
         plotfile_prefix = (optional) prefix for plots

    returns:
         List of true barcodes
    '''

    def getKneeDistance(values):
        '''
        This function is based on
        https://stackoverflow.com/questions/2018178/finding-the-best-trade-off-point-on-a-curve

        and https://dataplatform.cloud.ibm.com/analytics/notebooks/54d79c2a-f155-40ec-93ec-ed05b58afa39/view?access_token=6d8ec910cf2a1b3901c721fcb94638563cd646fe14400fecbb76cea6aaae2fb1

        The idea is to draw a line from the first to last point on the
        cumulative counts curve and then find the point on the curve
        which is the maximum distance away from this line
        '''

        # get coordinates of all the points
        nPoints = len(values)
        allCoord = np.vstack((range(nPoints), values)).T

        # get the first point
        firstPoint = allCoord[0]
        # get vector between first and last point - this is the line
        lineVec = allCoord[-1] - allCoord[0]
        lineVecNorm = lineVec / np.sqrt(np.sum(lineVec**2))

        # find the distance from each point to the line:
        # vector between all points and first point
        vecFromFirst = allCoord - firstPoint

        # To calculate the distance to the line, we split vecFromFirst into two
        # components, one that is parallel to the line and one that is perpendicular
        # Then, we take the norm of the part that is perpendicular to the line and
        # get the distance.
        # We find the vector parallel to the line by projecting vecFromFirst onto
        # the line. The perpendicular vector is vecFromFirst - vecFromFirstParallel
        # We project vecFromFirst by taking the scalar product of the vector with
        # the unit vector that points in the direction of the line (this gives us
        # the length of the projection of vecFromFirst onto the line). If we
        # multiply the scalar product by the unit vector, we have vecFromFirstParallel

        scalarProduct = np.sum(
            vecFromFirst * npm.repmat(lineVecNorm, nPoints, 1), axis=1)
        vecFromFirstParallel = np.outer(scalarProduct, lineVecNorm)
        vecToLine = vecFromFirst - vecFromFirstParallel

        # distance to line is the norm of vecToLine
        distToLine = np.sqrt(np.sum(vecToLine ** 2, axis=1))

        # knee/elbow is the point with max distance value
        idxOfBestPoint = np.argmax(distToLine)

        return(distToLine, idxOfBestPoint)

    counts = [x[1] for x in cell_barcode_counts.most_common()]
    values = list(np.cumsum(counts))

    # We need to perform the distance knee iteratively with reduced
    # number of CBs since it's sensitive to the number of CBs input
    # and overestimates if too many CBs are used
    previous_idxOfBestPoint = 0
    distToLine, idxOfBestPoint = getKneeDistance(values)
    if idxOfBestPoint == 0:
        raise ValueError("Something's gone wrong here!!")

    max_iterations = 100
    iterations = 0
    while idxOfBestPoint - previous_idxOfBestPoint != 0:
        previous_idxOfBestPoint = idxOfBestPoint
        iterations += 1
        if iterations > max_iterations:
            break
        distToLine, idxOfBestPoint = getKneeDistance(values[:idxOfBestPoint*3])

    knee_final_barcodes = [x[0] for x in cell_barcode_counts.most_common()[
        :idxOfBestPoint+1]]

    if cell_number:
        threshold = counts[cell_number]
        final_barcodes = set([
            x for x, y in cell_barcode_counts.items() if y > threshold])
    else:
        final_barcodes = knee_final_barcodes

    if plotfile_prefix:

        # colour-blind friendly colours - https://gist.github.com/thriveth/8560036
        CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                          '#f781bf', '#a65628', '#984ea3',
                          '#999999', '#e41a1c', '#dede00']

        user_line = mlines.Line2D(
            [], [], color=CB_color_cycle[2], ls="dashed",
            markersize=15, label='User-defined')
        selected_line = mlines.Line2D(
            [], [], color=CB_color_cycle[0], ls="dashed", markersize=15, label='Knee')

        # plot of the original curve and its corresponding distances
        plt.figure(figsize=(12, 6))
        plt.plot(distToLine, label='Distance', color='r')
        plt.plot(values, label='Cumulative', color='b')
        plt.plot([idxOfBestPoint], values[idxOfBestPoint], marker='o',
                 markersize=8, color="red", label='Knee')

        if cell_number:
            plt.axvline(x=cell_number, ls="dashed",
                        color=CB_color_cycle[2], label="User-defined")

        plt.legend()
        plt.savefig("%s_cell_barcode_knee.png" % plotfile_prefix)

        colours_selected = [CB_color_cycle[0] for x in range(0, len(final_barcodes))]
        colours_rejected = ["black" for x in range(0, len(counts)-len(final_barcodes))]
        colours = colours_selected + colours_rejected

        fig = plt.figure()
        fig3 = fig.add_subplot(111)
        fig3.scatter(x=range(1, len(counts)+1), y=counts,
                     c=colours, s=10, linewidths=0)
        fig3.loglog()
        fig3.set_xlim(0, len(counts)*1.25)
        fig3.set_xlabel('Barcode index')
        fig3.set_ylabel('Count')
        fig3.axvline(x=len(knee_final_barcodes), ls="dashed", color=CB_color_cycle[0])

        if cell_number:
            fig3.axvline(x=cell_number, ls="dashed", color=CB_color_cycle[2])

            lgd = fig3.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,
                              handles=[selected_line, user_line],
                              title="User threshold")
        else:
            lgd = fig3.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,
                              handles=[selected_line],
                              title="Knee threshold")

        fig.savefig("%s_cell_barcode_counts.png" % plotfile_prefix,
                    bbox_extra_artists=(lgd,), bbox_inches='tight')

        if not cell_number:
            with U.openFile("%s_cell_thresholds.tsv" % plotfile_prefix, "w") as outf:
                outf.write("count\n")
                outf.write("%s\n" % idxOfBestPoint)

    return(final_barcodes)


# Function contributed by https://github.com/redst4r
def getErrorCorrectMapping(cell_barcodes, whitelist, threshold=1):
    ''' Find the mappings between true and false cell barcodes based
    on an edit distance threshold.
    Any cell barcode within the threshold to more than one whitelist
    barcode will be excluded'''

    true_to_false = collections.defaultdict(set)

    # Unexpected results with cythonise hamming distance so redefine in python here
    def hamming_distance(first, second):
        ''' returns the edit distance/hamming distances between
        its two arguements '''

        # We only want to define hamming distance for barcodes with the same length
        if len(first) != len(second):
            return np.inf

        dist = sum([not a == b for a, b in zip(first, second)])
        return dist

    whitelist = set([str(x) for x in whitelist])

    U.info('building bktree')
    tree2 = pybktree.BKTree(hamming_distance, whitelist)
    U.info('done building bktree')

    for cell_barcode in cell_barcodes:

        if cell_barcode in whitelist:
            # if the barcode is already whitelisted, no need to add
            continue

        # get all members of whitelist that are at distance 1
        candidates = [white_cell for
                      d, white_cell in
                      tree2.find(cell_barcode, threshold) if
                      d > 0]

        if len(candidates) == 0:
            # the cell doesnt match to any whitelisted barcode,
            # hence we have to drop it
            # (as it cannot be asscociated with any frequent barcde)
            continue

        elif len(candidates) == 1:
            white_cell_str = candidates[0]
            true_to_false[white_cell_str].add(cell_barcode)

        else:
            # more than on whitelisted candidate:
            # we drop it as its not uniquely assignable
            continue
    return true_to_false


def getCellWhitelist(cell_barcode_counts,
                     knee_method="distance",
                     expect_cells=False,
                     cell_number=False,
                     error_correct_threshold=0,
                     plotfile_prefix=None):

    if knee_method == "distance":
            cell_whitelist = getKneeEstimateDistance(
                cell_barcode_counts, cell_number, plotfile_prefix)

    elif knee_method == "density":
        cell_whitelist = getKneeEstimateDensity(
            cell_barcode_counts, expect_cells, cell_number, plotfile_prefix)

    else:
        raise ValueError("knee_method must be 'distance' or 'density'")

    U.info("Finished - whitelist determination")

    true_to_false_map = None

    if cell_whitelist and error_correct_threshold > 0:
        U.info("Starting - finding putative error cell barcodes")
        true_to_false_map = getErrorCorrectMapping(
            cell_barcode_counts.keys(), cell_whitelist,
            error_correct_threshold)
        U.info("Finished - finding putative error cell barcodes")

    return cell_whitelist, true_to_false_map


def getUserDefinedBarcodes(whitelist_tsv, whitelist_tsv2=None,
                           getErrorCorrection=False,
                           deriveErrorCorrection=False,
                           threshold=1):
    '''
    whitelist_tsv: tab-separated file with whitelisted barcodes. First
    field should be whitelist barcodes. Second field [optional] should
    be comma-separated barcodes which are to be corrected to the
    barcode in the first field.

    whitelist_tsv2: as above but for read2s
    getErrorCorrection: extract the second field in whitelist_tsv and
    return a map of non-whitelist:whitelist

    deriveErrorCorrection: return a map of non-whitelist:whitelist
    using a simple edit distance threshold
    '''

    base2errors = {"A": ["T", "C", "G", "N"],
                   "T": ["A", "C", "G", "N"],
                   "C": ["T", "A", "G", "N"],
                   "G": ["T", "C", "A", "N"]}

    whitelist = []

    if getErrorCorrection or deriveErrorCorrection:
        false_to_true_map = {}
    else:
        false_to_true_map = None

    def singleBarcodeGenerator(whitelist_tsv):
        with U.openFile(whitelist_tsv, "r") as inf:
            for line in inf:
                if line.startswith('#'):
                    continue
                line = line.strip().split("\t")
                yield(line[0])

    def pairedBarcodeGenerator(whitelist_tsv, whitelist_tsv2):

        whitelist1 = []
        whitelist2 = []

        with U.openFile(whitelist_tsv, "r") as inf:
            for line in inf:
                if line.startswith('#'):
                    continue

                line = line.strip().split("\t")
                whitelist1.append(line[0])

        with U.openFile(whitelist_tsv2, "r") as inf2:
            for line in inf2:
                if line.startswith('#'):
                    continue

                line = line.strip().split("\t")
                whitelist2.append(line[0])

        for w1, w2 in itertools.product(whitelist1, whitelist2):
            yield(w1 + w2)

    if deriveErrorCorrection:

        if whitelist_tsv2:
            whitelist_barcodes = pairedBarcodeGenerator(whitelist_tsv, whitelist_tsv2)
        else:
            whitelist_barcodes = singleBarcodeGenerator(whitelist_tsv)

        for whitelist_barcode in whitelist_barcodes:
            whitelist.append(whitelist_barcode)

            # for every possible combination of positions for error(s)
            for positions in itertools.product(
                    range(0, len(whitelist_barcode)), repeat=threshold):

                m_bases = [base2errors[whitelist_barcode[x]] for x in positions]

                # for every possible combination of errors
                for m in itertools.product(*m_bases):
                    error_barcode = list(whitelist_barcode)

                    # add errors
                    for pos, error_base in zip(positions, m):
                        error_barcode[pos] = error_base

                    error_barcode = "".join(error_barcode)

                    # if error barcode has already been seen, must be within
                    # threshold edit distance of >1 whitelisted barcodes
                    if error_barcode in false_to_true_map:
                        # don't report multiple times for the same barcode
                        if false_to_true_map[error_barcode]:
                            U.info("Error barcode %s can be assigned to more than "
                                   "one possible true barcode: %s or %s" % (
                                       error_barcode,
                                       false_to_true_map[error_barcode],
                                       whitelist_barcode))
                        false_to_true_map[error_barcode] = None
                    else:
                        false_to_true_map[error_barcode] = whitelist_barcode

    elif getErrorCorrection:
        assert not whitelist_tsv2, ("Can only extract errors from the whitelist "
                                    "if a single whitelist is given")
        with U.openFile(whitelist_tsv, "r") as inf:

            for line in inf:

                if line.startswith('#'):
                    continue

                line = line.strip().split("\t")
                whitelist_barcode = line[0]
                whitelist.append(whitelist_barcode)

                if getErrorCorrection:
                    for error_barcode in line[1].split(","):
                        false_to_true_map[error_barcode] = whitelist_barcode

    else:  # no error correction
        if whitelist_tsv2:
            whitelist_barcodes = pairedBarcodeGenerator(whitelist_tsv, whitelist_tsv2)
        else:
            whitelist_barcodes = singleBarcodeGenerator(whitelist_tsv)

        whitelist = [x for x in whitelist_barcodes]

    return set(whitelist), false_to_true_map


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
