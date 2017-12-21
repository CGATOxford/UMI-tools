'''
umi_methods.py - Methods for dealing with UMIs
=========================================================

:Author: Ian Sudbery, Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Python UMI

'''

from __future__ import absolute_import
import itertools
import collections
import random
import pysam
import re
from scipy.stats import gaussian_kde
from scipy.signal import argrelextrema
import matplotlib
# require to run on systems with no X11
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
from functools import partial

# required to make iteritems python2 and python3 compatible
from future.utils import iteritems
from builtins import dict

import umi_tools.Utilities as U
from umi_tools._dedup_umi import edit_distance


RANGES = {
    'phred33': (33, 77),
    'solexa': (59, 106),
    'phred64': (64, 106),
}

###############################################################################
# The code for Record and fastqIterate are taken from CGAT.Fastq:
# https://github.com/CGATOxford/cgat/blob/master/CGAT/Fastq.py
###############################################################################


class Record:
    """A record representing a :term:`fastq` formatted record.

    Attributes
    ----------
    identifier : string
       Sequence identifier
    seq : string
       Sequence
    quals : string
       String representation of quality scores.
    format : string
       Quality score format. Can be one of ``sanger``,
       ``phred33``, ``phred64`` or ``solexa``.

    """
    def __init__(self, identifier, seq, quals, entry_format=None):
        self.identifier, self.seq, self.quals, entry_format = (
            identifier, seq, quals, entry_format)

    def guessFormat(self):
        '''return quality score format -
        might return several if ambiguous.'''

        c = [ord(x) for x in self.quals]
        mi, ma = min(c), max(c)
        r = []
        for entry_format, v in iteritems(RANGES):
            m1, m2 = v
            if mi >= m1 and ma < m2:
                r.append(entry_format)
        return r

    def __str__(self):
        return "@%s\n%s\n+\n%s" % (self.identifier, self.seq, self.quals)


def fastqIterate(infile):
    '''iterate over contents of fastq file.'''

    def convert2string(b):
        if type(b) == str:
            return b
        else:
            return b.decode("utf-8")

    while 1:
        line1 = convert2string(infile.readline())
        if not line1:
            break
        if not line1.startswith('@'):
            U.error("parsing error: expected '@' in line %s" % line1)
        line2 = convert2string(infile.readline())
        line3 = convert2string(infile.readline())
        if not line3.startswith('+'):
            U.error("parsing error: expected '+' in line %s" % line3)
        line4 = convert2string(infile.readline())
        # incomplete entry
        if not line4:
            U.error("incomplete entry for %s" % line1)

        yield Record(line1[1:-1], line2[:-1], line4[:-1])


def joinedFastqIterate(fastq_iterator1, fastq_iterator2, strict=True):
    '''This will return an iterator that returns tuples of fastq records.
    At each step it will confirm that the first field of the read name
    (before the first whitespace character) is identical between the
    two reads. The response if it is not depends on the value of
    :param:`strict`. If strict is true an error is returned. If strict
    is `False` the second file is advanced until a read that matches
    is found.

    This allows for protocols where read one contains cell barcodes, and these
    reads have been filtered and corrected before processing without regard
    to read2

    '''

    for read1 in fastq_iterator1:
        read2 = next(fastq_iterator2)
        pair_id = read1.identifier.split()[0]
        if not strict:
            while read2.identifier.split()[0] != pair_id:
                read2 = next(fastq_iterator2)
        if not read2.identifier.split()[0] == pair_id:
            raise ValueError("\nRead pairs do not match\n%s != %s" %
                             (pair_id, read2.identifier.split()[0]))
        yield (read1, read2)
# End of FastqIterate()
###############################################################################
###############################################################################


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

    if cell_whitelist is None:
        U.error("No local minima was accepted. Recommend checking the plot "
                "output and counts per local minima "
                "(requires `--plot-prefix` option) and then re-running with "
                "manually selected threshold (`--set-cell-number` option)")

    if error_correct_threshold > 0:
        U.info("Starting - finding putative error cell barcodes")
        true_to_false_map = getErrorCorrectMapping(
            cell_barcode_counts.keys(), cell_whitelist,
            error_correct_threshold)
        U.info("Finished - finding putative error cell barcodes")
    else:
        true_to_false_map = None

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


def get_barcode_read_id(read, cell_barcode=False, sep="_"):
    ''' extract the umi +/- cell barcode from the read id using the
    specified separator '''

    try:
        if cell_barcode:
            umi = read.qname.split(sep)[-1].encode('utf-8')
            cell = read.qname.split(sep)[-2].encode('utf-8')
        else:
            umi = read.qname.split(sep)[-1].encode('utf-8')
            cell = None

        return umi, cell

    except:
        raise ValueError(
            "Could not extract UMI +/- cell barcode from the read"
            "ID, please check UMI is encoded in the read name")


def get_barcode_tag(read,
                    cell_barcode=False,
                    umi_tag='RX',
                    cell_tag=None,
                    umi_tag_split=None,
                    umi_tag_delim=None,
                    cell_tag_split="-",
                    cell_tag_delim=None):
    ''' extract the umi +/- cell barcode from the specified tag '''

    try:
        if cell_barcode:
            umi = read.get_tag(umi_tag)
            cell = read.get_tag(cell_tag)
        else:
            umi = read.get_tag(umi_tag)
            cell = None

        if umi_tag_split:
            umi = umi.split(umi_tag_split)[0]

        if umi_tag_delim:
            umi = "".join(umi.split(umi_tag_delim))

        # 10X pipelines append a 'GEM' tag to the UMI, e.g
        # GATAGATACCTAGATA-1, hence default to split by "-
        if cell and cell_tag_split:
            cell = cell.split(cell_tag_split)[0]

        if cell and cell_tag_delim:
            cell = "".join(cell.split(cell_tag_delim))

        if cell:
            cell = cell.encode('utf-8')

        return umi.encode('utf-8'), cell

    except IndexError:
        raise ValueError("Could not extract UMI +/- cell barcode from the "
                         "read tag")


def get_barcode_umis(read, cell_barcode=False):
    ''' extract the umi +/- cell barcode from the read name where the barcodes
    were extracted using umis'''

    umi, cell = None, None
    try:
        read_name_elements = read.qname.split(":")
        for element in read_name_elements:
            if element.startswith("UMI_"):
                umi = element[4:].encode('utf-8')
            elif element.startswith("CELL_") and cell_barcode:
                cell = element[5:].encode('utf-8')

        if umi is None:
            raise ValueError()

        return umi, cell

    except:
        raise ValueError("Could not extract UMI +/- cell barcode from the "
                         "read tag")


def get_umi_read_string(read_id, sep="_"):
    ''' extract the umi from the read id (input as a string) using the
    specified separator '''

    try:
        return read_id.split(sep)[-1].encode('utf-8')
    except IndexError:
        raise ValueError(
            "Could not extract UMI from the read ID, please"
            "check UMI is encoded in the read name")


def get_average_umi_distance(umis):

    if len(umis) == 1:
        return -1

    dists = [edit_distance(x, y) for
             x, y in itertools.combinations(umis, 2)]
    return float(sum(dists))/(len(dists))


def addBarcodesToIdentifier(read, UMI, cell):
    '''extract the identifier from a read and append the UMI and
    cell barcode before the first space'''

    read_id = read.identifier.split(" ")

    if cell == "":
        read_id[0] = read_id[0] + "_" + UMI
    else:
        read_id[0] = read_id[0] + "_" + cell + "_" + UMI

    identifier = " ".join(read_id)

    return identifier


def extractSeqAndQuals(seq, quals, umi_bases, cell_bases, discard_bases,
                       retain_umi=False):
    '''Remove selected bases from seq and quals
    '''

    new_seq = ""
    new_quals = ""
    umi_quals = ""
    cell_quals = ""

    ix = 0
    for base, qual in zip(seq, quals):
        if ((ix not in discard_bases) and
            (ix not in cell_bases)):

            # if we are retaining the umi, this base is both seq and umi
            if retain_umi:
                new_quals += qual
                new_seq += base
                umi_quals += qual

            else:  # base is either seq or umi
                if ix not in umi_bases:
                    new_quals += qual
                    new_seq += base
                else:
                    umi_quals += qual

        elif ix in cell_bases:
            cell_quals += qual

        ix += 1

    return new_seq, new_quals, umi_quals, cell_quals


def ExtractBarcodes(read, match,
                    extract_umi=False,
                    extract_cell=False,
                    discard=False,
                    retain_umi=False):
    '''Extract the cell and umi barcodes using a regex.match object

    inputs:

    - read 1 and read2 = Record objects
    - match = regex.match object
    - extract_umi and extract_cell = switches to determine whether these
                                     barcodes should be extracted
    - discard = is there a region(s) of the sequence which should be
      discarded entirely?
    - retain_umi = Should UMI sequence be retained on the read sequence

    returns:

        - cell_barcode = Cell barcode string
        - cell_barcode_quals = Cell barcode quality scores
        - umi = UMI barcode string.
        - umi_quals = UMI barcode quality scores
        - new_seq = Read1 sequence after extraction
        - new_quals = Read1 qualities after extraction

    Barcodes and qualities default to empty strings where extract_cell
    or extract_umi are false.

    '''
    cell_barcode, umi, cell_barcode_quals, umi_quals, new_seq, new_quals = ("",)*6

    if not extract_cell and not extract_umi:
        U.error("must set either extract_cell and/or extract_umi to true")

    groupdict = match.groupdict()
    cell_bases = set()
    umi_bases = set()
    discard_bases = set()
    for k in sorted(list(groupdict)):
        span = match.span(k)
        if extract_cell and k.startswith("cell_"):
            cell_barcode += groupdict[k]
            cell_bases.update(range(span[0], span[1]))
        elif extract_umi and k.startswith("umi_"):
            umi += groupdict[k]
            umi_bases.update(range(span[0], span[1]))
        elif discard and k.startswith("discard_"):
            discard_bases.update(range(span[0], span[1]))

    new_seq, new_quals, umi_quals, cell_quals = extractSeqAndQuals(
        read.seq, read.quals, umi_bases, cell_bases, discard_bases, retain_umi)

    return (cell_barcode, cell_barcode_quals,
            umi, umi_quals,
            new_seq, new_quals)


def get_below_threshold(umi_quals, quality_encoding,  quality_filter_threshold):
    '''test whether the umi_quals are below the threshold'''
    umi_quals = [x - RANGES[quality_encoding][0] for x in map(ord, umi_quals)]
    below_threshold = [x < quality_filter_threshold for x in umi_quals]
    return below_threshold


def umi_below_threshold(umi_quals, quality_encoding,  quality_filter_threshold):
    ''' return true if any of the umi quals is below the threshold'''
    below_threshold = get_below_threshold(
        umi_quals, quality_encoding, quality_filter_threshold)
    return any(below_threshold)


def mask_umi(umi, umi_quals, quality_encoding,  quality_filter_threshold):
    ''' Mask all positions where quals < threshold with "N" '''
    below_threshold = get_below_threshold(
        umi_quals, quality_encoding, quality_filter_threshold)
    new_umi = ""

    for base, test in zip(umi, below_threshold):
        if test:
            new_umi += "N"
        else:
            new_umi += base

    return new_umi


class ExtractFilterAndUpdate:
    ''' A functor which extracts barcodes from a read(s), filters the
    read(s) and updates the read(s). Keeps track of events in
    read_counts Counter
    '''

    def _extract_5prime(self, sequence, read=1):
        if read == 1:
            return (sequence[:self.pattern_length],
                    sequence[self.pattern_length:])
        elif read == 2:
            return (sequence[:self.pattern_length2],
                    sequence[self.pattern_length2:])

    def _extract_3prime(self, sequence, read=1):
        if read == 1:
            return (sequence[-self.pattern_length:],
                    sequence[:-self.pattern_length])
        if read == 2:
            return (sequence[-self.pattern_length2:],
                    sequence[:-self.pattern_length2])

    def _joiner_5prime(self, sequence, sample):
        return sample + sequence

    def _joiner_3prime(self, sequence, sample):
        return sequence + sample

    def _getBarcodesString(self, read1, read2=None):

        if self.pattern:
            bc1, sequence1 = self.extract(read1.seq)
            bc_qual1, seq_qual1 = self.extract(read1.quals)
            umi_quals = [bc_qual1[x] for x in self.umi_bases]

            umi = "".join([bc1[x] for x in self.umi_bases])
            cell = "".join([bc1[x] for x in self.cell_bases])
            sample1 = "".join([bc1[x] for x in self.bc_bases])
            sample_qual1 = "".join([bc_qual1[x] for x in self.bc_bases])
            new_seq = self.joiner(sequence1, sample1)
            new_quals = self.joiner(seq_qual1, sample_qual1)

        else:
            cell, umi, umi_quals, new_seq, new_quals = ("",)*5

        if self.pattern2:
            bc2, sequence2 = self.extract(read2.seq, read=2)
            bc_qual2, seq_qual2 = self.extract(read2.quals, read=2)
            umi_quals2 = [bc_qual2[x] for x in self.umi_bases2]

            umi2 = "".join([bc2[x] for x in self.umi_bases2])
            cell2 = "".join([bc2[x] for x in self.cell_bases2])
            sample2 = "".join([bc2[x] for x in self.bc_bases2])
            sample_qual2 = "".join([bc_qual2[x] for x in self.bc_bases2])
            new_seq2 = self.joiner(sequence2, sample2)
            new_quals2 = self.joiner(seq_qual2, sample_qual2)

            cell += cell2
            umi += umi2
            umi_quals += umi_quals2
        else:
            new_seq2, new_quals2 = "", ""

        return cell, umi, umi_quals, new_seq, new_quals, new_seq2, new_quals2

    def _getBarcodesRegex(self, read1, read2=None):
        ''' '''
        # first check both regexes for paired end samples to avoid uneccessarily
        # extracting barcodes from read1 where regex2 doesn't match read2
        if self.pattern:
            match = self.pattern.match(read1.seq)
            if not match:
                self.read_counts['regex does not match read1'] += 1
                return None
            else:
                self.read_counts['regex matches read1'] += 1

        if read2 and self.pattern2:
            match2 = self.pattern2.match(read2.seq)
            if not match2:
                self.read_counts['regex does not match read1'] += 1
                return None
            else:
                self.read_counts['regex matches read2'] += 1

        # now extract barcodes
        if self.pattern:
            (cell, cell_quals,
             umi, umi_quals,
             new_seq, new_quals) = ExtractBarcodes(
                 read1, match, extract_cell=self.extract_cell,
                 extract_umi=True, discard=True, retain_umi=self.retain_umi)
        else:
            cell, cell_quals, umi, umi_quals, new_seq, new_quals = ("",)*6

        if read2 and self.pattern2:
            (cell2, cell_quals2,
             umi2, umi_quals2,
             new_seq2, new_quals2) = ExtractBarcodes(
                 read2, match2, extract_cell=self.extract_cell,
                 extract_umi=True, discard=True)

            cell += cell2
            cell_quals += cell_quals2
            umi += umi2
            umi_quals += umi_quals2
        else:
            new_seq2, new_quals2 = "", ""

        return cell, umi, umi_quals, new_seq, new_quals, new_seq2, new_quals2

    def _getCellBarcodeString(self, read1, read2=None):

        if self.pattern:
            bc1, sequence1 = self.extract(read1.seq)
            cell = "".join([bc1[x] for x in self.cell_bases])
        else:
            cell = ""

        if self.pattern2:
            bc2, sequence2 = self.extract(read2.seq, read=2)
            cell2 = "".join([bc2[x] for x in self.cell_bases2])

            cell += cell2

        return cell

    def _getCellBarcodeRegex(self, read1, read2=None):

        if read2 is None:
            match = self.pattern.match(read1.seq)
            if match:
                cell_barcode = ExtractBarcodes(
                    read1, match, extract_cell=True, extract_umi=False)[0]
                return cell_barcode
            else:
                return None

        else:

            match1, match2 = None, None

            if self.pattern:
                match1 = self.pattern.match(read1.seq)

            if self.pattern2:
                match2 = self.pattern2.match(read2.seq)

            # check matches have been made
            if not ((self.pattern and not match1) or
                    (self.pattern2 and not match2)):
                cell_barcode1, cell_barcode2 = "", ""

                if self.pattern:
                    cell_barcode1 = ExtractBarcodes(
                        read1, match1, extract_cell=True, extract_umi=False)[0]
                if self.pattern2:
                    cell_barcode2 = ExtractBarcodes(
                        read2, match2, extract_cell=True, extract_umi=False)[0]

                cell_barcode = cell_barcode1 + cell_barcode2

                return cell_barcode
            else:
                return None

    def filterQuality(self, umi_quals):
        if umi_below_threshold(
                umi_quals, self.quality_encoding,
                self.quality_filter_threshold):
            self.read_counts['filtered: umi quality'] += 1
            return True
        else:
            return False

    def maskQuality(self, umi, umi_quals):
        '''mask low quality bases and return masked umi'''
        masked_umi = mask_umi(umi, umi_quals,
                              self.quality_encoding,
                              self.quality_filter_mask)
        if masked_umi != umi:
            self.read_counts['UMI masked'] += 1
            return masked_umi
        else:
            return umi

    def filterCellBarcode(self, cell):
        '''Filter out cell barcodes not in the whitelist, with
        optional cell barcode error correction'''

        if self.cell_blacklist and cell in self.cell_blacklist:
            self.read_counts['Cell barcode in blacklist'] += 1
            return None

        if cell not in self.cell_whitelist:
            if self.false_to_true_map:
                if cell in self.false_to_true_map:
                    cell = self.false_to_true_map[cell]
                    self.read_counts['False cell barcode. Error-corrected'] += 1
                else:
                    self.read_counts['Filtered cell barcode. Not correctable'] += 1
                    return None
            else:
                self.read_counts['Filtered cell barcode'] += 1
                return None

        if self.cell_blacklist and cell in self.cell_blacklist:
            self.read_counts['Cell barcode corrected to barcode blacklist'] += 1
            return None

        return cell

    def __init__(self,
                 method="string",
                 pattern=None,
                 pattern2=None,
                 prime3=False,
                 extract_cell=False,
                 quality_encoding=None,
                 quality_filter_threshold=False,
                 quality_filter_mask=False,
                 filter_cell_barcode=False,
                 retain_umi=False):

        self.read_counts = collections.Counter()
        self.method = method
        self.pattern = pattern
        self.pattern2 = pattern2
        self.extract_cell = extract_cell
        self.quality_encoding = quality_encoding
        self.quality_filter_threshold = quality_filter_threshold
        self.quality_filter_mask = quality_filter_mask
        self.filter_cell_barcodes = filter_cell_barcode
        self.retain_umi = retain_umi

        self.cell_whitelist = None  # These will be updated if required
        self.false_to_true_map = None  # These will be updated if required
        self.cell_blacklist = None  # These will be updated if required

        # If the pattern is a string we can identify the position of
        # the cell and umi bases at instantiation
        if method == "string":
            if prime3:
                self.extract = self._extract_3prime
                self.joiner = self._joiner_3prime
            else:
                self.extract = self._extract_5prime
                self.joiner = self._joiner_5prime

            if pattern:
                self.pattern_length = len(pattern)
                self.umi_bases = [x for x in range(len(pattern)) if pattern[x] is "N"]
                self.bc_bases = [x for x in range(len(pattern)) if pattern[x] is "X"]
                self.cell_bases = [x for x in range(len(pattern)) if pattern[x] is "C"]

            if pattern2:
                self.pattern_length2 = len(pattern2)
                self.umi_bases2 = [x for x in range(len(pattern2))
                                   if pattern2[x] is "N"]
                self.bc_bases2 = [x for x in range(len(pattern2))
                                  if pattern2[x] is "X"]
                self.cell_bases2 = [x for x in range(len(pattern2))
                                    if pattern2[x] is "C"]

            self.getCellBarcode = self._getCellBarcodeString
            self.getBarcodes = self._getBarcodesString

        elif method == "regex":
            self.getCellBarcode = self._getCellBarcodeRegex
            self.getBarcodes = self._getBarcodesRegex

    def getReadCounts(self):
        return self.read_counts

    def __call__(self, read1, read2=None):

        self.read_counts['Input Reads'] += 1

        umi_values = self.getBarcodes(read1, read2)
        if umi_values is None:
            return None
        else:
            cell, umi, umi_quals, new_seq, new_quals, new_seq2, new_quals2 = umi_values

        if self.quality_filter_threshold:
            if self.filterQuality(umi_quals):
                return None

        if self.quality_filter_mask:
            umi = self.maskQuality(umi, umi_quals)

        if self.filter_cell_barcodes:
            cell = self.filterCellBarcode(cell)
            if cell is None:
                return None

        self.read_counts['Reads output'] += 1

        new_identifier = addBarcodesToIdentifier(
            read1, umi, cell)
        read1.identifier = new_identifier
        if self.pattern:  # seq and quals need to be updated
            read1.seq = new_seq
            read1.quals = new_quals

        if read2:
            new_identifier2 = addBarcodesToIdentifier(
                read2, umi, cell)
            read2.identifier = new_identifier2
            if self.pattern2:   # seq and quals need to be updated
                read2.seq = new_seq2
                read2.quals = new_quals2

        if read2 is None:
            return read1
        else:
            return read1, read2


class TwoPassPairWriter:
    '''This class makes a note of reads that need their pair outputting
    before outputting.  When the chromosome changes, the reads on that
    chromosome are read again, and any mates of reads already output
    are written and removed from the list of mates to output. When
    close is called, this is performed for the last chormosome, and
    then an algorithm identicate to pysam's mate() function is used to
    retrieve any remaining mates.

    This means that if close() is not called, at least as contigs
    worth of mates will be missing. '''

    def __init__(self, infile, outfile, tags=False):
        self.infile = infile
        self.outfile = outfile
        self.read1s = set()
        self.chrom = None

    def write(self, read, unique_id=None, umi=None, unmapped=False):
        '''Check if chromosome has changed since last time. If it has, scan
        for mates. Write the read to outfile and save the identity for paired
        end retrieval'''

        if unmapped or read.mate_is_unmapped:
            self.outfile.write(read)
            return

        if not self.chrom == read.reference_name:
            self.write_mates()
            self.chrom = read.reference_name

        key = read.query_name, read.next_reference_name, read.next_reference_start
        self.read1s.add(key)

        self.outfile.write(read)

    def write_mates(self):
        '''Scan the current chromosome for matches to any of the reads stored
        in the read1s buffer'''
        if self.chrom is not None:
            U.debug("Dumping %i mates for contig %s" % (
                len(self.read1s), self.chrom))

        for read in self.infile.fetch(reference=self.chrom, multiple_iterators=True):
            if any((read.is_unmapped, read.mate_is_unmapped, read.is_read1)):
                continue

            key = read.query_name, read.reference_name, read.reference_start
            if key in self.read1s:
                self.outfile.write(read)
                self.read1s.remove(key)

        U.debug("%i mates remaining" % len(self.read1s))

    def close(self):
        '''Write mates for remaining chromsome. Search for matches to any
        unmatched reads'''

        self.write_mates()
        U.info("Searching for mates for %i unmatched alignments" %
               len(self.read1s))

        found = 0
        for read in self.infile.fetch(until_eof=True, multiple_iterators=True):

            if read.is_unmapped:
                continue

            key = read.query_name, read.reference_name, read.reference_start
            if key in self.read1s:
                self.outfile.write(read)
                self.read1s.remove(key)
                found += 1
                continue

        U.info("%i mates never found" % len(self.read1s))
        self.outfile.close()


def getMetaContig2contig(bamfile, gene_transcript_map):
    ''' '''
    references = set(bamfile.references)
    metacontig2contig = collections.defaultdict(set)
    for line in U.openFile(gene_transcript_map, "r"):

        if line.startswith("#"):
            continue

        if len(line.strip()) == 0:
            break

        gene, transcript = line.strip().split("\t")
        if transcript in references:
            metacontig2contig[gene].add(transcript)

    return metacontig2contig


def metafetcher(bamfile, metacontig2contig, metatag):
    ''' return reads in order of metacontigs'''
    for metacontig in metacontig2contig:
        for contig in metacontig2contig[metacontig]:
            for read in bamfile.fetch(contig):
                read.tags += [(metatag, metacontig)]
                yield read


def find_splice(cigar):
    '''Takes a cigar string and finds the first splice position as
    an offset from the start. To find the 5' end (read coords) of
    the junction for a reverse read, pass in the reversed cigar tuple'''

    offset = 0
    # a soft clip at the end of the read is taken as splicing
    # where as a soft clip at the start is not.
    if cigar[0][0] == 4:
        offset = cigar[0][1]
        cigar = cigar[1:]

    for op, bases in cigar:
        if op in (3, 4):
            # N or S: found the splice
            return offset
        elif op in (0, 2, 7, 8):
            # M, D, = or X: reference consuming
            offset += bases
        elif op in (1, 5, 6):
            # I, H, P: non-reference consuming
            continue
        else:
            raise ValueError("Bad Cigar operation: %i" % op)

    return False


def get_read_position(read, soft_clip_threshold):
    ''' get the read position (taking account of clipping) '''
    is_spliced = False

    if read.is_reverse:
        pos = read.aend
        if read.cigar[-1][0] == 4:
            pos = pos + read.cigar[-1][1]
        start = read.pos

        if ('N' in read.cigarstring or
            (read.cigar[0][0] == 4 and
             read.cigar[0][1] > soft_clip_threshold)):

            cigar = read.cigar[::-1]
            is_spliced = find_splice(cigar)
    else:
        pos = read.pos
        if read.cigar[0][0] == 4:
            pos = pos - read.cigar[0][1]
        start = pos

        if ('N' in read.cigarstring or
            (read.cigar[-1][0] == 4 and
             read.cigar[-1][1] > soft_clip_threshold)):
            is_spliced = find_splice(read.cigar)

    return start, pos, is_spliced


class get_bundles:

    ''' A functor - When called returns a dictionary of dictionaries,
    representing the unique reads at a position/spliced/strand
    combination. The key to the top level dictionary is a umi. Each
    dictionary contains a "read" entry with the best read, and a count
    entry with the number of reads with that
    position/spliced/strand/umi combination

    initiation arguments:

    options: script options

    all_reads: if true, return all reads in the dictionary. Else,
    return the 'best' read (using MAPQ +/- multimapping) for each key

    return_read2: Return read2s immediately as a single read

    metacontig_contig: Maps metacontigs to the consistuent contigs
    '''

    def __init__(self,
                 options,
                 only_count_reads=False,
                 all_reads=False,
                 return_unmapped=False,
                 return_read2=False,
                 metacontig_contig=None):

        self.options = options
        self.only_count_reads = only_count_reads
        self.all_reads = all_reads
        self.return_unmapped = return_unmapped
        self.return_read2 = return_read2
        self.metacontig_contig = metacontig_contig

        self.contig_metacontig = {}
        if self.metacontig_contig:
            for metacontig in metacontig_contig:
                for contig in metacontig_contig[metacontig]:
                    self.contig_metacontig[contig] = metacontig

        # set the method with which to extract umis from reads
        if self.options.get_umi_method == "read_id":
            self.barcode_getter = partial(
                get_barcode_read_id,
                cell_barcode=self.options.per_cell,
                sep=self.options.umi_sep)

        elif self.options.get_umi_method == "tag":
            self.barcode_getter = partial(
                get_barcode_tag,
                umi_tag=self.options.umi_tag,
                cell_barcode=self.options.per_cell,
                cell_tag=self.options.cell_tag,
                umi_tag_split=self.options.umi_tag_split,
                umi_tag_delim=self.options.umi_tag_delim,
                cell_tag_split=self.options.cell_tag_split,
                cell_tag_delim=self.options.cell_tag_delim)
            #umi_tag_split="-",
            #    umi_tag_delim=None,
            #    cell_tag_split="-",
            #    cell_tag_delim=None)

        elif self.options.get_umi_method == "umis":
            self.barcode_getter = partial(
                get_barcode_umis,
                cell_barcode=self.options.per_cell)

        else:
            raise ValueError("Unknown UMI extraction method")

        self.read_events = collections.Counter()
        self.observed_contigs = collections.defaultdict(set)

        self.last_pos = 0
        self.last_chr = None
        self.start = 0
        self.current_chr = None

        self.reads_dict = collections.defaultdict(
            lambda: collections.defaultdict(
                lambda: collections.defaultdict(dict)))
        self.read_counts = collections.defaultdict(
            lambda: collections.defaultdict(dict))

    def update_dicts(self, read, pos, key, umi):

        # The content of the reads_dict depends on whether all reads
        # are being retained

        if self.all_reads:
            # retain all reads per key
            try:
                self.reads_dict[pos][key][umi]["count"] += 1
            except KeyError:
                self.reads_dict[pos][key][umi]["read"] = [read]
                self.reads_dict[pos][key][umi]["count"] = 1
            else:
                self.reads_dict[pos][key][umi]["read"].append(read)

        elif self.only_count_reads:
            # retain all reads per key
            try:
                self.reads_dict[pos][key][umi]["count"] += 1
            except KeyError:
                self.reads_dict[pos][key][umi]["count"] = 1

        else:
            # retain just a single read per key
            try:
                self.reads_dict[pos][key][umi]["count"] += 1
            except KeyError:
                self.reads_dict[pos][key][umi]["read"] = read
                self.reads_dict[pos][key][umi]["count"] = 1
                self.read_counts[pos][key][umi] = 0
            else:
                if self.reads_dict[pos][key][umi]["read"].mapq > read.mapq:
                    return

                if self.reads_dict[pos][key][umi]["read"].mapq < read.mapq:
                    self.reads_dict[pos][key][umi]["read"] = read
                    self.read_counts[pos][key][umi] = 0
                    return

                # TS: implemented different checks for multimapping here
                if self.options.detection_method in ["NH", "X0"]:
                    tag = self.options.detection_method
                    if (self.reads_dict[pos][key][umi]["read"].opt(tag) <
                        read.opt(tag)):
                        return
                    elif (self.reads_dict[pos][key][umi]["read"].opt(tag) >
                          read.opt(tag)):
                        self.reads_dict[pos][key][umi]["read"] = read
                        self.read_counts[pos][key][umi] = 0

                elif self.options.detection_method == "XT":
                    if self.reads_dict[pos][key][umi]["read"].opt("XT") == "U":
                        return
                    elif read.opt("XT") == "U":
                        self.reads_dict[pos][key][umi]["read"] = read
                        self.read_counts[pos][key][umi] = 0

                self.read_counts[pos][key][umi] += 1
                prob = 1.0/self.read_counts[pos][key][umi]

                if random.random() < prob:
                    self.reads_dict[pos][key][umi]["read"] = read

    def check_output(self):

        do_output = False
        out_keys = None

        if self.options.per_gene:

            if self.metacontig_contig:

                if (self.current_chr != self.last_chr and
                    (self.observed_contigs[self.last_pos] ==
                     self.metacontig_contig[self.last_pos])):
                    do_output = True
                    out_keys = [self.last_pos]

            else:
                if self.current_chr != self.last_chr:
                    do_output = True
                    out_keys = sorted(self.reads_dict.keys())

        elif self.options.whole_contig:

            if self.current_chr != self.last_chr:
                do_output = True
                out_keys = sorted(self.reads_dict.keys())

        else:

            if (self.start > (self.last_pos+1000) or
                self.current_chr != self.last_chr):

                do_output = True
                out_keys = sorted(self.reads_dict.keys())

                if self.current_chr == self.last_chr:
                    out_keys = [x for x in out_keys if x <= self.start-1000]

        return do_output, out_keys

    def __call__(self, inreads):

        for read in inreads:

            if read.is_read2:
                if self.return_read2:
                    if not read.is_unmapped or (
                            read.is_unmapped and self.return_unmapped):
                        yield read, None, "single_read"
                continue
            else:
                self.read_events['Input Reads'] += 1

            if read.is_unmapped:
                if self.options.paired:
                    if read.mate_is_unmapped:
                        self.read_events['Both unmapped'] += 1
                    else:
                        self.read_events['Read 1 unmapped'] += 1
                else:
                    self.read_events['Single end unmapped'] += 1

                if self.return_unmapped:
                    self.read_events['Input Reads'] += 1
                    yield read, None, "single_read"
                continue

            if read.mate_is_unmapped and self.options.paired:
                if not read.is_unmapped:
                    self.read_events['Read 2 unmapped'] += 1
                if self.return_unmapped:
                    yield read, None, "single_read"
                continue

            if self.options.paired:
                self.read_events['Paired Reads'] += 1

            if self.options.subset:
                if random.random() >= self.options.subset:
                    self.read_events['Randomly excluded'] += 1
                    continue

            if self.options.mapping_quality:
                if read.mapq < self.options.mapping_quality:
                    self.read_events['< MAPQ threshold'] += 1
                    continue

            self.current_chr = read.reference_name

            if self.options.per_gene:

                if self.options.per_contig:

                    if self.metacontig_contig:
                        transcript = read.reference_name
                        gene = self.contig_metacontig[transcript]
                    else:
                        gene = read.reference_name

                elif self.options.gene_tag:

                    try:
                        gene = read.get_tag(self.options.gene_tag)
                    except KeyError:
                        self.read_events['Read skipped, no tag'] += 1
                        continue

                    if re.search(self.options.skip_regex, gene):
                        self.read_events['Gene skipped - matches regex'] += 1
                        continue

                pos = gene
                key = pos

                if self.last_chr:
                    do_output, out_keys = self.check_output()
                else:
                    do_output = False

                if do_output:
                    for p in out_keys:
                        for k in sorted(self.reads_dict[p].keys()):
                            yield self.reads_dict[p][k], k, "bundle"

                        del self.reads_dict[p]

                self.last_chr = self.current_chr
                self.last_pos = pos

            else:

                start, pos, is_spliced = get_read_position(
                    read, self.options.soft_clip_threshold)

                do_output, out_keys = self.check_output()

                if do_output:
                    for p in out_keys:
                        for k in sorted(self.reads_dict[p].keys()):
                            yield self.reads_dict[p][k], k, "bundle"

                        del self.reads_dict[p]
                        if p in self.read_counts:
                            del self.read_counts[p]

                self.last_pos = self.start
                self.last_chr = self.current_chr

                if self.options.read_length:
                    r_length = read.query_length
                else:
                    r_length = 0

                key = (read.is_reverse, self.options.spliced & is_spliced,
                       self.options.paired*read.tlen, r_length)

            # get the umi +/- cell barcode and update dictionaries
            if self.options.ignore_umi:
                if self.options.per_cell:
                    umi, cell = self.barcode_getter(read)
                    umi = ""
                else:
                    umi, cell = "", ""
            else:
                umi, cell = self.barcode_getter(read)

            key = (key, cell)
            self.update_dicts(read, pos, key, umi)

            if self.metacontig_contig:
                # keep track of observed contigs for each gene
                self.observed_contigs[gene].add(transcript)

        # yield remaining bundles
        for p in sorted(self.reads_dict.keys()):
            for k in sorted(self.reads_dict[p].keys()):
                yield self.reads_dict[p][k], k, "bundle"


def get_gene_count_tab(infile,
                       umi_getter=None):

    ''' Yields the counts per umi for each gene

    umi_getter: method to get umi from read, e.g get_umi_read_id or get_umi_tag


    TODO: ADD FOLLOWING OPTION

    skip_regex: skip genes matching this regex. Useful to ignore
                unassigned reads where the 'gene' is a descriptive tag
                such as "Unassigned"

    '''

    gene = None
    counts = collections.Counter()

    for line in infile:

        values = line.strip().split("\t")

        assert len(values) == 2, "line: %s does not contain 2 columns" % line

        read_id, assigned_gene = values

        # only output when the contig changes to avoid problems with
        # overlapping genes
        if assigned_gene != gene:
            if gene:
                yield gene, counts

            gene = assigned_gene
            counts = collections.Counter()

        umi = umi_getter(read_id)
        counts[umi] += 1

    # yield final gene
    yield gene, counts


class random_read_generator:
    ''' class to generate umis at random based on the
    distributon of umis in a bamfile '''

    def __init__(self, bamfile, chrom, barcode_getter):
        inbam = pysam.Samfile(bamfile)

        if chrom:
            self.inbam = inbam.fetch(reference=chrom)
        else:
            self.inbam = inbam.fetch()

        self.umis = collections.defaultdict(int)
        self.barcode_getter = barcode_getter
        self.random_fill_size = 100000  # Higher = faster, more memory
        self.fill()

    def refill_random(self):
        ''' refill the list of random_umis '''
        self.random_umis = np.random.choice(
            list(self.umis.keys()), self.random_fill_size, p=self.prob)
        self.random_ix = 0

    def fill(self):
        ''' parse the BAM to obtain the frequency for each UMI'''
        self.frequency2umis = collections.defaultdict(list)

        for read in self.inbam:

            if read.is_unmapped:
                continue

            if read.is_read2:
                continue

            self.umis[self.barcode_getter(read)[0]] += 1

        self.umis_counter = collections.Counter(self.umis)
        total_umis = sum(self.umis_counter.values())
        U.info("total_umis %i" % total_umis)
        U.info("#umis %i" % len(self.umis_counter))

        self.prob = self.umis_counter.values()
        sum_prob = sum(self.prob)
        self.prob = [float(x) / sum_prob for x in self.prob]
        self.refill_random()

    def getUmis(self, n):
        ''' return n umis from the random_umis atr.'''
        if n < (self.random_fill_size - self.random_ix):
            barcodes = self.random_umis[self.random_ix: self.random_ix+n]
        else:
            # could use the end of the random_umis but
            # let's just make a new random_umis
            if n > self.random_fill_size:  # ensure random_umis is long enough
                self.random_fill_size = n * 2
            self.refill_random()
            barcodes = self.random_umis[self.random_ix: self.random_ix+n]

        self.random_ix += n
        return barcodes
