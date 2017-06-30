'''
umi_methods.py - Methods for dealing with UMIs
=========================================================

:Author: Ian Sudbery, Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Python UMI

'''

import itertools
import collections
import random
import numpy as np
import pysam
import re
from scipy.stats import gaussian_kde
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt
import numpy as np

# required to make iteritems python2 and python3 compatible
from future.utils import iteritems
from builtins import dict

try:
    import umi_tools.Utilities as U
except:
    import Utilities as U

try:
    from umi_tools._dedup_umi import edit_distance
except:
    from _dedup_umi import edit_distance

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

# End of FastqIterate()
###############################################################################
###############################################################################


def getKneeEstimate(cell_barcode_counts, plotfile_prefix=None):
    ''' estimate the number of "true" cell barcodes
    input:
         cell_barcode_counts = dict(key = barcode, value = count)
         plotfile_prefix = (optional) prefix for plots

    returns:
         List of true barcodes
    '''

    # very low abundance cell barcodes are filtered out (< 0.0001 *
    # the most abundant)
    threshold = 0.0001 * cell_barcode_counts.most_common(1)[0][1]

    counts = sorted(cell_barcode_counts.values(), reverse=True)
    counts_thresh = [x for x in counts if x > threshold]
    log_counts = np.log10(counts_thresh)

    # guassian density with default optimal bw estimation
    density = gaussian_kde(log_counts)

    xx_values = 1000  # how many x values for density plot
    xx = np.linspace(log_counts.min(), log_counts.max(), xx_values)

    local_mins = argrelextrema(density(xx), np.less)[0]
    local_min = None
    for poss_local_min in local_mins:
        if poss_local_min >= 0.2 * xx_values:
            local_min = poss_local_min
            break
    if local_min is None:
        return None

    threshold = np.power(10, xx[local_min])
    final_barcodes = set([
        x for x, y in cell_barcode_counts.items() if y > threshold])

    if plotfile_prefix:

        # colour-blind friendly colours - https://gist.github.com/thriveth/8560036
        CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                          '#f781bf', '#a65628', '#984ea3',
                          '#999999', '#e41a1c', '#dede00']

        # make density plot
        fig = plt.figure()
        fig1 = fig.add_subplot(111)
        fig1.plot(xx, density(xx), 'k')
        fig1.set_xlabel("Reads (log10)")
        fig1.set_ylabel("Density")

        for pos in xx[local_mins]:
            if pos == xx[local_min]:  # selected local minima
                fig1.axvline(x=xx[local_min], color=CB_color_cycle[0],
                             label="selected local minima")
            else:
                fig1.axvline(x=pos, color=CB_color_cycle[3],
                             label="other local minima")

        fig1.legend(loc=0)

        fig.savefig("%s_cell_barcode_count_desnity.png" % plotfile_prefix)

        # make knee plot
        fig = plt.figure()
        fig2 = fig.add_subplot(111)
        fig2.plot(range(0, len(counts)), np.cumsum(counts), c="black")
        xmax = min(local_min * 5, len(counts))  # reasonable maximum x-axis value
        fig2.set_xlim((0 - 0.01 * xmax, xmax))
        fig2.set_xlabel("Rank")
        fig2.set_ylabel("Cumulative sum of reads")

        for pos in xx[local_mins]:
            threshold = np.power(10, pos)
            passing_threshold = sum([y > threshold
                                     for x, y in cell_barcode_counts.items()])

            if passing_threshold == len(final_barcodes):  # selected local minima
                fig2.axvline(x=passing_threshold, color=CB_color_cycle[0],
                             label="selected local minima")
            else:
                fig2.axvline(x=passing_threshold, color=CB_color_cycle[3],
                             label="other local minima")

        fig2.legend(loc=0)
        fig.savefig("%s_cell_barcode_knee.png" % plotfile_prefix)

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


# TS: This is not being used since the network methods are not
# appropriate to identify cell barcodes as they currently stand

# def getNetworkEstimate(cell_barcode_counts):
#     ''' Use the directional network method to identify the true cell
#     barcodes '''

#     clusterer = network.UMIClusterer("directional")
#     barcode_groups = clusterer(
#         [bytes(x, encoding="utf-8") for x in cell_barcode_counts.keys()],
#         {bytes(x, encoding="utf-8"): y for x, y in cell_barcode_counts.items()},
#         1)

#     false_to_true = {}
#     true_to_false = collections.defaultdict(set)

#     true_barcodes = []
#     for group in barcode_groups:
#         true_barcode = group[0].decode("utf-8")
#         true_barcodes.append(true_barcode)
#         if len(group) > 1:
#             error_barcodes = group[1:]
#             for error_barcode in error_barcodes:
#                 error_barcode = error_barcode.decode("utf-8")
#                 false_to_true[error_barcode] = true_barcode
#                 true_to_false[true_barcode].add(error_barcode)

#     return true_barcodes, (false_to_true, true_to_false)


def getCellWhitelist(cell_barcode_counts,
                     error_correct_threshold=0,
                     plotfile_prefix=None):

    cell_whitelist = getKneeEstimate(
        cell_barcode_counts, plotfile_prefix=plotfile_prefix)
    if error_correct_threshold > 0:
        true_to_false_map = getErrorCorrectMapping(
            cell_barcode_counts.keys(), cell_whitelist,
            error_correct_threshold)
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


def get_barcode_tag(read, cell_barcode=False, umi_tag='RX', cell_tag=None):
    ''' extract the umi +/- cell barcode from the specified tag '''

    try:
        # 10X pipelines append a 'GEM' tag to the UMI, e.g
        # AGAGSGATAGATA-1
        if cell_barcode:
            umi = read.get_tag(umi_tag).split("-")[0].encode('utf-8')
            cell = read.get_tag(cell_tag).split("-")[0].encode('utf-8')
        else:
            umi = read.get_tag(umi_tag).split("-")[0].encode('utf-8')
            cell = None
        return umi, cell

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


def extractSeqAndQuals(seq, quals, umi_bases, cell_bases, discard_bases):
    '''Remove selected bases from seq and quals'''

    new_seq = ""
    new_quals = ""
    umi_quals = ""
    cell_quals = ""

    ix = 0
    for base, qual in zip(seq, quals):
        if ((ix not in discard_bases) and
            (ix not in cell_bases) and
            (ix not in umi_bases)):
            new_quals += qual
            new_seq += base
        elif ix in cell_bases:
            cell_quals += qual
        elif ix in umi_bases:
            umi_quals += qual

        ix += 1

    return new_seq, new_quals, umi_quals, cell_quals


def ExtractBarcodes(read, match,
                    extract_umi=False,
                    extract_cell=False,
                    discard=False):
    '''Extract the cell and umi barcodes using a regex.match object

    inputs:

    - read 1 and read2 = Record objects
    - match = regex.match object
    - extract_umi and extract_cell = switches to determine whether these
                                     barcodes should be extracted

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
        read.seq, read.quals, umi_bases, cell_bases, discard_bases)

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
            bc2, sequence2 = self.extract(read2.seq)
            bc_qual2, seq_qual2 = self.extract(read2.quals)
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
                 extract_umi=True, discard=True)
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
            bc2, sequence2 = self.extract(read2.seq)
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
                 filter_cell_barcode=False):

        self.read_counts = collections.Counter()
        self.method = method
        self.pattern = pattern
        self.pattern2 = pattern2
        self.extract_cell = extract_cell
        self.quality_encoding = quality_encoding
        self.quality_filter_threshold = quality_filter_threshold
        self.quality_filter_mask = quality_filter_mask
        self.filter_cell_barcodes = filter_cell_barcode

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

    def getReadCounts(self):
        return self.read_counts


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


def get_read_position(read, soft_clip_threshold):
    ''' '''
    is_spliced = False

    if read.is_reverse:
        pos = read.aend
        if read.cigar[-1][0] == 4:
            pos = pos + read.cigar[-1][1]
        start = read.pos

        if ('N' in read.cigarstring or
            (read.cigar[0][0] == 4 and
             read.cigar[0][1] > soft_clip_threshold)):
            is_spliced = True
    else:
        pos = read.pos
        if read.cigar[0][0] == 4:
            pos = pos - read.cigar[0][1]
        start = pos

        if ('N' in read.cigarstring or
            (read.cigar[-1][0] == 4 and
             read.cigar[-1][1] > soft_clip_threshold)):
            is_spliced = True

    return start, pos, is_spliced


def getMetaContig2contig(bamfile, gene_transcript_map):
    ''' '''
    references = bamfile.references
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


def get_bundles(inreads,
                ignore_umi=False,
                subset=None,
                quality_threshold=0,
                paired=False,
                spliced=False,
                soft_clip_threshold=0,
                per_contig=False,
                gene_tag=None,
                skip_regex=None,
                whole_contig=False,
                read_length=False,
                detection_method=False,
                barcode_getter=None,
                all_reads=False,
                return_read2=False,
                return_unmapped=False):

    ''' Returns a dictionary of dictionaries, representing the unique reads at
    a position/spliced/strand combination. The key to the top level dictionary
    is a umi. Each dictionary contains a "read" entry with the best read, and a
    count entry with the number of reads with that position/spliced/strand/umi
    combination

    ignore_umi: don't include the umi in the dict key

    subset: randomly exclude 1-subset fraction of reads

    quality_threshold: exclude reads with MAPQ below this

    paired: input is paired

    spliced: include the spliced/not-spliced status in the dict key

    soft_clip_threshold: reads with less than this 3' soft clipped are
    treated as spliced

    per_contig: use just the umi and contig as the dict key

    gene_tag: use just the umi and gene as the dict key. Get the gene
    id from the this tag

    skip_regex: skip genes matching this regex. Useful to ignore
    unassigned reads where the 'gene' is a descriptive tag such as
    "Unassigned"

    whole_contig: read the whole contig before yielding a bundle

    read_length: include the read length in the dict key

    detection_method: which method to use to detect multimapping
    reads. options are NH", "X0", "XT". defaults to False (just select
    the 'best' read by MAPQ)

    barcode_getter: method to get umi from read, e.g get_umi_read_id or get_umi_tag

    all_reads: if true, return all reads in the dictionary. Else,
    return the 'best' read (using MAPQ +/- multimapping) for each key

    return_read2: Return read2s immediately as a single read

    return_unmapped: Return unmapped reads immediately as a single read
    '''

    last_pos = 0
    last_chr = ""
    reads_dict = collections.defaultdict(
        lambda: collections.defaultdict(
            lambda: collections.defaultdict(dict)))
    read_counts = collections.defaultdict(
        lambda: collections.defaultdict(dict))

    read_events = collections.Counter()

    for read in inreads:

        if read.is_read2:
            if return_read2:
                if not read.is_unmapped or (read.is_unmapped and return_unmapped):
                    yield read, read_events, 'single_read'
            continue
        else:
            read_events['Input Reads'] += 1

        if read.is_unmapped:
            if paired:
                if read.mate_is_unmapped:
                    read_events['Both unmapped'] += 1
                else:
                    read_events['Read 1 unmapped'] += 1
            else:
                read_events['Single end unmapped'] += 1

            if return_unmapped:
                read_events['Input Reads'] += 1
                yield read, read_events, 'single_read'
            continue

        if read.mate_is_unmapped and paired:
            if not read.is_unmapped:
                read_events['Read 2 unmapped'] += 1
            if return_unmapped:
                yield read, read_events, 'single_read'
            continue

        if paired:
            read_events['Paired Reads'] += 1

        if subset:
            if random.random() >= subset:
                read_events['Randomly excluded'] += 1
                continue

        if quality_threshold:
            if read.mapq < quality_threshold:
                read_events['< MAPQ threshold'] += 1
                continue

        # TS - some methods require deduping on a per contig or per
        # gene basis. To fit in with current workflow, simply assign
        # pos and key as contig

        if per_contig or gene_tag:

            if per_contig:
                pos = read.tid
                key = pos
            elif gene_tag:
                pos = read.get_tag(gene_tag)
                key = pos
                if re.search(skip_regex, pos):
                    continue

            if not pos == last_chr:

                out_keys = list(reads_dict.keys())

                for p in out_keys:
                    for bundle in reads_dict[p].values():
                        yield bundle, read_events, 'mapped'
                    del reads_dict[p]
                    del read_counts[p]

                last_chr = pos

        else:

            start, pos, is_spliced = get_read_position(
                read, soft_clip_threshold)

            if whole_contig:
                do_output = not read.tid == last_chr
            else:
                do_output = start > (last_pos+1000) or not read.tid == last_chr

            if do_output:
                if not read.tid == last_chr:
                    out_keys = list(reads_dict.keys())
                else:
                    out_keys = [x for x in reads_dict.keys() if x <= start-1000]

                for p in out_keys:
                    for bundle in reads_dict[p].values():
                        yield bundle, read_events, 'mapped'
                    del reads_dict[p]
                    del read_counts[p]

                last_pos = start
                last_chr = read.tid

            if read_length:
                r_length = read.query_length
            else:
                r_length = 0

            key = (read.is_reverse, spliced & is_spliced,
                   paired*read.tlen, r_length)

        if ignore_umi:
            umi = ""
        else:
            umi, cell = barcode_getter(read)  # cell is always None here

        # The content of the reads_dict depends on whether all reads
        # are being retained

        if all_reads:
            # retain all reads per key
            try:
                reads_dict[pos][key][umi]["count"] += 1
            except KeyError:
                reads_dict[pos][key][umi]["read"] = [read]
                reads_dict[pos][key][umi]["count"] = 1
                read_counts[pos][key][umi] = 0
            else:
                reads_dict[pos][key][umi]["read"].append(read)

        else:
            # retain just a single read per key
            try:
                reads_dict[pos][key][umi]["count"] += 1
            except KeyError:
                reads_dict[pos][key][umi]["read"] = read
                reads_dict[pos][key][umi]["count"] = 1
                read_counts[pos][key][umi] = 0
            else:
                if reads_dict[pos][key][umi]["read"].mapq > read.mapq:
                    continue

                if reads_dict[pos][key][umi]["read"].mapq < read.mapq:
                    reads_dict[pos][key][umi]["read"] = read
                    read_counts[pos][key][umi] = 0
                    continue

                # TS: implemented different checks for multimapping here
                if detection_method in ["NH", "X0"]:
                    tag = detection_method
                    if reads_dict[pos][key][umi]["read"].opt(tag) < read.opt(tag):
                        continue
                    elif reads_dict[pos][key][umi]["read"].opt(tag) > read.opt(tag):
                        reads_dict[pos][key][umi]["read"] = read
                        read_counts[pos][key][umi] = 0

                elif detection_method == "XT":
                    if reads_dict[pos][key][umi]["read"].opt("XT") == "U":
                        continue
                    elif read.opt("XT") == "U":
                        reads_dict[pos][key][umi]["read"] = read
                        read_counts[pos][key][umi] = 0

                read_counts[pos][key][umi] += 1
                prob = 1.0/read_counts[pos][key][umi]

                if random.random() < prob:
                    reads_dict[pos][key][umi]["read"] = read

    # yield remaining bundles
    for p in reads_dict:
        for bundle in reads_dict[p].values():
            yield bundle, read_events, 'mapped'


def get_gene_count(inreads,
                   subset=None,
                   quality_threshold=0,
                   paired=False,
                   per_contig=False,
                   gene_tag=None,
                   metacontig2contig=None,
                   skip_regex=None,
                   barcode_getter=None):

    ''' Yields the counts per umi for each gene

    ignore_umi: don't include the umi in the dict key

    subset: randomly exclude 1-subset fraction of reads

    quality_threshold: exclude reads with MAPQ below this

    paired: input is paired

    per_contig: use just the umi and contig as the dict key

    gene_tag: use just the umi and gene as the dict key. Get the gene
              id from the this tag

    skip_regex: skip genes matching this regex. Useful to ignore
                unassigned reads where the 'gene' is a descriptive tag
                such as "Unassigned"

    barcode_getter: method to get umi from read, e.g get_umi_read_id or get_umi_tag
    '''

    previous_chr = None
    previous_gene = None
    gene = ""

    # make an empty counts_dict counter
    counts_dict = collections.defaultdict(
        lambda: collections.defaultdict(
            lambda: collections.defaultdict(dict)))

    read_events = collections.Counter()

    if metacontig2contig:
        observed_contigs = collections.defaultdict(set)

    for read in inreads:

        if read.is_read2:
            continue
        else:
            read_events['Input Reads'] += 1

        if read.is_unmapped:
            read_events['Skipped - Unmapped Reads'] += 1
            continue

        if paired:
            read_events['Paired Reads'] += 1

        if subset:
            if random.random() >= subset:
                read_events['Skipped - Randomly excluded'] += 1
                continue

        if quality_threshold:
            if read.mapq < quality_threshold:
                read_events['Skipped - < MAPQ threshold'] += 1
                continue

        if per_contig:
            gene = read.tid

        elif gene_tag:
            gene = read.get_tag(gene_tag)

            if metacontig2contig:
                # keep track of observed contigs for each gene
                observed_contigs[gene].add(read.reference_name)

            if re.search(skip_regex, gene):
                read_events['Skipped - matches --skip-tags-regex'] += 1
                continue

        # TS: We can safely yield when the contig changes to keep the
        # size of the count_dict smaller and avoid problems with
        # overlapping genes

        if read.reference_name != previous_chr and previous_chr:

            # TS: However, when the BAM contains reads aligned to
            # transcripts (i.e no gene_transcript_map), we also need
            # to check all the transcripts for a gene have been
            # observed, otherwise we may yield for a gene more than once
            if metacontig2contig:

                if observed_contigs[previous_gene] == metacontig2contig[previous_gene]:
                    # yield only the counts for this gene
                    for cell in counts_dict[previous_gene]:
                        count = counts_dict[previous_gene][cell]
                        yield previous_gene, cell, count, read_events

                    del counts_dict[previous_gene]

            else:

                # yield all gene counts
                for gene in counts_dict:
                    for cell in counts_dict[gene]:
                        count = counts_dict[gene][cell]
                        yield gene, cell, count, read_events

                # make a new empty counts_dict counter
                counts_dict = collections.defaultdict(
                    lambda: collections.defaultdict(
                        lambda: collections.defaultdict(dict)))

        umi, cell = barcode_getter(read)
        try:
            counts_dict[gene][cell][umi]["count"] += 1
        except KeyError:
            counts_dict[gene][cell][umi]["count"] = 1

        previous_chr = read.reference_name
        previous_gene = gene

    # yield gene counts
    for gene in counts_dict:
        for cell in counts_dict[gene]:
            count = counts_dict[gene][cell]
            yield gene, cell, count, read_events


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
        self.fill()

    def fill(self):

        self.frequency2umis = collections.defaultdict(list)

        for read in self.inbam:

            if read.is_unmapped:
                continue

            if read.is_read2:
                continue

            self.umis[self.barcode_getter(read)[0]] += 1

        self.umis_counter = collections.Counter(self.umis)
        total_umis = sum(self.umis_counter.values())

        for observed_umi, freq in iteritems(self.umis_counter):
            self.frequency2umis[freq+0.0/total_umis].append(observed_umi)

        self.frequency_counter = collections.Counter(self.umis_counter.values())
        self.frequency_prob = [(float(x)/total_umis)*y for x, y in
                               iteritems(self.frequency_counter)]

    def getUmis(self, n):
        '''get n umis at random'''

        umi_sample = []

        frequency_sample = np.random.choice(
            list(self.frequency_counter.keys()), n, p=self.frequency_prob)

        for frequency in frequency_sample:
            umi_sample.append(np.random.choice(self.frequency2umis[frequency]))

        return umi_sample
