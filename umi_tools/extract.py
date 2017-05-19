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
barcodes and append these to the read name also.


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

       (?P<cell_1>.{8,12})(?P<discard_1>GAGTGATTGCTTGTGACGCCTT)(?P<cell_3>.{8})(?P<umi_1>.{6})T{3}.*

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
user-supplied whitelist (--whitelist-tsv) or the whitelist can be
generated using the cell barcode counts by the 'knee' method or
'network' method (--whitelist-method). The 'knee' method uses the
distribution of counts per UMI to identify the cut-off for 'true' UMIs
(the 'knee'). The 'network' method uses the 'directional' network
method. This is the same method used in dedup with UMI barcodes
and is explained more fully in the documentation for the dedup
command. A note of caution, the 'network' method assumes all false
barcodes are the result of errors. For droplet-based single cell
RNA-Sequencing, this assumption does not hold true and the 'network'
method will therefore identify far too many 'true' cell barcodes. We
recommend using the 'knee' method unless you are absolutely certain
all 'false' barcodes have been generated through errors from the
'true' barcodes. See this blog post for a more detailed explanation:

https://cgatoxford.wordpress.com/2017/05/18/estimating-the-number-of-true-cell-barcodes-in-single-cell-rna-seq/

If you are using the 'knee' method, you can supply the --plot-prefix
option to visualise the threshold set for true cell barcodes.

Cell barcodes which do not match the whitelist (user-generated or
automatically generated) can also be optionally corrected using the
--error-correct-cell option. If using the 'network' method, all cell
barcodes in the same network will be "corrected" to the most abundant
UMI. If using the 'knee' method or a user-supplied whitelist, all UMIs
which do not match the whitelist but are within
--error-correct-threshold (default 1) of a single whitelisted UMI will
be "corrected" to this UMI.


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


Command line options
--------------------

'''
import sys
import re
import regex
import collections
from six import iteritems
from scipy.stats import gaussian_kde
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt
import numpy as np
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
    import network

try:
    from umi_tools._dedup_umi import edit_distance
except:
    from _dedup_umi import edit_distance


###############################################################################
# The code for Record and fastqIterate are taken from CGAT.Fastq:
# https://github.com/CGATOxford/cgat/blob/master/CGAT/Fastq.py
###############################################################################

RANGES = {
    'phred33': (33, 77),
    'solexa': (59, 106),
    'phred64': (64, 106),
}


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


def getErrorCorrectMappings(cell_barcodes, whitelist, threshold=1):
    ''' Find the mappings between true and false cell barcodes based
    on an edit distance threshold '''

    false_to_true = {}
    true_to_false = collections.defaultdict(set)

    whitelist = set([bytes(x, encoding="utf-8") for x in whitelist])

    for cell_barcode in cell_barcodes:
        match = None
        barcode_in_bytes = bytes(cell_barcode, encoding="utf-8")
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
            false_to_true[cell_barcode] = match
            true_to_false[match].add(cell_barcode)

    return false_to_true, true_to_false


def getNetworkEstimate(cell_barcode_counts):
    ''' Use the directional network method to identify the true cell
    barcodes '''

    clusterer = network.UMIClusterer("directional")
    barcode_groups = clusterer(
        [bytes(x, encoding="utf-8") for x in cell_barcode_counts.keys()],
        {bytes(x, encoding="utf-8"): y for x, y in cell_barcode_counts.items()},
        1)

    false_to_true = {}
    true_to_false = collections.defaultdict(set)

    true_barcodes = []
    for group in barcode_groups:
        true_barcode = group[0].decode("utf-8")
        true_barcodes.append(true_barcode)
        if len(group) > 1:
            error_barcodes = group[1:]
            for error_barcode in error_barcodes:
                error_barcode = error_barcode.decode("utf-8")
                false_to_true[error_barcode] = true_barcode
                true_to_false[true_barcode].add(error_barcode)

    return true_barcodes, (false_to_true, true_to_false)


def getCellWhitelist(cell_barcode_counts,
                     method,
                     error_correct_threshold=0,
                     plotfile_prefix=None):

    if method == 'knee':
        cell_whitelist = getKneeEstimate(
            cell_barcode_counts, plotfile_prefix=plotfile_prefix)
        if error_correct_threshold > 0:
            error_correct_mappings = getErrorCorrectMappings(
                cell_barcode_counts.keys(), cell_whitelist,
                error_correct_threshold)
        else:
            error_correct_mappings = None

    elif method == 'network':
        cell_whitelist, error_correct_mappings = getNetworkEstimate(
            cell_barcode_counts)

    else:
        raise ValueError("%s is not a recognised method" % method)

    return cell_whitelist, error_correct_mappings


def getUserDefinedBarcodes(whitelist_tsv):
    cell_whitelist = []
    with U.openFile(whitelist_tsv, "r") as inf:
        for line in inf:
            cell_whitelist.append(line.strip())
    return set(cell_whitelist)


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
    parser.add_option("--whitelist-method",
                      dest="whitelist_method", type="choice",
                      choices=["knee", "network"],
                      help=("What method to use to derive the cell barcode "
                            "whitelist"))
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

    parser.set_defaults(extract_method="string",
                        filter_cell_barcodes=False,
                        whitelist_method="knee",
                        whitelist_tsv=None,
                        blacklist_tsv=None,
                        error_correct_cell=False,
                        error_correct_threshold=0,
                        pattern=None,
                        pattern2=None,
                        read2_in=None,
                        read2_out=False,
                        read2_out_only=False,
                        quality_filter_threshold=None,
                        quality_encoding=None,
                        plot_prefix=None,
                        output_whitelist=None)

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

    if (options.whitelist_method == "network" and
        options.error_correct_threshold > 0):
        U.info("The --error-correct-threshold does not change which barcodes "
               "are error corrected when using the 'network' whitelist method "
               "and --error-correct-cell. All barcodes in the same network as "
               "the true barcode are corrected to the true barcode")

    if options.pattern2:
        if not options.read2_in:
            U.error("must specify a paired fastq ``--read2-in``")

        if not options.pattern2:
            options.pattern2 = options.pattern

    extract_cell = False
    extract_umi = False

    # If the pattern is a regex we can compile the regex(es) prior to
    # ExtractFilterAndUpdateinstantiation
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

    if options.output_whitelist or options.plot_prefix:
        if not options.filter_cell_barcode or options.whitelist_tsv:
            U.error(
                "To output the automatically generated cell barcode "
                "(--output-whitelist) or plot this whitelist (--plot-prefix), "
                "you must supply the --filter-cell-barcode option and cannot "
                "supply you own whitelist with --whitelist_tsv")

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
        read1s = fastqIterate(U.openFile(options.stdin))
    else:
        read1s = fastqIterate(U.openFile(options.stdin.name))

    # set up read extractor
    ReadExtractor = ExtractFilterAndUpdate(
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

            if not options.read2_in:
                for read1 in read1s:
                    cell_barcode = ReadExtractor.getCellBarcode(read1)
                    if cell_barcode:
                        cell_barcode_counts[cell_barcode] += 1
            else:
                read2s = fastqIterate(U.openFile(options.read2_in))
                for read1, read2 in izip(read1s, read2s):
                    cell_barcode = ReadExtractor.getCellBarcode(read1, read2)
                    if cell_barcode:
                        cell_barcode_counts[cell_barcode] += 1

            if options.blacklist_tsv:
                cell_blacklist = getUserDefinedBarcodes(options.blacklist_tsv)
                for cell in cell_blacklist:
                    del cell_barcode_counts[cell]

            if options.whitelist_tsv:
                cell_whitelist = getUserDefinedBarcodes(options.whitelist_tsv)
                error_correct_mappings = getErrorCorrectMappings(
                    cell_barcode_counts.keys(), cell_whitelist,
                    options.error_correct_threshold)
            else:
                # getCellWhitelist has not been properly defined yet!
                cell_whitelist, error_correct_mappings = getCellWhitelist(
                    cell_barcode_counts,
                    options.whitelist_method,
                    options.error_correct_threshold,
                    options.plot_prefix)

            # re-make the reads1s iterator
            read1s = fastqIterate(U.openFile(options.stdin.name))

        else:
            cell_whitelist = getUserDefinedBarcodes(options.whitelist_tsv)
            error_correct_mappings = None

        if options.error_correct_cell:
            false_to_true_map, true_to_false_map = error_correct_mappings
        else:
            false_to_true_map, true_to_false_map = None, None

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

            if not new_read:
                continue

            options.stdout.write(str(new_read) + "\n")

    else:
        read2s = fastqIterate(U.openFile(options.read2_in))

        if options.read2_out:
            read2_out = U.openFile(options.read2_out, "w")

        for read1, read2 in izip(read1s, read2s):
            reads = ReadExtractor(read1, read2)

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
