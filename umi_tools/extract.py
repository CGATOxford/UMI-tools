'''
extract.py - Extract UMI from fastq
====================================================

:Author: Ian Sudbery
:Release: $Id$
:Date: |today|
:Tags: Python UMI

Purpose
-------

Extract UMI barcode from a read and add it to the read name, leaving
any sample barcode in place. Can deal with paired end reads and UMIs
split across the paired ends

Options
-------



-L (string, filename)
       Specify a log file to retain logging information and final statistics

--supress-stats
Supress logging of summary statistics

Usage:
------

For single ended reads:
        umi_tools extract --barcode-regex=[REGEX] -L extract.log [OPTIONS]

reads from stdin and outputs to stdout.

For paired end reads:
    xumi_tools extract --barcode-regex=[REGEX] --barcode-regex2=[REGEX] --read2-in=[FASTQIN] --read2-out=[FASTQOUT] -L extract.log [OPTIONS]

reads end one from stdin and end two from FASTQIN and outputs end one to stdin
and end two to FASTQOUT.

Command line options
--------------------

'''
import sys
import re
import regex
import collections
from six import iteritems

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
            raise ValueError("parsing error: expected '@' in line %s" % line1)
        line2 = convert2string(infile.readline())
        line3 = convert2string(infile.readline())
        if not line3.startswith('+'):
            raise ValueError("parsing error: expected '+' in line %s" % line3)
        line4 = convert2string(infile.readline())
        # incomplete entry
        if not line4:
            raise ValueError("incomplete entry for %s" % line1)

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
        raise ValueError("must set either extract_cell and/or extract_umi to true")

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


def getCellWhitelist(cell_barcode_counts,
                     method,
                     error_correct_cell=False,
                     hamming_threshold=1):

    # we need to define cell whitelist using method.
    # line below is placeholder.
    cell_whitelist = list(cell_barcode_counts.keys())
    error_correct_mapping = {}

    # error_correct_mapping is a dict mapping barcodes not in the
    # final whitelist to cell barcodes in the whitelist. For 'knee'
    # method we use a Hamming distance threshold. For network-based
    # methods, we use the networks to derive this

    # # This should work to create the error_correct_mapping post-knee method
    # if error_correct_cell:
    #     error_correct_mapping = {}
    #     for cell_barcode in cell_barcode_counts.keys():
    #         if cell_barcode not in cell_whitelist:
    #             near_matches = []
    #             for w_cell in cell_whitelist:
    #                 if edit_distance(cell_barcode.encode('utf-8'),
    #                                  w_cell.encode('utf-8')) < hamming_threshold:
    #                     near_matches.append(w_cell)
    #                     if len(near_matches) > 1:
    #                         break
    #             if len(near_matches) == 1:
    #                 error_correct_mapping[cell_barcode] = near_matches[0]
    # else:
    #     error_correct_mapping = {}

    return cell_whitelist, error_correct_mapping


def getUserDefinedWhitelist(whitelist_tsv):
    cell_whitelist = []
    with U.openFile(whitelist_tsv, "r") as inf:
        for line in inf:
            cell_whitelist.append(line.strip())
    return cell_whitelist


def get_below_threshold(umi_quals, quality_encoding,  quality_filter_threshold):
    '''test whether the umi_quals are below the threshold'''
    umi_quals = [x - RANGES[quality_encoding][0] for x in map(ord, umi_quals)]
    below_threshold = [x < quality_filter_threshold for x in umi_quals]
    return (below_threshold)


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

    def __init__(self,
                 pattern=None,
                 pattern2=None,
                 extract_cell=False,
                 quality_encoding=None,
                 quality_filter_threshold=False,
                 quality_filter_mask=False,
                 filter_cell_barcode=False,
                 error_correct_cell=False,
                 cell_whitelist=None,
                 error_correct_mapping=None):
        self.read_counts = collections.Counter()
        self.pattern = pattern
        self.pattern2 = pattern2
        self.extract_cell = extract_cell
        self.quality_encoding = quality_encoding
        self.quality_filter_threshold = quality_filter_threshold
        self.quality_filter_mask = quality_filter_mask
        self.filter_cell_barcodes = filter_cell_barcode
        self.error_correct_cell = error_correct_cell
        self.cell_whitelist = cell_whitelist
        self.error_correct_mapping = error_correct_mapping

    def __call__(self, read1, read2=None):

        self.read_counts['Input Reads'] += 1

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

        if self.pattern:
            (cell_barcode, cell_barcode_quals,
             umi, umi_quals,
             new_seq, new_quals) = ExtractBarcodes(
                 read1, match, extract_cell=self.extract_cell,
                 extract_umi=True, discard=True)
        else:
            cell_barcode, cell_barcode_quals, umi, umi_quals = ("",)*4

        if read2 and self.pattern2:
            (cell_barcode2, cell_barcode_quals2,
             umi2, umi_quals2,
             new_seq2, new_quals2) = ExtractBarcodes(
                 read2, match2, extract_cell=self.extract_cell,
                 extract_umi=True, discard=True)

            cell_barcode += cell_barcode2
            cell_barcode_quals += cell_barcode_quals2
            umi += umi2
            umi_quals += umi_quals2

        if self.quality_filter_threshold:
            if umi_below_threshold(
                    umi_quals, self.quality_encoding,
                    self.quality_filter_threshold):
                self.read_counts['filtered: umi-quality'] += 1
                return None

        if self.quality_filter_mask:
            masked_umi = mask_umi(umi, umi_quals,
                                  self.quality_encoding,
                                  self.quality_filter_mask)
            if masked_umi != umi:
                umi = masked_umi
                self.read_counts['UMI masked'] += 1

        if self.filter_cell_barcodes:
            if cell_barcode not in self.cell_whitelist:
                if self.error_correct_cell:
                    if self.error_correct_mapping:
                        if cell_barcode in self.error_correct_mapping:
                            cell_barcode = self.error_correct_mapping[cell_barcode]
                            self.read_counts['False cell barcode. Error-corrected'] += 1
                        else:
                            self.read_counts['Filtered cell barcode. Not correctable'] += 1
                            return None
                    else:
                        pass  # need to implement error correction on the fly
                else:
                    self.read_counts['Filtered cell barcode'] += 1
                    return None

        self.read_counts['Reads output'] += 1

        new_identifier = addBarcodesToIdentifier(
            read1, umi, cell_barcode)
        read1.identifier = new_identifier
        if self.pattern:  # seq and quals need to be updated
            read1.seq = new_seq
            read1.quals = new_quals

        if read2:
            new_identifier2 = addBarcodesToIdentifier(
                read2, umi, cell_barcode)
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


def getCellBarcode(read1, read2=None, pattern=None, pattern2=None):

    if read2 is None:
        match = pattern.match(read1.seq)
        if match:
            cell_barcode = ExtractBarcodes(
                read1, match, extract_cell=True)[0]
            return cell_barcode
        else:
            return None

    else:

        match1, match2 = None, None

        if pattern:
            match1 = pattern.match(read1.seq)

        if pattern2:
            match2 = pattern2.match(read2.seq)

        # check matches have been made
        if not ((pattern and not match1) or
                (pattern2 and not match2)):
            cell_barcode1, cell_barcode2 = "", ""

            if pattern:
                cell_barcode1 = ExtractBarcodes(
                    read1, match1, extract_cell=True)[0]
            if pattern2:
                cell_barcode2 = ExtractBarcodes(
                    read2, match2, extract_cell=True)[0]

            cell_barcode = cell_barcode1 + cell_barcode2

            return cell_barcode
        else:
            return None


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = U.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("--barcode-regex", dest="barcode_regex", type="string",
                      help="barcode regex")
    parser.add_option("--barcode-regex2", dest="barcode_regex2", type="string",
                      help="barcode regex for paired read")
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
    parser.add_option("--filter-cell-barcode",
                      dest="filter_cell_barcode",
                      action="store_true",
                      help="Filter the cell barcodes")
    parser.add_option("--error-correct-cell",
                      dest="error_correct_cell",
                      action="store_true",
                      help="Error correct the cell barcode")
    parser.add_option("--error-correct-threshold",
                      dest="error_correct_threshold",
                      type="int",
                      help="Hamming distance allowed for correction")
    parser.add_option("--whitelist-method",
                      dest="whitelist_method", type="choice",
                      choices=["knee", "directional", "adjacency", "cluster"],
                      help=("What method to use to derive the cell barcode whitelist"))
    parser.add_option("--whitelist-tsv",
                      dest="whitelist_tsv", type="string",
                      help=("A whitelist of accepted cell barcodes"))
    parser.add_option("--supress-stats", dest="stats", action="store_false",
                      help="Suppress the writing of stats to the log")

    parser.set_defaults(stats=True,
                        filter_cell_barcodes=False,
                        whitelist_method="directional",
                        whitelist_tsv=None,
                        error_correct_cell=False,
                        error_correct_threshold=1,
                        barcode_regex=None,
                        barcode_regex2=None,
                        read2_in=None,
                        read2_out=False,
                        read2_out_only=False,
                        quality_filter_threshold=None,
                        quality_encoding=None)

    # add common options (-h/--help, ...) and parse command line

    (options, args) = U.Start(parser, argv=argv)

    if options.quality_filter_threshold or options.quality_filter_mask:
        if not options.quality_encoding:
            raise ValueError("must provide a quality encoding (--quality-"
                             "encoding) to filter UMIs by quality (--quality"
                             "-filter-threshold) or mask low quality bases "
                             "with (--quality-filter-mask)")

    if options.barcode_regex:
        try:
            pattern = regex.compile(options.barcode_regex)
        except regex.error:
            raise ValueError("barcode_regex '%s' is not a "
                             "valid regex" % options.barcode_regex)
    else:
        pattern = None
        if not options.read2_in:
            raise ValueError("Must supply --barcode-regex if single-end")
        if not options.barcode_regex2:
            raise ValueError("Must supply --barcode-regex2 if paired-end "
                             "and no --barcode-regex")

    if options.barcode_regex2:
        try:
            pattern2 = regex.compile(options.barcode_regex2)
        except regex.Error:
            raise ValueError("barcode_regex2 '%s' is not a "
                             "valid regex" % options.barcode_regex2)
    else:
        pattern2 = None

    # check whether the regex contains a umi group(s) and cell groups(s)
    extract_cell = False
    extract_umi = False
    if pattern:
        for group in pattern.groupindex:
            if group.startswith("cell_"):
                extract_cell = True
            elif group.startswith("umi_"):
                extract_umi = True
    if pattern2:
        for group in pattern2.groupindex:
            if group.startswith("cell_"):
                extract_cell = True
            elif group.startswith("umi_"):
                extract_umi = True

    if not extract_umi:
        raise ValueError("barcode regex(es) does not include any umi groups "
                         "(starting with 'umi_') %s" % options.barcode_regex)

    if options.stdin == sys.stdin:
        if not options.whitelist_tsv and options.filter_cell_barcode:
            raise ValueError("cannot support reading from stdin if correcting "
                             "cell barcode without providing a whitelist")
        read1s = fastqIterate(U.openFile(options.stdin))
    else:
        read1s = fastqIterate(U.openFile(options.stdin.name))

    if options.filter_cell_barcode:
        if not options.whitelist_tsv:
            cell_barcode_counts = collections.Counter()

            if not options.read2_in:
                for read1 in read1s:
                    cell_barcode = getCellBarcode(read1, pattern=pattern)
                    if cell_barcode:
                        cell_barcode_counts[cell_barcode] += 1
            else:
                read2s = fastqIterate(U.openFile(options.read2_in))
                for read1, read2 in izip(read1s, read2s):
                    cell_barcode = getCellBarcode(read1, read2, pattern, pattern2)
                    if cell_barcode:
                        cell_barcode_counts[cell_barcode] += 1

            # getCellWhitelist has not been properly defined yet!
            cell_whitelist, error_correct_mapping = getCellWhitelist(
                cell_barcode_counts,
                options.whitelist_method,
                options.error_correct_cell,
                options.error_correct_threshold)

            # re-make the reads1s iterator
            read1s = fastqIterate(U.openFile(options.stdin.name))

        else:
            cell_whitelist = getUserDefinedWhitelist(options.whitelist_tsv)
            error_correct_mapping = None
    else:
        cell_whitelist, error_correct_mapping = None, None

    # set up read extractor
    ReadExtractor = ExtractFilterAndUpdate(
        pattern,
        pattern2,
        extract_cell,
        options.quality_encoding,
        options.quality_filter_threshold,
        options.quality_filter_mask,
        options.filter_cell_barcode,
        options.error_correct_cell,
        cell_whitelist,
        error_correct_mapping)

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

    if options.stats:
        pass
        # Need to define stats.
        # Would be especially good to have some stats regarding the cell barcode filtering
        # options.stdlog.write("\t".join(["Barcode", "UMI", "Sample", "Count"]) + "\n")
        # for id in processor.bc_count:
        #     options.stdlog.write("\t".join(id+(str(processor.bc_count[id]),)) + "\n")

    for k, v in ReadExtractor.getReadCounts().most_common():
        U.info("%s: %s" % (k, v))

    U.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
