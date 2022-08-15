'''
extract_methods.py - Methods for extracting UMIs and cell barcodes
===================================================================

'''
import collections
import umi_tools.Utilities as U
from umi_tools.umi_methods import RANGES


def addBarcodesToIdentifier(read, UMI, cell, umi_separator):
    '''extract the identifier from a read and append the UMI and
    cell barcode before the first space'''

    read_id = read.identifier.split(" ")

    if cell == "":
        read_id[0] = read_id[0] + umi_separator + UMI
    else:
        read_id[0] = read_id[0] + umi_separator + cell + umi_separator + UMI

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
                if not self.either_read:
                    return None
            else:
                self.read_counts['regex matches read1'] += 1

        if read2 and self.pattern2:
            match2 = self.pattern2.match(read2.seq)
            if not match2:
                self.read_counts['regex does not match read2'] += 1
                # if barcodes can be on either read, check if read1
                # match also failed
                if self.either_read:
                    if not match:
                        return None
                else:
                    return None
            else:
                self.read_counts['regex matches read2'] += 1

        # now extract barcodes
        if self.either_read:  # extract depending on match and match2

            if match and match2 and self.either_read_resolve == "discard":
                self.read_counts['regex matches both. discarded'] += 1
                return None

            if match:
                (cell, cell_quals, umi, umi_quals,
                 new_seq, new_quals) = ExtractBarcodes(
                     read1, match, extract_cell=self.extract_cell,
                     extract_umi=True, discard=True, retain_umi=self.retain_umi)
            if match2:
                (cell2, cell_quals2, umi2, umi_quals2,
                 new_seq2, new_quals2) = ExtractBarcodes(
                     read2, match2, extract_cell=self.extract_cell,
                     extract_umi=True, discard=True, retain_umi=self.retain_umi)

            if match and match2:
                if self.either_read_resolve == "quality":
                    read1_min_quals = min(
                        [x - RANGES[self.quality_encoding][0] for
                         x in map(ord, umi_quals)])

                    read2_min_quals = min(
                        [x - RANGES[self.quality_encoding][0] for
                         x in map(ord, umi_quals2)])

                    # Note that UMI and UMI quals are defined by best
                    # match but UMI is removed from both reads
                    if read1_min_quals >= read2_min_quals:
                        self.read_counts['regex matches both. read1 used'] += 1
                        return(cell, umi, umi_quals, new_seq, new_quals, new_seq2, new_quals2)
                    else:
                        self.read_counts['regex matches both. read2 used'] += 1
                        return(cell2, umi2, umi_quals2, new_seq, new_quals, new_seq2, new_quals2)

                else:
                    raise ValueError("unexpected value for either_read_resolve")

            elif match:
                self.read_counts['regex only matches read1'] += 1
                return(cell, umi, umi_quals,
                       new_seq, new_quals, "", "")
            else:
                self.read_counts['regex only matches read2'] += 1
                return(cell2, umi2, umi_quals2, "", "", new_seq2, new_quals2)

        else:  # extract dependening on patterns/reads provided
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
                     extract_umi=True, discard=True, retain_umi=self.retain_umi)

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

    def filterUMIBarcode(self, umi):
        ''' Filter out umi barcodes not in the whitelist, with optional
        umi barcode error correction and update read_counts and
        umi_whitelist_counts counters'''

        if umi not in self.umi_whitelist:  # if umi not in whitelist

            corrected_umi = None  # need to try and correct UMI

            if self.umi_false_to_true_map:  # if there is a error correction map

                if umi in self.umi_false_to_true_map:  # and umi in map

                    # Get corrected UMI
                    # Will be None if not correctable to single whitelist UMI
                    corrected_umi = self.umi_false_to_true_map[umi]

                    if corrected_umi:  # if correctable
                        self.read_counts['False UMI barcode. Error-corrected'] += 1
                        self.umi_whitelist_counts[corrected_umi]["error"] += 1

                    else:  # Not correctable to single whitelist UMI
                        self.read_counts[
                            ("False UMI barcode. Not correctable - "
                             "within threshold to more than one whitelisted UMI")] += 1

                else:  # Not correctable to any whitelist UMI
                    self.read_counts[
                        ("False UMI barcode. Not correctable - not within "
                         "threshold to whitelisted UMI")] += 1

            else:  # no error correction map so log simply as filtered
                self.read_counts['Filtered umi barcode'] += 1

            return corrected_umi

        else:  # umi in whitelist
            self.umi_whitelist_counts[umi]["no_error"] += 1

            return umi

    def __init__(self,
                 method="string",
                 pattern=None,
                 pattern2=None,
                 prime3=False,
                 extract_cell=False,
                 quality_encoding=None,
                 quality_filter_threshold=False,
                 quality_filter_mask=False,
                 filter_umi_barcode=False,
                 filter_cell_barcode=False,
                 retain_umi=False,
                 either_read=False,
                 either_read_resolve="discard",
                 umi_separator="_"):

        self.method = method
        self.read_counts = collections.Counter()
        self.pattern = pattern
        self.pattern2 = pattern2
        self.extract_cell = extract_cell
        self.quality_encoding = quality_encoding
        self.quality_filter_threshold = quality_filter_threshold
        self.quality_filter_mask = quality_filter_mask
        self.filter_umi_barcodes = filter_umi_barcode
        self.filter_cell_barcodes = filter_cell_barcode
        self.retain_umi = retain_umi
        self.either_read = either_read
        self.either_read_resolve = either_read_resolve

        self.umi_whitelist = None  # These will be updated if required
        self.umi_false_to_true_map = None  # These will be updated if required
        self.umi_whitelist_counts = None  # These will be updated if required
        self.umi_separator = umi_separator

        self.cell_whitelist = None  # These will be updated if required
        self.false_to_true_map = None  # These will be updated if required
        self.cell_blacklist = None  # These will be updated if required

        # If the pattern is a string we can identify the position of
        # the cell and umi bases at instantiation
        if self.method == "string":
            if prime3:
                self.extract = self._extract_3prime
                self.joiner = self._joiner_3prime
            else:
                self.extract = self._extract_5prime
                self.joiner = self._joiner_5prime

            if pattern:
                if len(pattern.replace("N", "").replace("X", "").replace("C", "")) > 0:
                    raise ValueError("barcode pattern (%s) should only contain "
                                     "N/X/C characters" % pattern)
                self.pattern_length = len(pattern)
                self.umi_bases = [x for x in range(len(pattern)) if pattern[x] == "N"]
                self.bc_bases = [x for x in range(len(pattern)) if pattern[x] == "X"]
                self.cell_bases = [x for x in range(len(pattern)) if pattern[x] == "C"]

            if pattern2:
                if len(pattern2.replace("N", "").replace("X", "").replace("C", "")) > 0:
                    raise ValueError("barcode pattern2 (%s) should only contain "
                                     "N/X/C characters" % pattern2)
                self.pattern_length2 = len(pattern2)
                self.umi_bases2 = [x for x in range(len(pattern2))
                                   if pattern2[x] == "N"]
                self.bc_bases2 = [x for x in range(len(pattern2))
                                  if pattern2[x] == "X"]
                self.cell_bases2 = [x for x in range(len(pattern2))
                                    if pattern2[x] == "C"]

            self.getCellBarcode = self._getCellBarcodeString
            self.getBarcodes = self._getBarcodesString

        elif self.method == "regex":
            self.getCellBarcode = self._getCellBarcodeRegex
            self.getBarcodes = self._getBarcodesRegex

    def getReadCounts(self):
        return self.read_counts

    def __call__(self, read1, read2=None):

        # Check for string ensures sensible error message rather
        # than out of bounds exception (#424)
        # Check for regex is more complex.
        # Reads too short with regex pattern will silently fail to match
        if self.method == "string":
            if(len(read1.seq) < len(self.pattern)):
                raise ValueError('Read sequence: %s is shorter than pattern: %s' % (
                    read1.seq, self.pattern))

            if(read2 is not None and self.pattern2 is not None and
               len(read2.seq) < len(self.pattern2)):
                raise ValueError('Read2 sequence: %s is shorter than pattern2: %s' % (
                    read2.seq, self.pattern2))

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

        if self.filter_umi_barcodes:
            umi = self.filterUMIBarcode(umi)
            if umi is None:
                return None

        if self.filter_cell_barcodes:
            cell = self.filterCellBarcode(cell)
            if cell is None:
                return None

        if self.umi_separator:
            umi_separator = self.umi_separator

        self.read_counts['Reads output'] += 1

        # if UMI could be on either read, use umi_values to identify
        # which read(s) it was on
        if self.either_read:
            new_identifier = addBarcodesToIdentifier(
                read1, umi, cell, umi_separator)
            read1.identifier = new_identifier

            new_identifier2 = addBarcodesToIdentifier(
                read2, umi, cell, umi_separator)
            read2.identifier = new_identifier

            # UMI was on read 1
            if new_seq2 == "" and new_quals2 == "":
                read1.seq = new_seq
                read1.quals = new_quals

            # UMI was on read 2
            if new_seq == "" and new_quals == "":
                read2.seq = new_seq2
                read2.quals = new_quals2

        # Otherwise, use input from user to identiy which reads need updating
        else:
            new_identifier = addBarcodesToIdentifier(
                read1, umi, cell, umi_separator)
            read1.identifier = new_identifier
            if self.pattern:  # seq and quals need to be updated
                read1.seq = new_seq
                read1.quals = new_quals

            if read2:
                new_identifier2 = addBarcodesToIdentifier(
                    read2, umi, cell, umi_separator)
                read2.identifier = new_identifier2
                if self.pattern2:   # seq and quals need to be updated
                    read2.seq = new_seq2
                    read2.quals = new_quals2

        if read2 is None:
            return read1
        else:
            return read1, read2
