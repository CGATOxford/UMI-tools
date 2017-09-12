'''
whitelist_methods.py - Methods for whitelisting cell barcodes
=============================================================

:Author: Tom Smith, Ian Sudbery
:Release: $Id$
:Date: |today|
:Tags: Python UMI

'''
import collections
import re
import random
from functools import partial
import umi_tools.Utilities as U


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


def get_umi_read_string(read_id, sep="_"):
    ''' extract the umi from the read id (input as a string) using the
    specified separator '''

    try:
        return (None, read_id.split(sep)[-1].encode('utf-8'))
    except IndexError:
        raise ValueError(
            "Could not extract UMI from the read ID, please"
            "check UMI is encoded in the read name")


def get_cell_umi_read_string(read_id, sep="_"):
    ''' extract the umi and cell barcode from the read id (input as a
    string) using the specified separator '''

    try:
        return (read_id.split(sep)[-1].encode('utf-8'),
                read_id.split(sep)[-2].encode('utf-8'))
    except IndexError:
        raise ValueError(
            "Could not extract UMI or CB from the read ID, please"
            "check UMI and CB are encoded in the read name:"
            "%s" % read_id)


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


class get_readpairs:
    # Used in place of a real read as the mate of reads that aren't properly paired
    NoMate = object()

    def incomplete_pair_to_singleton(self, queue_entry):
        for i in [i for i, v in enumerate(queue_entry) if v is get_readpairs.NoMate]:
            queue_entry[i] = None

    def incomplete_pairs_as_singletons(self):
        self.incomplete_queue_entries.clear()
        while self.queue:
            queue_entry = self.queue.popleft()
            if get_readpairs.NoMate in queue_entry:
                # Incomplete entry. Turn into singleton entry
                read = next(r for r in queue_entry if r is not get_readpairs.NoMate)
                U.warn("inconsistent BAM file: only one of the two mates of proper pair %s found on chromosome %s, treating as unpaired" %
                       (read.qname, read.reference_name))
                incomplete_pair_to_singleton(queue_entry)
            yield queue_entry

    def __call__(self, inreads):
        self.queue = collections.deque()
        self.incomplete_queue_entries = dict()
        self.current_chr = None

        for read in inreads:
            read_i = int(read.is_read2)
            read_chr = read.reference_name

            # Output leftover incomplete read pairs if the chrosomome changes
            if self.current_chr is not None and self.current_chr != read_chr:
                for queue_entry in self.incomplete_pairs_as_singletons():
                    yield queue_entry
            self.current_chr = read_chr

            # Queue read, either as a singleton, or as part of a read pair
            if read.is_unmapped or read.mate_is_unmapped or not read.is_proper_pair:
                # Read is not part of a proper pair. Set mate to None
                queue_entry = [None, None]
                queue_entry[read_i] = read
                self.queue.append(queue_entry)
            elif read.qname in self.incomplete_queue_entries:
                # Read is part of an already queued (incomplete) read pair. Complete it.
                queue_entry = self.incomplete_queue_entries.pop(read.qname)
                if queue_entry[read_i] is not get_readpairs.NoMate:
                    U.warn("inconsistent BAM file: both mates %s flagged to be read%d" %
                           (read.qname, read_i+1))
                    read_i = 1 - read_i
                queue_entry[read_i] = read
            else:
                # Read is part of a new read pair. Create & queue incomplete pair.
                queue_entry = [get_readpairs.NoMate, get_readpairs.NoMate]
                queue_entry[read_i] = read
                self.queue.append(queue_entry)
                self.incomplete_queue_entries[read.qname] = queue_entry

            # Output queued reads and readpairs up to the first incomplete pair
            while self.queue and sum(1 for r in self.queue[0] if r is get_readpairs.NoMate) == 0:
                yield self.queue.popleft()

        # Output leftover incomplete read pairs at the end of the input
        for queue_entry in self.incomplete_pairs_as_singletons():
            yield queue_entry


class get_singletons:
    def __call__(self, inreads):
        for read in inreads:
            if read.is_read2:
                yield None, read
            else:
                yield read, None


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

    def update_dicts(self, read, read2, pos, key, umi):

        # The content of the reads_dict depends on whether all reads
        # are being retained

        if self.all_reads:
            # retain all reads per key
            try:
                self.reads_dict[pos][key][umi]["count"] += 1
            except KeyError:
                self.reads_dict[pos][key][umi]["read"] = [read]
                self.reads_dict[pos][key][umi]["read2"] = [read2]
                self.reads_dict[pos][key][umi]["count"] = 1
            else:
                self.reads_dict[pos][key][umi]["read"].append(read)
                self.reads_dict[pos][key][umi]["read2"].append(read2)

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
                self.reads_dict[pos][key][umi]["read2"] = read2
                self.reads_dict[pos][key][umi]["count"] = 1
                self.read_counts[pos][key][umi] = 0
            else:
                dict_r1 = self.reads_dict[pos][key][umi]["read"]
                dict_r2 = self.reads_dict[pos][key][umi]["read2"]
                dict_mapq = dict_r1.mapq + (dict_r2.mapq if dict_r2 is not None else 0)
                read_mapq = read.mapq + (read2.mapq if read2 is not None else 0)
                if dict_mapq > read_mapq:
                    return

                if read_mapq > dict_mapq:
                    self.reads_dict[pos][key][umi]["read"] = read
                    self.reads_dict[pos][key][umi]["read2"] = read2
                    self.read_counts[pos][key][umi] = 0
                    return

                # TS: implemented different checks for multimapping here
                # XXX: What to do in paired mode? Probably just complain...
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
                    self.reads_dict[pos][key][umi]["read2"] = read2

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
        if self.options.paired:
            generator = get_readpairs()
        else:
            generator = get_singletons()

        for read, read2 in generator(inreads):

            # Unpaired second reads are returned if the return_read2 option is set
            if read is None:
                if self.return_read2:
                    if not read2.is_unmapped or (
                            read2.is_unmapped and self.return_unmapped):
                        yield read2, None, "single_read"
                continue
            else:
                self.read_events['Input Reads'] += 1

            # only ever dealing with read1s from here

            if self.options.paired:
                if read.is_paired:
                    self.read_events['Read pairs'] += 1
                else:
                    self.read_events['Unpaired reads'] += 1

                    # if paired end input and read1 is unpaired...

                    # skip, or
                    if self.options.unpaired_reads == "discard":
                        continue

                    # yield without grouping, or
                    elif self.options.unpaired_reads == "output":
                        yield read, None, "single_read"

                    # Use read pair; TLEN will be 0
                    elif self.options.unpaired_reads == "use":
                        pass

            if read.is_unmapped:
                if self.options.paired:
                    if read.mate_is_unmapped:
                        self.read_events['Both unmapped'] += 1
                    else:
                        self.read_events['Read 1 unmapped'] += 1
                else:
                    self.read_events['Single end unmapped'] += 1

                # if read1 is unmapped, yield immediately or skip read
                if self.return_unmapped:
                    yield read, None, "single_read"
                continue

            if self.options.paired and read.mate_is_unmapped:
                self.read_events['Read 2 unmapped'] += 1

                # if paired end input and read2 is unmapped, skip unless
                # options.unmapped_reads == "use", in which case TLEN will be 0
                if self.options.unmapped_reads != "use":
                    if self.return_unmapped:
                        yield read, None, "single_read"
                        continue

            if read.is_paired and (
                    read.reference_name != read.next_reference_name):
                self.read_events['Chimeric read pair'] += 1

                # if paired end input and read2 is mapped to another contig...

                # skip, or
                if self.options.chimeric_pairs == "discard":
                    continue

                # yield without grouping, or
                elif self.options.chimeric_pairs == "output":
                    yield read, None, "single_read"
                    continue

                # Use read pair; TLEN will be 0
                elif self.options.chimeric_pairs == "use":
                    pass

            if self.options.subset:
                if random.random() >= self.options.subset:
                    self.read_events['Randomly excluded'] += 1
                    continue

            if self.options.mapping_quality:
                if read.mapq < self.options.mapping_quality:
                    self.read_events['< MAPQ threshold'] += 1
                    continue

            # get the umi +/- cell barcodes
            if self.options.ignore_umi:
                if self.options.per_cell:
                    umi, cell = self.barcode_getter(read)
                    umi = ""
                else:
                    umi, cell = "", ""
            else:
                try:
                    umi, cell = self.barcode_getter(read)
                except KeyError:
                    error_msg = "Read skipped, missing umi and/or cell tag"
                    if self.read_events[error_msg] == 0:

                        # pysam renamed .tostring -> to_string in 0.14
                        # .tostring requies access to the parent AlignmentFile
                        try:
                            formatted_read = read.to_string()
                        except AttributeError:
                            formatted_read = read.query_name

                        U.warn("At least one read is missing UMI and/or "
                               "cell tag(s): %s" % formatted_read)
                    self.read_events[error_msg] += 1
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
                        assigned = read.get_tag(self.options.assigned_tag)
                        gene = read.get_tag(self.options.gene_tag)
                    except KeyError:
                        self.read_events['Read skipped, no tag'] += 1
                        continue

                    if gene == "":
                        if self.read_events['Read skipped - gene string is empty'] == 0:
                            U.warn("Assigned gene is empty string. First such "
                                   "read:\n%s" % read.to_string())
                        self.read_events['Read skipped - gene string is empty'] += 1
                        continue

                    if re.search(self.options.skip_regex, assigned):
                        self.read_events['Read skipped - assigned tag matches skip_regex'] += 1
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

                key = (read.is_reverse, self.options.spliced and is_spliced,
                       self.options.paired*read.tlen, r_length)

            # update dictionaries
            key = (key, cell)
            self.update_dicts(read, read2, pos, key, umi)

            if self.metacontig_contig:
                # keep track of observed contigs for each gene
                self.observed_contigs[gene].add(transcript)

        # yield remaining bundles
        for p in sorted(self.reads_dict.keys()):
            for k in sorted(self.reads_dict[p].keys()):
                yield self.reads_dict[p][k], k, "bundle"


def get_gene_count_tab(infile,
                       bc_getter=None):

    ''' Yields the counts per umi for each gene

    bc_getter: method to get umi (plus optionally, cell barcode) from
    read, e.g get_umi_read_id or get_umi_tag


    TODO: ADD FOLLOWING OPTION

    skip_regex: skip genes matching this regex. Useful to ignore
                unassigned reads (as per get_bundles class above)

    '''

    gene = None
    counts = collections.Counter()

    for line in infile:

        values = line.strip().split("\t")

        assert len(values) == 2, "line: %s does not contain 2 columns" % line

        read_id, assigned_gene = values

        if assigned_gene != gene:
            if gene:
                yield gene, counts

            gene = assigned_gene
            counts = collections.defaultdict(collections.Counter)

        cell, umi = bc_getter(read_id)
        counts[cell][umi] += 1

    # yield final values
    yield gene, counts


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

            if any((read.is_unmapped, read.mate_is_unmapped, read.is_read1)):
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
                read.set_tag(metatag, metacontig)
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
