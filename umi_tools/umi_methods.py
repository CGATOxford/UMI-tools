'''
umi_methods.py - Methods for dealing with UMIs
=========================================================

:Author: Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Python UMI

'''

import itertools
import collections
import random
import numpy as np
import pysam

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


def get_umi_read_id(read, sep="_"):
    ''' extract the umi from the read id using the specified separator '''

    try:
        return read.qname.split(sep)[-1].encode('utf-8')
    except IndexError:
        raise ValueError(
            "Could not extract UMI from the read ID, please"
            "check UMI is encoded in the read name")


def get_umi_tag(read, tag='RX'):
    ''' extract the umi from the specified tag '''

    try:
        # 10X pipelines append a 'GEM' tag to the UMI, e.g
        # AGAGSGATAGATA-1
        return read.get_tag(tag).split("-")[0].encode('utf-8')
    except IndexError:
        raise ValueError(
            "Could not extract UMI from the read tags, please"
            "check UMI is encoded in the read tag: %s" % tag)


def get_average_umi_distance(umis):

    if len(umis) == 1:
        return -1

    dists = [edit_distance(x, y) for
             x, y in itertools.combinations(umis, 2)]
    return float(sum(dists))/(len(dists))


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
        if tags:
            self.read2tags = {}
        else:
            self.read2tags = None

    def write(self, read, unique_id=None, umi=None, unmapped=False):
        '''Check if chromosome has changed since last time. If it has, scan
        for mates. Write the read to outfile and save the identity for paired
        end retrieval'''

        if unmapped:
            self.outfile.write(read)
            return

        if not self.chrom == read.reference_id:
            self.write_mates()
            self.chrom = read.reference_id

        key = read.query_name, read.next_reference_id, read.next_reference_start
        self.read1s.add(key)

        if self.read2tags is not None:
            self.read2tags[key] = (unique_id, umi)
            read.tags += [('UG', unique_id)]
            read.tags += [('FU', umi)]

        self.outfile.write(read)

    def write_mates(self):
        '''Scan the current chromosome for matches to any of the reads stored
        in the read1s buffer'''

        if self.chrom is not None:
            U.debug("Dumping %i mates for contig %s" % (
                len(self.read1s), self.infile.get_reference_name(self.chrom)))

        for read in self.infile.fetch(tid=self.chrom, multiple_iterators=True):
            if any((read.is_unmapped, read.mate_is_unmapped, read.is_read1)):
                continue

            key = read.query_name, read.reference_id, read.reference_start
            if key in self.read1s:
                if self.read2tags is not None:
                    unique_id, umi = self.read2tags[key]
                    self.read2tags.pop(key)

                    read.tags += [('UG', unique_id)]
                    read.tags += [('FU', umi)]

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
        for name, chrom, pos in self.read1s:
            for read in self.infile.fetch(start=pos, end=pos+1, tid=chrom):
                if (read.query_name, read.pos) == (name, pos):
                    if self.read2tags:
                        unique_id, umi = self.read2tags[name, chrom, pos]
                        read.tags += [('UG', unique_id)]
                        read.tags += [('FU', umi)]

                    self.outfile.write(read)
                    found += 1
                    break

        U.info("%i mates never found" % (len(self.read1s) - found))
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


def getMetaContig2contig(gene_transcript_map):
    ''' '''
    metacontig2contig = collections.defaultdict(set)
    for line in U.openFile(gene_transcript_map, "r"):

        if line.startswith("#"):
            continue

        if len(line.strip()) == 0:
            break

        gene, transcript = line.strip().split("\t")
        metacontig2contig[gene].add(transcript)

    return metacontig2contig


def metafetcher(bamfile, metacontig2contig):
    ''' return reads in order of metacontigs'''
    for metacontig in metacontig2contig:
        for contig in metacontig2contig[metacontig]:
            for read in bamfile.fetch(contig):
                read.tags += [('MC', metacontig)]
                yield read


def get_bundles(inreads, ignore_umi=False, subset=None,
                quality_threshold=0, paired=False, spliced=False,
                soft_clip_threshold=0, per_contig=False,
                per_gene=False, gene_tag=None,
                whole_contig=False, read_length=False,
                detection_method="MAPQ", umi_getter=None,
                all_reads=False, return_unmapped=False):
    ''' Returns a dictionary of dictionaries, representing the unique reads at
    a position/spliced/strand combination. The key to the top level dictionary
    is a umi. Each dictionary contains a "read" entry with the best read, and a
    count entry with the number of reads with that position/spliced/strand/umi
    combination'''

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
            if read.is_unmapped and return_unmapped:
                yield read, read_events, 'unmapped'
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
                yield read, read_events, 'unmapped'
            continue

        if read.mate_is_unmapped and paired:
            if not read.is_unmapped:
                read_events['Read 2 unmapped'] += 1
            if return_unmapped:
                yield read, read_events, 'unmapped'
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

        if per_contig or per_gene or gene_tag:

            if per_contig:
                pos = read.tid
                key = pos
            elif per_gene:
                pos = read.get_tag('MC')
                key = pos
            elif gene_tag:
                pos = read.get_tag(gene_tag)
                key = pos

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
            umi = umi_getter(read)

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


class random_read_generator:
    ''' class to generate umis at random based on the
    distributon of umis in a bamfile '''

    def __init__(self, bamfile, chrom, umi_getter):
        inbam = pysam.Samfile(bamfile)

        if chrom:
            self.inbam = inbam.fetch(reference=chrom)
        else:
            self.inbam = inbam.fetch()

        self.umis = collections.defaultdict(int)
        self.umi_getter = umi_getter
        self.fill()

    def fill(self):

        self.frequency2umis = collections.defaultdict(list)

        for read in self.inbam:

            if read.is_unmapped:
                continue

            if read.is_read2:
                continue

            self.umis[self.umi_getter(read)] += 1

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
