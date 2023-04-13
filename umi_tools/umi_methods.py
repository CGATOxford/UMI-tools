'''
umi_methods.py - Methods for dealing with UMIs
=========================================================

'''

from __future__ import absolute_import
import itertools
import collections
import pysam
import numpy as np

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


def fastqIterate(infile, remove_suffix=False):
    '''iterate over contents of fastq file. If remove_suffic is True then /1 
    and /2 is removed from the first field of the identifier. Raises ValueError
    if remove_suffix is true and /1 or /2 is not present at the end of the
    first field of any reads'''

    def convert2string(b):
        if type(b) == str:
            return b
        else:
            return b.decode("utf-8")
 
    def removeReadIDSuffix(line):
        '''Removes /1 or /2 from the a provided identifier :param:line.
        Identifier must be in the first field of the identifier (with ' ' as
        the field seperator) and sepreated with a '/'. Raises ValueError if this
        is not the case'''
        
        components = line.split(' ')
        read_id = components[0]
        
        if not read_id[-2:] in ['/1', '/2']:
            raise ValueError(
                'read suffix must be /1 or /2. Observed: %s' % read_id[-2:])

        read_id = read_id[:-2]
        components[0] = read_id
        line = ' '.join(components)
        
        return(line)
     
    while 1:
        line1 = convert2string(infile.readline()).strip()
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

        if remove_suffix:
            line1 = removeReadIDSuffix(line1)
            
        yield Record(line1[1:], line2[:-1], line4[:-1])

# End of FastqIterate()
###############################################################################
###############################################################################


def joinedFastqIterate(fastq_iterator1, fastq_iterator2,
                       strict=True):
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

    If provided iterators were created with `remove_suffix=True`, /1 and /2 will
    be removed from the end of read1 and read2, respectively before
    checking their names are identical
    '''

    def getReadID(read):
        return(read.identifier.split(' ')[0])

    for read1 in fastq_iterator1:
        read2 = next(fastq_iterator2)

        read1_id = getReadID(read1)
        read2_id = getReadID(read2)

        if not strict:
            while read2_id != read1_id:
                read2 = next(fastq_iterator2)
                read2_id = getReadID(read2)

        if not read2_id == read1_id:
            raise ValueError("\nRead pairs do not match\n%s != %s" %
                             (read1_id, read2_id))

        yield (read1, read2)


def get_average_umi_distance(umis):

    if len(umis) == 1:
        return -1

    dists = [edit_distance(x, y) for
             x, y in itertools.combinations(umis, 2)]
    return float(sum(dists))/(len(dists))


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
        self.first_kerror = 1
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

            try:
                self.umis[self.barcode_getter(read)[0]] += 1
            except KeyError:
                continue

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
