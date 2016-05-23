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


--read2-in, --read2-out
       Optionally a file with read pairs can be provided. The UMI will
       be added to the read2 read name as well.  Assumes that reads
       are in the same order. read2-out specifies the output file.

--split-barcode
       By default the UMI is assumed to be on the first read. Use this
       option if the UMI is contained on both reads and specify the
       pattern of the barcode/UMI on the second read using the option
       ``--bc-pattern2``

--bc-pattern
       Use this option to specify the format of the UMI/barcode. Use Ns to
       represent the random positions and Xs to indicate the bc positions.
       Bases with Ns will be extracted and added to the read name. Remaining
       bases, marked with an X will be reattached to the read.

       E.g. If the pattern is NNXXNN,
       Then the read:

       @HISEQ:87:00000000 read1
       AAGGTTGCTGATTGGATGGGCTAG
       DA1AEBFGGCG01DFH00B1FF0B
       +

       will become:
       @HISEQ:87:00000000_AATT read1
       GGGCTGATTGGATGGGCTAG
       1AFGGCG01DFH00B1FF0B
       +

--bc-pattern2
       Use this option to specify the format of the UMI/barcode for
       the second read pair if required. If --bc-pattern2 is not
       supplied, this defaults to the same pattern as --bc-pattern

--3prime
       By default the barcode is assumed to be on the 5' end of the read, but
       use this option to sepecify that it is on the 3' end instead

-L
       Specify a log file to retain logging information and final statistics

--supress-stats
Supress logging of summary statistics

Usage:
------

For single ended reads:
        python extract_umi.py --bc-pattern=[PATTERN] -L extract.log [OPTIONS]

reads from stdin and outputs to stdout.

For paired end reads:
        python extract_umi.py --bc-pattern=[PATTERN] --read2-in=[FASTQIN] --read2-out=[FASTQOUT] -L extract.log [OPTIONS]

reads end one from stdin and end two from FASTQIN and outputs end one to stdin
and end two to FASTQOUT.

Command line options
--------------------

'''
import sys

# python 3 doesn't require izip
try:
    # Python 2
    from itertools import izip
except ImportError:
    # Python 3
    izip = zip

import umi_tools.Utilities as U


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
       ``illumina-1.8``, ``solexa`` or ``phred64``.

    """
    def __init__(self, identifier, seq, quals, format=None):
        self.identifier, self.seq, self.quals, format = (
            identifier, seq, quals, format)
        self.format = None

    def __str__(self):
        return "@%s\n%s\n+\n%s" % (self.identifier, self.seq, self.quals)


def fastqIterate(infile):
    '''iterate over contents of fastq file.'''

    while 1:
        line1 = infile.readline()
        if not line1:
            break
        if not line1.startswith('@'):
            raise ValueError("parsing error: expected '@' in line %s" % line1)
        line2 = infile.readline()
        line3 = infile.readline()
        if not line3.startswith('+'):
            raise ValueError("parsing error: expected '+' in line %s" % line3)
        line4 = infile.readline()
        # incomplete entry
        if not line4:
            raise ValueError("incomplete entry for %s" % line1)

        yield Record(line1[1:-1], line2[:-1], line4[:-1])

# End of FastqIterate()
###############################################################################
###############################################################################


def addUMItoIdentifier(read, UMI):
    '''extract the identifier from a read and append the UMI before
    the first space'''

    read_id = read.identifier.split(" ")
    read_id[0] = read_id[0] + "_" + UMI
    identifier = " ".join(read_id)

    return identifier


class Extractor:
    ''' A functor that extracts the UMI and the barcode from the read(s),
    adds it to the identifier and pastes back on the sample barcode'''

    bc_count = {}

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

    def __init__(self, pattern, pattern2, prime3=False):

        if prime3:
            self.extract = self._extract_3prime
            self.joiner = self._joiner_3prime
        else:
            self.extract = self._extract_5prime
            self.joiner = self._joiner_5prime

        self.pattern_length = len(pattern)
        self.umi_bases = [x for x in range(len(pattern)) if pattern[x] is "N"]
        self.bc_bases = [x for x in range(len(pattern)) if pattern[x] is "X"]
        self.split = False

        if pattern2:
            self.pattern_length2 = len(pattern2)
            self.split = True
            self.umi_bases2 = [x for x in range(len(pattern2))
                               if pattern2[x] is "N"]
            self.bc_bases2 = [x for x in range(len(pattern2))
                              if pattern2[x] is "X"]

    def __call__(self, read1, read2=None):

        bc1, sequence1 = self.extract(read1.seq)
        bc_qual1, seq_qual1 = self.extract(read1.quals)

        umi1 = "".join([bc1[x] for x in self.umi_bases])
        sample1 = "".join([bc1[x] for x in self.bc_bases])
        sample_qual1 = "".join([bc_qual1[x] for x in self.bc_bases])

        if self.split:
            bc2, sequence2 = self.extract(read2.seq, read=2)
            bc_qual2, seq_qual2 = self.extract(read2.quals, read=2)

            umi2 = "".join([bc2[x] for x in self.umi_bases2])
            sample2 = "".join([bc2[x] for x in self.bc_bases2])
            sample_qual2 = "".join([bc_qual2[x] for x in self.bc_bases2])

            umi = umi1 + umi2
            bc = "_".join((bc1, bc2))
            sample = "_".join((sample1, sample2))

            read2.seq = self.joiner(sequence2, sample2)
            read2.quals = self.joiner(seq_qual2, sample_qual2)

        else:
            umi = umi1
            bc = bc1
            sample = sample1

        id = (bc, umi, sample)
        self.bc_count[id] = self.bc_count.get(id, 0) + 1

        read1.seq = self.joiner(sequence1, sample1)
        read1.quals = self.joiner(seq_qual1, sample_qual1)
        read1.identifier = addUMItoIdentifier(read1, umi)

        if read2 is not None:
            read2.identifier = addUMItoIdentifier(read2, umi)

            return (read1, read2)
        else:
            return read1


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = U.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("--split-barcode", dest="split", action="store_true",
                      help="barcode is split across read pair")
    parser.add_option("-p", "--bc-pattern", dest="pattern", type="string",
                      help="Barcode pattern. Ns are random bases X's fixed")
    parser.add_option("--bc-pattern2", dest="pattern2", type="string",
                      help="Barcode pattern. Ns are random bases X's fixed")
    parser.add_option("--read2-in", dest="read2_in", type="string",
                      help="file name for read pairs")
    parser.add_option("--3prime", dest="prime3", action="store_true",
                      help="barcode is on 3' end of read")
    parser.add_option("--read2-out", dest="read2_out", type="string",
                      help="file to output processed paired read to")
    parser.add_option("--supress-stats", dest="stats", action="store_false",
                      help="Suppress the writing of stats to the log")

    parser.set_defaults(split=False,
                        pattern=None,
                        pattern2=None,
                        read2_in=None,
                        read2_out=None,
                        prime3=False,
                        stats=True)

    # add common options (-h/--help, ...) and parse command line

    (options, args) = U.Start(parser, argv=argv)

    # check options
    if not options.pattern:
        raise ValueError("must specify a pattern using ``--bc-pattern``")

    if options.split:
        if not options.read2_in:
            raise ValueError("must specify a paired fastq ``--read2-in``")

        if not options.pattern2:
            options.pattern2 = options.pattern

    if options.read2_in:
        if not options.read2_out:
            raise ValueError("must specify an output for the paired end "
                             "``--read2-out``")

    # Initialise the processor
    processor = Extractor(options.pattern, options.pattern2, options.prime3)
    read1s = fastqIterate(options.stdin)

    if options.read2_in is None:

        for read in read1s:
            options.stdout.write(str(processor(read)) + "\n")

    else:

        read2s = fastqIterate(U.openFile(options.read2_in))
        read2_out = U.openFile(options.read2_out, "w")

        for read1, read2 in izip(read1s, read2s):
            new_1, new_2 = processor(read1, read2)
            options.stdout.write(str(new_1) + "\n")
            read2_out.write(str(new_2) + "\n")

    # write footer and output benchmark information.

    if options.stats:

        options.stdlog.write("\t".join(["Barcode", "UMI", "Sample", "Count"]) + "\n")
        for id in processor.bc_count:
            options.stdlog.write("\t".join(id+(str(processor.bc_count[id]),)) + "\n")

    U.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
