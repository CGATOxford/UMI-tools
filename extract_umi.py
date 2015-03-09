'''
cgat_script_template.py - template for CGAT scripts
====================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Extract UMI barcode from a read and add it to the read name, leaving any sample
barcode in place.

Options
-------

--bc-pattern
       Use this option to specify the format of the UMI/barcode. Use Ns to
       represent the random positions and Xs to indicate the bc positions.
       Bases with Ns will be extracted and added to the read name. Remaining
       bases, marked with an X will be reattached to the read.
     
       E.g. If the pattern is NNXXNN,
       Then the read:

       @HISEQ:87:00000000 read1
       AAGGAAGCTGATTGGATGGGCTAG
       DA1AEBFGGCG01DFH00B1FF0B
       +

       will become:
       @HISEQ:87:00000000 read1 AAAA
       GGGCTGATTGGATGGGCTAG
       1AFGGCG01DFH00B1FF0B
       +
--read2-in, --read2-out
       Optionally a file with read pairs can be provided. The UMI will be
       extracted from the first read, but added to the read2 read name as well.
       Assumes that reads are in the same order. read2-out specifies the output
       file.
--3prime
       By default the barcode is assumed to be on the 5' end of the read, but
       use this option to sepecify that it is on the 3' end instead

Usage:
------

For single ended reads

      python extract_umi.py --bc-pattern=[PATTERN] [OPTIONS]

reads from stdin and outputs to stdout.

For paired end reads:

      python extract_umi.py --bc-pattern=[PATTERN] --read2-in=[FASTQIN] --read2-out=[FASTQOUT] [OPTIONS]

reads end one from stdin and end two from FASTQIN and outputs end one to stdin
and end two to FASTQOUT.

Command line options
--------------------

'''


import sys


import CGAT.Experiment as E
import CGAT.Fastq as Fastq
import CGAT.IOTools as IOTools


class Extractor:
    ''' A functor that extracts the UMI and the barcode from the read,
    adds it to the identifier and pastes back on the sample barcode'''

    bc_count = {}

    def _extract_5prime(self, sequence):
        return (sequence[:self.pattern_length], sequence[self.pattern_length:])

    def _extract_3prime(self, sequence):
        return (sequence[-self.pattern_length:],
                sequence[:-self.pattern_length])

    def _joiner_5prime(self, sequence, sample):
        return  sample + sequence
        
    def _joiner_3prime(self, read, sequence, sample):
        return sequence+sample

    def __init__(self, pattern, prime3=False):

        if prime3:
            self.extract = self._extract_3prime
            self.joiner = self._joiner_3prime
        else:
            self.extract = self._extract_5prime
            self.joiner = self._joiner_5prime

        self.pattern_length = len(pattern)
        self.umi_bases = [x for x in range(len(pattern)) if pattern[x] is "N"]
        self.bc_bases = [x for x in range(len(pattern)) if pattern[x] is "X"]

    def __call__(self, read1, read2=None):

        bc, sequence = self.extract(read1.seq)
        bc_qual, seq_qual = self.extract(read1.quals)

        umi = "".join([bc[x] for x in self.umi_bases])
        sample = "".join([bc[x] for x in self.bc_bases])
        sample_qual = "".join([bc_qual[x] for x in self.bc_bases])

        id = (bc, umi, sample)
        self.bc_count[id] = self.bc_count.get(id, 0) + 1

        read1.seq = self.joiner(sequence, sample)
        read1.quals = self.joiner(seq_qual, sample_qual)
        read1.identifier = read1.identifier.split("/")[0] + " " + umi

        if not read2 is None:
            read2.identifier = read2.identifier.split("/")[0] + " " + umi
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
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-p", "--bc-pattern", dest="pattern", type="string",
                      help="Barcode pattern. Ns are random bases X's fixed")
    parser.add_option("--read2-in", dest="read2_in", type="string",
                      help="file name for read pairs")
    parser.add_option("--3prime", dest="prime3", action="store_true",
                      help="barcode is on 3' end of read")
    parser.add_option("--read2-out", dest="read2_out", type="string",
                      help="file to output processed paired read to")
    parser.add_option("--supress-stats", dest="stats", action="store_false",
                      help="Suppress the writing of stats to the log")

    parser.set_defaults(pattern=None,
                        read2_in=None,
                        read2_out=None,
                        prime3=False,
                        stats=True)

    # add common options (-h/--help, ...) and parse command line

    (options, args) = E.Start(parser, argv=argv)

    #Initialise the processor
    processor = Extractor(options.pattern, options.prime3)
    read1s = Fastq.iterate(options.stdin)

    if options.read2_in is None:

        for read in read1s:
            options.stdout.write(str(processor(read)) + "\n")

    else:
        
        read2s = Fastq.iterate(IOTools.openFile(options.read2_in))
        read2_out = IOTools.openFile(options.read2_out)

        for read1, read2 in zip(read1s, read2s):
            new_1, new_2 = processor(read1, read2)
            options.stdout.write(str(new_1) + "\n")
            read2_out.write(str(new_2) + "\n")

    # write footer and output benchmark information.

    if options.stats:
     
        options.stdlog.write("\t".join(["Barcode", "UMI", "Sample", "Count"]) + "\n")
        for id in processor.bc_count:
            options.stdlog.write("\t".join(id+(str(processor.bc_count[id]),)) + "\n")
                             
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
