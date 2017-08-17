'''
count_tab.py - Count reads per gene from flatfile using UMIs
=================================================================

:Author: Ian Sudbery, Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Python UMI

Purpose
-------

The purpose of this command is to count the number of reads per gene
based on the read's gene assignment and UMI. Note this command is not
currently able to perform per-cell counting. See the count command if
you want to perform per-cell counting.

The input must be in the following format (tab separated), where the
first column is the read identifier and the second column is the
assigned gene. The input must be sorted by the gene identifier:

NS500668:144:H5FCJBGXY:2:22309:18356:15843_TCTAA    ENSG00000279457.3
NS500668:144:H5FCJBGXY:3:23405:3971:19716_CGATG     ENSG00000225972.1

You can perform any required file transformation and pipe the output
directly to count_tab. For example to pipe output from featureCounts
with the '-R' option you can do the following:

    awk '$2=="Assigned" {print $1"\t"$3}' my.bam.featureCounts| sort -k2 |
    umi_tools count_tab -S gene_counts.tsv -L count.log

The tab file is assumed to contain each read id once only. For paired
end reads with featureCounts you must include the "-p" option so each
read id is included once only.

'''

import sys

# required to make iteritems python2 and python3 compatible
from builtins import dict

from functools import partial

import umi_tools.Utilities as U
import umi_tools.network as network
import umi_tools.umi_methods as umi_methods

# add the generic docstring text
__doc__ = __doc__ + U.GENERIC_DOCSTRING_GDC


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = U.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    # add common options (-h/--help, ...) and parse command line
    (options, args) = U.Start(parser, argv=argv, add_group_dedup_options=False)

    nInput, nOutput = 0, 0

    # set the method with which to extract umis from reads
    umi_getter = partial(
        umi_methods.get_umi_read_string, sep=options.umi_sep)

    options.stdout.write("%s\t%s\n" % ("gene", "count"))

    # set up UMIClusterer functor with methods specific to
    # specified options.method
    processor = network.UMIClusterer(options.method)

    for gene, counts in umi_methods.get_gene_count_tab(
            options.stdin,
            umi_getter=umi_getter):

        umis = counts.keys()

        nInput += sum(counts.values())

        # group the umis
        groups = processor(
            umis,
            counts,
            threshold=options.threshold)

        gene_count = len(groups)
        options.stdout.write("%s\t%i\n" % (gene, gene_count))
        nOutput += gene_count

    U.info("Number of reads counted: %i" % nOutput)

    U.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
