'''
count_tab - Count reads per gene from flatfile using UMIs
=================================================================

Purpose
-------

The purpose of this command is to count the number of reads per gene
based on the read's gene assignment and UMI. See the count command if
you want to perform per-cell counting using a BAM file input.

The input must be in the following format (tab separated), where the
first column is the read identifier (including UMI) and the second
column is the assigned gene. The input must be sorted by the gene
identifier.

Input template::

    read_id[SEP]_UMI    gene

Example::

    NS500668:144:H5FCJBGXY:2:22309:18356:15843_TCTAA     ENSG00000279457.3
    NS500668:144:H5FCJBGXY:3:23405:39715:19716_CGATG     ENSG00000225972.1

You can perform any required file transformation and pipe the output
directly to count_tab. For example to pipe output from featureCounts
with the '-R CORE' option you can do the following::

    awk '$2=="Assigned" {print $1"\t"$4}' my.bam.featureCounts | sort -k2 |
    umi_tools count_tab -S gene_counts.tsv -L count.log

The tab file is assumed to contain each read id once only. For paired
end reads with featureCounts you must include the "-p" option so each
read id is included once only.

Per-cell counting can be enable with ``--per-cell``. For per-cell
counting, the input must be in the following format (tab separated),
where the first column is the read identifier (including UMI and Cell
Barcode) and the second column is the assigned gene. The input must be
sorted by the gene identifier:

Input template::

    read_id[SEP]_UMI_CB    gene

Example::

    NS500668:144:H5FCJBGXY:2:22309:18356:15843_TCTAA_AGTCGA     ENSG00000279457.3
    NS500668:144:H5FCJBGXY:3:23405:39715:19716_CGATG_GGAGAA     ENSG00000225972.1

'''

import sys

# required to make iteritems python2 and python3 compatible
from builtins import dict

from functools import partial

import umi_tools.Utilities as U
import umi_tools.Documentation as Documentation
import umi_tools.network as network
import umi_tools.umi_methods as umi_methods
import umi_tools.sam_methods as sam_methods

# add the generic docstring text
__doc__ = __doc__ + Documentation.GENERIC_DOCSTRING_GDC

usage = '''
count_tab - Count reads per gene from flatfile using UMIs

Usage: umi_tools count_tab [OPTIONS] [--stdin=IN_TSV[.gz]] [--stdout=OUT_TSV[.gz]]

       note: If --stdin/--stdout are ommited standard in and standard
             out are used for input and output. Input/Output will be
             (de)compressed if a filename provided to --stdin/--stdout
             ends in .gz '''


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = U.OptionParser(version="%prog version: $Id$",
                            usage=usage,
                            description=globals()["__doc__"])

    if len(argv) == 1:
        parser.print_usage()
        print ("Required options missing, see --help for more details")
        return 1

    group = U.OptionGroup(parser, "count_tab-specific options")

    group.add_option("--barcode-separator", dest="bc_sep",
                     type="string", help="separator between read id and UMI "
                     " and (optionally) the cell barcode", default="_")

    group.add_option("--per-cell", dest="per_cell",
                     action="store_true",
                     help="Readname includes cell barcode as well as UMI in "
                     "format: read[sep]UMI[sep]CB")

    parser.add_option_group(group)

    # add common options (-h/--help, ...) and parse command line
    (options, args) = U.Start(parser, argv=argv, add_group_dedup_options=False,
                              add_sam_options=False)

    nInput, nOutput = 0, 0

    # set the method with which to extract umis from reads
    if options.per_cell:
        bc_getter = partial(
            sam_methods.get_cell_umi_read_string, sep=options.bc_sep)
    else:
        bc_getter = partial(
            sam_methods.get_umi_read_string, sep=options.bc_sep)

    if options.per_cell:
        options.stdout.write("%s\t%s\t%s\n" % ("cell", "gene", "count"))
    else:
        options.stdout.write("%s\t%s\n" % ("gene", "count"))

    # set up UMIClusterer functor with methods specific to
    # specified options.method
    processor = network.UMIClusterer(options.method)

    for gene, counts in sam_methods.get_gene_count_tab(
            options.stdin,
            bc_getter=bc_getter):

        for cell in counts.keys():
            umis = counts[cell].keys()

            nInput += sum(counts[cell].values())

            # group the umis
            groups = processor(
                counts[cell],
                threshold=options.threshold)

            gene_count = len(groups)
            if options.per_cell:
                options.stdout.write("%s\t%s\t%i\n" % (cell, gene, gene_count))
            else:
                options.stdout.write("%s\t%i\n" % (gene, gene_count))
                nOutput += gene_count

    U.info("Number of reads counted: %i" % nOutput)

    U.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
