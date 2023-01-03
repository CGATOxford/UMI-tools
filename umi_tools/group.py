'''
==============================================================
Group - Group reads based on their UMI and mapping coordinates
==============================================================

*Identify groups of reads based on their genomic coordinate and UMI*

The group command can be used to create two types of outfile: a tagged
BAM or a flatfile describing the read groups

To generate the tagged-BAM file, use the option ``--output-bam`` and
provide a filename with the ``--stdout``/``-S`` option. Alternatively,
if you do not provide a filename, the bam file will be outputted to
the stdout. If you have provided the ``--log``/``-L`` option to send
the logging output elsewhere, you can pipe the output from the group
command directly to e.g samtools view like so::

    umi_tools group -I inf.bam --group-out=grouped.tsv --output-bam
    --log=group.log --paired | samtools view - |less

The tagged-BAM file will have two tagged per read:

 - UG
   Unique_id. 0-indexed unique id number for each group of reads
   with the same genomic position and UMI or UMIs inferred to be
   from the same true UMI + errors
 - BX
   Final UMI. The inferred true UMI for the group

To generate the flatfile describing the read groups, include the
``--group-out=<filename>`` option. The columns of the read groups file are
below. The first five columns relate to the read. The final 3 columns
relate to the group.

  - read_id
      read identifier

  - contig
      alignment contig

  - position
      Alignment position. Note that this position is not the start
      position of the read in the BAM file but the start of the read
      taking into account the read strand and cigar

  - gene
      The gene assignment for the read. Note, this will be NA unless the
      --per-gene option is specified

  - umi
      The read UMI

  - umi_count
      The number of times this UMI is observed for reads at the same
      position

  - final_umi
      The inferred true UMI for the group

  - final_umi_count
      The total number of reads within the group

  - unique_id
      The unique id for the group


group-specific options
----------------------

"""""""""""
--group-out
"""""""""""
   Outfile name for file mapping read id to read group

"""""""""
--out-bam
"""""""""
   Output a bam file with read groups tagged using the UG tag

"""""""""""""""
--umi-group-tag
"""""""""""""""
   BAM tag for the error corrected UMI selected for the group. Default=BX


'''
import sys
import collections
import os

# required to make iteritems python2 and python3 compatible
from builtins import dict
from future.utils import iteritems

import pysam

import umi_tools.Utilities as U
import umi_tools.Documentation as Documentation
import umi_tools.network as network
import umi_tools.sam_methods as sam_methods

# add the generic docstring text
__doc__ = __doc__ + Documentation.GENERIC_DOCSTRING_GDC
__doc__ = __doc__ + Documentation.GROUP_DEDUP_GENERIC_OPTIONS

usage = '''
group - Group reads based on their UMI

Usage: umi_tools group --output-bam [OPTIONS] [--stdin=INFILE.bam] [--stdout=OUTFILE.bam]

       note: If --stdout is ommited, standard out is output. To
             generate a valid BAM file on standard out, please
             redirect log with --log=LOGFILE or --log2stderr '''


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

    group = U.OptionGroup(parser, "group-specific options")

    group.add_option("--group-out", dest="tsv", type="string",
                     help="Outfile name for file mapping read id to read group",
                     default=None)

    group.add_option("--output-bam", dest="output_bam", action="store_true",
                     default=False,
                     help=("output a bam file with read groups tagged using the UG tag"
                           "[default=%default]"))

    parser.add_option("--umi-group-tag", dest="umi_group_tag",
                      type="string", help="tag for the outputted umi group",
                      default='BX')

    parser.add_option_group(group)

    # add common options (-h/--help, ...) and parse command line
    (options, args) = U.Start(parser, argv=argv)

    U.validateSamOptions(options, group=True)

    if options.stdin != sys.stdin:
        in_name = options.stdin.name
        options.stdin.close()
    else:
        raise ValueError("Input on standard in not currently supported")

    if options.stdout != sys.stdout:
        if options.no_sort_output:
            out_name = options.stdout.name
        else:
            out_name = U.getTempFilename(dir=options.tmpdir)
            sorted_out_name = options.stdout.name
        options.stdout.close()
        assert options.output_bam, (
            "To output a bam you must include --output-bam option")
    else:
        if options.no_sort_output:
            out_name = "-"
        else:
            out_name = U.getTempFilename(dir=options.tmpdir)
            sorted_out_name = "-"

    if not options.no_sort_output:  # need to determine the output format for sort
        if options.out_sam:
            sort_format = "sam"
        else:
            sort_format = "bam"

    if options.in_sam:
        in_mode = "r"
    else:
        in_mode = "rb"

    if options.out_sam:
        out_mode = "wh"
    else:
        out_mode = "wb"

    infile = pysam.Samfile(in_name, in_mode)

    if options.output_bam:
        outfile = pysam.Samfile(out_name, out_mode, template=infile)
    else:
        outfile = None

    if options.tsv:
        mapping_outfile = U.openFile(options.tsv, "w")
        mapping_outfile.write("%s\n" % "\t".join(
            ["read_id", "contig", "position", "gene", "umi", "umi_count",
             "final_umi", "final_umi_count", "unique_id"]))

    nInput, nOutput, unique_id, input_reads, output_reads = 0, 0, 0, 0, 0

    gene_tag = options.gene_tag
    metacontig2contig = None

    if options.unmapped_reads in ["use", "output"]:
        output_unmapped = True
    else:
        output_unmapped = False

    if options.chrom:
        inreads = infile.fetch(reference=options.chrom)
    else:
        if options.per_gene and options.gene_transcript_map:
            metacontig2contig = sam_methods.getMetaContig2contig(
                infile, options.gene_transcript_map)
            metatag = "MC"
            inreads = sam_methods.metafetcher(infile, metacontig2contig, metatag)
            gene_tag = metatag

        else:
            inreads = infile.fetch(until_eof=output_unmapped)

    bundle_iterator = sam_methods.get_bundles(
        options,
        all_reads=True,
        return_read2=True,
        return_unmapped=output_unmapped,
        metacontig_contig=metacontig2contig)

    # set up UMIClusterer functor with methods specific to
    # specified options.method
    processor = network.UMIClusterer(options.method)

    for bundle, key, status in bundle_iterator(inreads):

        # write out read2s and unmapped/chimeric (if these options are set)
        if status == 'single_read':
            # bundle is just a single read here
            nInput += 1

            if outfile:
                outfile.write(bundle)

            nOutput += 1
            continue

        umis = bundle.keys()
        counts = {umi: bundle[umi]["count"] for umi in umis}

        nInput += sum(counts.values())

        while nOutput >= output_reads + 10000:
            output_reads += 10000
            U.info("Written out %i reads" % output_reads)

        while nInput >= input_reads + 1000000:
            input_reads += 1000000
            U.info("Parsed %i input reads" % input_reads)

        # group the umis
        groups = processor(
            counts,
            threshold=options.threshold)

        for umi_group in groups:
            top_umi = umi_group[0]

            group_count = sum(counts[umi] for umi in umi_group)

            for umi in umi_group:
                reads = bundle[umi]['read']
                for read in reads:
                    if outfile:
                        # Add the 'UG' tag to the read
                        read.set_tag('UG', unique_id)
                        read.set_tag(options.umi_group_tag, top_umi)
                        outfile.write(read)

                    if options.tsv:
                        if options.per_gene:
                            if options.per_contig:
                                gene = read.reference_name
                            else:
                                gene = read.get_tag(gene_tag)
                        else:
                            gene = "NA"
                        mapping_outfile.write("%s\n" % "\t".join(map(str, (
                            read.query_name, read.reference_name,
                            sam_methods.get_read_position(
                                read, options.soft_clip_threshold)[1],
                            gene,
                            umi.decode(),
                            counts[umi],
                            top_umi.decode(),
                            group_count,
                            unique_id))))

                    nOutput += 1

            unique_id += 1

    if outfile:
        outfile.close()
        if not options.no_sort_output:
            # sort the output
            pysam.sort("-o", sorted_out_name, "-O", sort_format, "--no-PG", out_name)
            os.unlink(out_name)  # delete the tempfile

    if options.tsv:
        mapping_outfile.close()

    # write footer and output benchmark information.
    U.info(
        "Reads: %s" % ", ".join(["%s: %s" % (x[0], x[1]) for x in
                                 bundle_iterator.read_events.most_common()]))
    U.info("Number of reads out: %i, Number of groups: %i" %
           (nOutput, unique_id))

    U.info("Total number of positions deduplicated: %i" %
           processor.positions)
    if processor.positions > 0:
        U.info("Mean number of unique UMIs per position: %.2f" %
               (float(processor.total_umis_per_position) /
                processor.positions))
        U.info("Max. number of unique UMIs per position: %i" %
               processor.max_umis_per_position)
    else:
        U.warn("The BAM did not contain any valid "
               "reads/read pairs for deduplication")

    U.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
