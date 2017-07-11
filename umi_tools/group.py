'''
group.py - Group reads based on their UMI
=========================================

:Author: Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Python UMI

Purpose
-------

The purpose of this command is to identify groups of reads based on
their genomic coordinate and UMI. It is assumed that the FASTQ files
were processed with extract_umi.py before mapping and thus the UMI is
the last word of the read name. e.g:

@HISEQ:87:00000000_AATT

where AATT is the UMI sequeuence.

If you have used an alternative method which does not separate the
read id and UMI with a "_", such as bcl2fastq which uses ":", you can
specify the separator with the option "--umi-separator=<sep>",
replacing <sep> with e.g ":".

Alternatively, if your UMIs are encoded in a tag, you can specify this
by setting the option --extract-umi-method=tag and set the tag name
with the --umi-tag option. For example, if your UMIs are encoded in
the 'UM' tag, provide the following options:
"--extract-umi-method=tag --umi-tag=UM"

Finally, if you have used umis to extract the UMI +/- cell barcode,
you can specify --extract-umi-method=umis

By default, reads are considered identical if they have the same start
coordinate, are on the same strand, and have the same UMI. Optionally,
splicing status can be considered (see below).

The start postion of a read is considered to be the start of its alignment
minus any soft clipped bases. A read aligned at position 500 with
cigar 2S98M will be assumed to start at postion 498.

Methods
-------

group can be run with multiple methods to identify group of reads with
the same (or similar) UMI(s). All methods start by identifying the
reads with the same mapping position.

The simpliest method, "unique", groups reads with the exact same
UMI. The network-based methods, "cluster", "adjacency" and
"directional", build networks where nodes are UMIs and edges connect
UMIs with an edit distance <= threshold (usually 1). The groups of
reads are then defined from the network in a method-specific manner.

Note that the "percentile" method used with the group command is not
available with group. This is because this method does not group
similar UMIs as per the network methods. Instead it applies a
threshold for inclusion of the UMI in the output and excluded UMIs are
not assigned to a "true" UMI.

  "unique"
      Reads group share the exact same UMI

  "cluster"
      Identify clusters of connected UMIs (based on hamming distance
      threshold). Each network is a read group

  "adjacency"
      Cluster UMIs as above. For each cluster, select the node(UMI)
      with the highest counts. Visit all nodes one edge away. If all
      nodes have been visted, stop. Otherise, repeat with remaining
      nodes until all nodes have been visted. Each step
      defines a read group.

  "directional"
      Identify clusters of connected UMIs (based on hamming distance
      threshold) and umi A counts >= (2* umi B counts) - 1. Each
      network is a read group.

The group command can be used to create two types of outfile: a tagged
BAM or a flatfile describing the read groups

To generate the tagged-BAM file, use the option --output-bam and
provide a filename with the -S option. Alternatively, if you do not
provide a filename, the bam file will be outputted to the stdout. If
you have provided the --log/-L option to send the logging output
elsewhere, you can pipe the output from the group command directly to
e.g samtools view like so:

umi_tools group -I inf.bam --group-out=grouped.tsv --output-bam
--log=group.log --paired | samtools view - |less

The tagged-BAM file will have two tagged per read:
UG = Unique_id. 0-indexed unique id number for each group of reads
     with the same genomic position and UMI or UMIs inferred to be
     from the same true UMI + errors
BX = Final UMI. The inferred true UMI for the group

To generate the flatfile describing the read groups, include the
--group-out=<filename> option. The columns of the read groups file are
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


Options
-------

--extract-umi-method (choice)
      How are the UMIs encoded in the read?

      Options are:

      - "read_id" (default)
            UMIs contained at the end of the read separated as
            specified with --umi-separator option

      - "tag"
            UMIs contained in a tag, see --umi-tag option

--umi-separator (string)
      Separator between read id and UMI. See --extract-umi-method above

--umi-tag (string)
      Tag which contains UMI. See --extract-umi-method above

--method (choice, string)
      Method used to identify PCR duplicates within reads. All methods
      start by identifying the reads with the same mapping position

      Options are:

      - "unique"

      - "cluster"

      - "adjacency"

      - "directional" (default)

--edit-distance-threshold (int)
       For the adjacency and cluster methods the threshold for the
       edit distance to connect two UMIs in the network can be
       increased. The default value of 1 works best unless the UMI is
       very long (>14bp)

--no-sort-output
       By default, output from UMI-tools are sorted. This involves the
       use of a temporary unsorted file since reads are considered in
       the order of their start position which is may not be the same
       as their alignment coordinate due to soft-clipping and reverse
       alignments. The temp file will be saved in $TMPDIR and deleted
       when it has been sorted to the outfile. Use this option to turn
       off sorting.

--paired
       BAM is paired end - output both read pairs. This will also
       force the use of the template length to determine reads with
       the same mapping coordinates.

--spliced-is-unique
       Causes two reads that start in the same position on the same
       strand and having the same UMI to be considered unique if one is spliced
       and the other is not. (Uses the 'N' cigar operation to test for
       splicing)

--soft-clip-threshold (int)
       Mappers that soft clip, will sometimes do so rather than mapping a
       spliced read if there is only a small overhang over the exon
       junction. By setting this option, you can treat reads with at least
       this many bases soft-clipped at the 3' end as spliced.

--read-length
      Use the read length as as a criteria when grouping, for e.g sRNA-Seq

--whole-contig (string)
      forces group to parse an entire contig before yielding any reads
      for grouping. This is the only way to absolutely guarantee
      that all reads with the same start position are grouped together
      for grouping since group uses the start position of the
      read, not the alignment coordinate on which the reads are
      sorted. However, by default, group reads for another 1000bp
      before outputting read groups which will avoid any reads being
      missed with short read sequencing (<1000bp)

--subset (float, [0-1])
      Only consider a fraction of the reads, chosen at random. This is useful
      for doing saturation analyses.

--chrom
      Only consider a single chromosome. This is useful for debugging purposes

--per-contig (string)
      Group per contig (field 3 in BAM; RNAME).

--per-gene (string)
      Group per gene. As above except with this option you can
      align to a reference transcriptome with more than one transcript
      per gene. You need to also provide --gene-transcript-map option.
      This will also add a metacontig ('MC') tag to the reads if used
      in conjunction with --output-bam

--gene-transcript-map (string)
      File mapping genes to transripts (tab separated), e.g:

      gene1   transcript1
      gene1   transcript2
      gene2   transcript3

--gene-tag (string)
      Group per gene. As per --per-gene except here the gene
      information is encoded in the bam read tag specified so you do
      not need to supply --gene-transcript-map

--group-out (string, filename)
      Output a flatfile describing the read groups

--output-bam (string, filename)
      Output a tagged bam file to stdout or -S <filename>

-i, --in-sam/-o, --out-sam
      By default, inputs are assumed to be in BAM format and output are output
      in BAM format. Use these options to specify the use of SAM format for
      inputs or outputs.

-I    (string, filename) input file name
      The input file must be sorted and indexed.

-S    (string, filename) output file name

-L    (string, filename) log file name


.. note::
   In order to get a valid sam/bam file you need to redirect logging
   information or turn it off logging via -v 0. You can redirect the
   logging to a file with -L <logfile> or use the --log2stderr option
   to send the logging to stderr.

'''
import sys
import collections
import os

# required to make iteritems python2 and python3 compatible
from builtins import dict
from future.utils import iteritems

import pysam
import numpy as np

try:
    import umi_tools.Utilities as U
except ImportError:
    import Utilities as U

try:
    import umi_tools.network as network
except ImportError:
    import network

try:
    import umi_tools.umi_methods as umi_methods
except ImportError:
    import umi_methods


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = U.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    group = U.OptionGroup(parser, "group-specific options")

    group.add_option("--group-out", dest="tsv", type="string",
                     help="Outfile name for file mapping read id to read group",
                     default=None)

    group.add_option("--output-bam", dest="output_bam", action="store_true",
                     default=False,
                     help=("output a bam file with read groups tagged using the UG tag"
                           "[default=%default]"))

    group.add_option("--output-unmapped", dest="output_unmapped", action="store_true",
                     default=False,
                     help=("Retain all unmapped reads in output[default=%default]"))

    parser.add_option("--umi-group-tag", dest="umi_group_tag",
                      type="string", help="tag for the outputted umi group",
                      default='BX')

    parser.add_option_group(group)

    # add common options (-h/--help, ...) and parse command line
    (options, args) = U.Start(parser, argv=argv)

    U.validateSamOptions(options)

    if options.stdin != sys.stdin:
        in_name = options.stdin.name
        options.stdin.close()
    else:
        raise ValueError("Input on standard in not currently supported")

    if options.stdout != sys.stdout:
        if options.no_sort_output:
            out_name = options.stdout.name
        else:
            out_name = U.getTempFilename()
            sorted_out_name = options.stdout.name
        options.stdout.close()
        assert options.output_bam, (
            "To output a bam you must include --output-bam option")
    else:
        if options.no_sort_output:
            out_name = "-"
        else:
            out_name = U.getTempFilename()
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

    if options.chrom:
        inreads = infile.fetch(reference=options.chrom)
    else:
        if options.per_gene and options.gene_transcript_map:
            metacontig2contig = umi_methods.getMetaContig2contig(
                infile, options.gene_transcript_map)
            metatag = "MC"
            inreads = umi_methods.metafetcher(infile, metacontig2contig, metatag)
            gene_tag = metatag

        else:
            inreads = infile.fetch(until_eof=options.output_unmapped)

    bundle_iterator = umi_methods.get_bundles(
        options,
        all_reads=True,
        return_read2=True,
        metacontig_contig=metacontig2contig)

    for bundle, key, status in bundle_iterator(inreads):

        # write out read2s and unmapped (if these options are set)
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
            U.debug("Outputted %i" % nOutput)

        while nInput >= output_reads + 1000000:
            input_reads += 1000000
            U.debug("Read %i input reads" % nInput)

        # set up UMIClusterer functor with methods specific to
        # specified options.method
        processor = network.UMIClusterer(options.method)

        # group the umis
        groups = processor(
            umis,
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
                        read.tags += [('UG', unique_id)]
                        read.tags += [(options.umi_group_tag, top_umi)]
                        outfile.write(read)

                    if options.tsv:
                        if options.per_gene:
                            gene = read.get_tag(gene_tag)
                        else:
                            gene = "NA"
                        mapping_outfile.write("%s\n" % "\t".join(map(str, (
                            read.query_name, read.reference_name,
                            umi_methods.get_read_position(read, options.soft)[1],
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
            pysam.sort("-o", sorted_out_name, "-O", sort_format, out_name)
            os.unlink(out_name)  # delete the tempfile

    if options.tsv:
        mapping_outfile.close()

    # write footer and output benchmark information.
    U.info(
        "Reads: %s" % ", ".join(["%s: %s" % (x[0], x[1]) for x in
                                 bundle_iterator.read_events.most_common()]))
    U.info("Number of reads out: %i, Number of groups: %i" %
           (nOutput, unique_id))
    U.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
