'''
==============================================================================
prepare_for_rsem - make the output from dedup or group compatible with RSEM
===============================================================================
The SAM format specification states that the mnext and mpos fields should point
to the primary alignment of a read's mate. However, not all aligners adhere to
this standard. In addition, the RSEM software requires that the mate of a read1
appears directly after it in its input BAM. This requires that there is exactly
one read1 alignment for every read2 and vice versa.

In general (except in a few edge cases) UMI tools outputs only the read2 to that 
corresponds to the read specified in the mnext and mpos positions of a selected
read1, and only outputs this read once, even if multiple read1s point to it. 
This makes UMI-tools outputs incompatible with RSEM. This script takes the output
from dedup or groups and ensures that each read1 has exactly one read2 (and vice
versa), that read2 always appears directly after read1,and that pairs point to 
each other (note this is technically not valid SAM format). Copy any specified
tags from read1 to read2 if they are present (by default, UG and BX, the unique
group and correct UMI tags added by _group_)

Input must to name sorted.

'''


from umi_tools import Utilities as U
from collections import defaultdict
import pysam
import sys

usage = '''
prepare_for_rsem - make output from dedup or group compatible with RSEM

Usage: umi_tools prepare_for_rsem [OPTIONS] [--stdin=IN_BAM] [--stdout=OUT_BAM]

       note: If --stdout is ommited, standard out is output. To
             generate a valid BAM file on standard out, please
             redirect log with --log=LOGFILE or --log2stderr '''

def chunk_bam(bamfile):
    '''Take in a iterator of pysam.AlignmentSegment entries and yeild
    lists of reads that all share the same name'''

    last_query_name = None
    output_buffer = list()

    for read in bamfile:
    
        if last_query_name is not None and last_query_name != read.query_name:
            yield(output_buffer)
            output_buffer = list()
        
        output_buffer.append(read)

    yield (output_buffer)

def copy_tags(tags, read1, read2):
    '''Given a  list of tags, copies the values of these tags from read1
    to read2, if the tag is set'''

    for tag in tags:
        
        try:
            read1_tag = read1.get_tag(tag, with_value_type=True)
            read2.set_tag(tag, value=read1_tag[0], value_type=read1_tag[1])
        except KeyError:
            pass
        
    return(read2)


def main(argv=None):

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = U.OptionParser(version="%prog version: $Id$",
                            usage=usage,
                            description=globals()["__doc__"])
    group = U.OptionGroup(parser, "RSEM preparation specific options")

    group.add_option("--tags", dest="tags", type="string",
                     default="UG,BX",
                     help="Comma-seperated list of tags to transfer from read1 to read2")
    group.add_option("--sam", dest="sam", action="store_true",
                     default=False,
                     help="input and output SAM rather than BAM")

    parser.add_option_group(group)

    # add common options (-h/--help, ...) and parse command line
    (options, args) = U.Start(parser, argv=argv,
                              add_group_dedup_options=False,
                              add_umi_grouping_options=False,
                              add_sam_options=False)

    if options.stdin != sys.stdin:
        in_name = options.stdin.name
        options.stdin.close()
    else:
        in_name = "-"

    if options.sam:
        mode = ""
    else:
        mode = "b"

    inbam = pysam.AlignmentFile(in_name, "r"+mode)

    if options.stdout != sys.stdout:
        out_name = options.stdout.name
        options.stdout.close()
    else:
        out_name = "-"

    outbam = pysam.AlignmentFile(out_name, "w" + mode, template = inbam)

    options.tags = options.tags.split(",")

    for template in chunk_bam(inbam):
        
        current_template = {True: defaultdict(list),
                            False: defaultdict(list)}

        for read in template:
            key = (read.reference_name, read.pos, not read.is_secondary)
            current_template[read.is_read1][key].append(read)

        output = set()

        for read in template:
           
            mate = None
           
            # if this read is a non_primary alignment, we first want to check if it has a mate
            # with the non-primary alignment flag set. 
            mate_key_secondary = (read.next_reference_name, read.next_reference_start, False)
            if read.is_secondary:

                # get a list of secondary reads at the correct alignment position
                potential_secondary_mates = current_template[not read.is_read1][mate_key_secondary]

                # search through one at a time to find a read that points to the current read
                # as its mate.
                for candidate_mate in potential_secondary_mates:
                    if candidate_mate.next_reference_name == read.reference_name and \
                       candidate_mate.next_reference_start == read.pos:
                        mate = candidate_mate
                
                # if no such read is found, then pick any old secondary alignment at that position
                # note: this happens when UMI-tools outputs the wrong read as somethings pair.
                if mate is None and len(potential_secondary_mates) >0:
                    mate = potential_secondary_mates[0]

            # The following happens if read is primary, or is secondary, but doesn't point to a secondary
            # alignment as its mate (or at least we didn't find one). If which case find the primary 
            # alignment at that location.
            if mate is None:
                mate_key_primary = (read.next_reference_name, read.next_reference_start, True)

                try:
                    # there should be exactly one primary alignment at this location
                    mate = current_template[not read.is_read1][mate_key_primary][0]
                except (KeyError, IndexError):
                    U.warn("Alignment {} has no mate -- skipped".format(
                        "\t".join([read.query_name, read.flag, read.reference_name, read.pos])
                    ))
                    continue
            
            # because we might want to make changes to the read, but not have those changes reflected
            # if we need the read again,we copy the read. This is only way I can find to do this.
            read = pysam.AlignedSegment().from_dict(read.to_dict(), read.header)
            mate = pysam.AlignedSegment().from_dict(mate.to_dict(), read.header)

            # Make it so that if our read is secondary, the mate is also secondary. We don't make the
            # mate primary if the read is primary because we would otherwise end up with mulitple
            # primary alignments. 
            if read.is_secondary:
                mate.is_secondary = True
            
            # In a situation where there is already one mate for each read, then we will come across
            # each pair twice - once when we scan read1 and once when we scan read2. Thus we need
            # to make sure we don't output something already output. 
            if read.is_read1:
                
                mate = copy_tags(options.tags, read, mate)  
                output_key = str(read) + str(mate)

                if output_key not in output:
                    output.add(output_key)
                    outbam.write(read)
                    outbam.write(mate)

            elif read.is_read2:

                read = copy_tags(options.tags, mate, read)
                output_key = str(mate) + str(read)

                if output_key not in output:
                    output.add(output_key)
                    outbam.write(mate)
                    outbam.write(read)

            else:
                U.warn("Alignment {} is neither read1 nor read2 -- skipped".format(
                    "\t".join([read.query_name, read.flag, read.reference_name, read.pos])
                ))
                continue

    if not out_name == "-":
        outbam.close()
    U.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))



