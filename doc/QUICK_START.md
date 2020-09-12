Quick start guide
=============================

This quick start guide uses an iCLIP dataset as an example. If you want to apply UMI-tools in a single cell RNA-Seq analysis, please see the [Single cell tutorial](Single_cell_tutorial.md)

- [UMI-Tools quick start guide](#umi-tools-quick-start-guide)
  * [Step 1: Install `UMI-Tools`](#step-1--install--umi-tools-)
  * [Step 2: Download the test data](#step-2--download-the-test-data)
  * [Step 3: Extract the UMIs](#step-3--extract-the-umis)
  * [Step 4: Mapping](#step-4--mapping)
  * [Step 5: Deduplication](#step-5--deduplication)
  * [Common variations](#common-variations)
    + [Paired-end sequencing](#paired-end-sequencing)
    + [Read grouping](#read-grouping)
    + [Other options](#other-options)

The following steps will guide you through a short example of how to
use the `UMI-tools` package to process data with UMIs added to them. The
data used comes from one of the control replicates from [Mueller-Mcnicoll et al, 
*Genes Dev* (2016) 30: 553](http://genesdev.cshlp.org/content/30/5/553). We
have adaptor trimmed and filtered the data to reduce its size.

The general pipeline is:

extract UMI from raw reads -> map reads -> deduplicate reads based on UMIs

The most computationally intensive part of this is the middle part -
mapping the reads. It is also the least interesting for us here. To aid
the ability of folks to follow along without having to worry if they
have the correct indexs install etc, we have provided a BAM file of
the mapped reads from this example. You can download it [here](https://github.com/CGATOxford/UMI-tools/releases/download/v0.2.3/example.bam>). It
will need indexing with `samtools index` before use. You can then skip
straight to Step 5 below.

Step 1: Install `UMI-Tools`
----------------------------

The easiest way to install `UMI-Tools` is using your favorite python
package manager, either `pip` or `conda`. If you don't know which you
have installed, we recommend trying both, starting with `conda`:

    $ conda install -c bioconda umi_tools

Alternatively, `pip`:

    $ pip install umi_tools

If neither of these work, see our installation guide.

Step 2: Download the test data
--------------------------------

The test data we are going to use is a control sample from a recent
iCLIP experiment. We have processed this data to remove the various
adaptors that you find in iCLIP data. You can download the trimmed
file here:

    $ wget https://github.com/CGATOxford/UMI-tools/releases/download/v0.2.3/example.fastq.gz

If you're using macOS use:

    $ curl -L "https://github.com/CGATOxford/UMI-tools/releases/download/v0.2.3/example.fastq.gz" -o "example.fastq.gz"

The file is about 100Mb, and takes a couple of minutes to download on
our system.

Step 3: Extract the UMIs
-------------------------

UMIs are strings of random nucleotides attached to the start of
reads. Before the reads are mapped the random nucleotides must be
removed from the reads, but the sequence must be kept. The `extract`
command of `UMI-Tools` moves the UMI from the read to the read name.

Several techniques that use UMIs mix the UMI sequence in with a
library barcode. In this case we want to remove the random part of the
barcodes, but leave the library part so that the reads can be
de-multiplexed. We specify this using the `--bc-pattern` parameter to
`extract`. Ns represent the random part of the barcode and Xs the
fixed part. For example, in a standard iCLIP experiment, the barcode
is made of 3 random bases, followed by a 4 base library barcode,
followed by 2 more random bases. Thus the `--bc-pattern` would be
"NNNXXXXNN".

      Read:          TAGCCGGCTTTGCCCAATTGCCAAATTTTGGGGCCCCTATGAGCTAG 
      Barcode:       NNNXXXXNN
                         |
                         v
                     TAGCCGGCT
                         |
                         V
           random-> TAG CCGG CT <- random
                         ^
                         |
                      library


     Processed read: CCGGTTGCCCAATTGCCAAATTTTGGGGCCCCTATGAGCTAG
                     ^^^^  

The processed reads could then be passed to a demultiplexing tool to
deal with the library part of the barcode.

Since the file we have downloaded contains only one library, here we
will treat the whole barcode as a UMI, and so the pattern will contain
only Ns.

    $ umi_tools extract --stdin=example.fastq.gz --bc-pattern=NNNNNNNNN --log=processed.log --stdout processed.fastq.gz 

Note that `extract` can output to a gzipped or uncompressed file
depending on the file extension. It can also output to `stdout` if no
output is specified. A log file is saved containing the parameters
with which `extract` was run and the frequency of each UMI
encountered. This can be redirected with `--log` or suppressed with
`--supress-stats` (run parameters are still output).


Step 4: Mapping
---------------

The next step is to map the reads (in real life, you might also want
to demultiplex, trim and quality filter the reads). Below we will use
`bowtie` to map the reads to the mouse genome and `samtools` to create
a BAM file from the results. If you don't wish to spend the time doing
this, or don't have access to `bowtie` or `samtools` (or suitable
alternatives), we provide a premapped BAM file (see command at the end of this step).

First map the reads with your favorite read mapper, here `bowtie`,
using parameters from the paper which we stole the sample from. This
assumes the mm9 `bowtie` index and fasta are in your current
directory.
 
    $ bowtie --threads 4 -v 2 -m 10 -a mm9 <( gunzip < processed.fastq.gz ) --sam > mapped.sam

Next we need to convert the SAM file to BAM (actually `dedup` will use
SAM files for single ended analysis, but it's much slower).

    $ samtools import mm9.fa mapped.sam mapped.bam  

The BAM now needs to be sorted and indexed:

    $ samtools sort mapped.bam -o example.bam
    $ samtools index example.bam

If you want to skip the mapping, you can get the file here.
It will still need indexing (see "samtools index" command above):

    $ wget https://github.com/CGATOxford/UMI-tools/releases/download/v0.2.3/example.bam

Again for macOS use the following to download:

    $ curl -L "https://github.com/CGATOxford/UMI-tools/releases/download/v0.2.3/example.bam" -o "example.bam"
    
Step 5: Deduplication
----------------------

Now that we have a mapped, sorted, indexed file, we can proceed to run
the deduplication procedure on it:

    $ umi_tools dedup -I example.bam --output-stats=deduplicated -S deduplicated.bam

The `--output-stats` option is optional, but selecting it will provide a range
of statistics about the run. One of the most interesting is the
distribution of edit distances (here named deduplicated_edit_distance.tsv).
The content of this file after running the above will look something like:

|directional-adjacency|directional-adjacency_null|edit_distance|unique|unique_null|
|---------------------|--------------------------|-------------|------|-----------|
|10465                |10465                     |Single_UMI   |8976  |8976       |
|0                    |1                         |0            |0     |2          |
|2                    |8                         |1            |1167  |33         |
|164                  |73                        |2            |496   |183        |
|211                  |258                       |3            |281   |566        |
|700                  |916                       |4            |730   |1663       |
|395                  |320                       |5            |317   |617        |
|89                   |4                         |6            |59    |5          |
|21                   |2                         |7            |21    |2          |
|0                    |0                         |8            |0     |0          |




The first two columns show the distribution of average edit distances
between UMIs found at a single base in the genome after deduplication
with the directional-adjacency method (the default). Thus in the third
line we see that there are 2 bases in the genome where the average
edit distance between the UMIs found at that base is 1. The second
column is what we would expect to see if UMIs were randomly
distributed between mapping locations (taking into account any biases
in the overall usage of particular UMI sequences). The last two
columns the same, but for the naive `unique` deduplication method
where every UMI is assumed to represent an independent molecule in the
biological sample. Looking at the third row, we see that there are
1167 positions where the average edit distance between
UMIs is 1, whereas in the random null (in the final column) we would
only expect to see 33 such bases.

The statistics options signficantly reduce the speed at which
deduplication is performed and increase the memory usage. If time or
memory usage is an issue, try running without the `--output-stats`
option.

Common variations
------------------


### Paired-end sequencing ###

If paired-end sequencing has been performed, it is necessary to make
sure that the UMI sequence is added to both reads. When processing,
provide the second read like so:

    $ umi_tools extract -I pair.1.fastq.gz --bc-pattern=NNNXXXXNN \ 
      --read2-in=pair.2.fastq.gz --stdout=processed.1.fastq.gz \
      --read2-out=processed.2.fastq.gz

This assumes the UMI is on the 5' end of read1. For other
possibilities (such as a UMI on both reads) see `umi_tools extract
--help`. After paired-end mapping, paired end deduplication can be
achieved by adding the `--paired` option to the call to `dedup`:

    $ umi_tools dedup -I mapped.bam --paired -S deduplicated.bam

Paired deduplicating is signficantly slower and more memory intensive
than single-ended.

### Read grouping ###
For some applications it may be neccessary to mark the duplicates but retain all reads, for example,
where the PCR duplicates are used to correct sequence errors by generating a consensus sequence. In
these cases, the _group_ command can be used to mark each read with its read group. Optionally a flatfile
detailing the read groups and read identifiers can also be output using the `--group-out` option:

    $ umi_tools group -I mapped.bam --paired --group-out=groups.tsv --output-bam -S mapped_grouped.bam
    
The output bam will contain two tags: UG = read group id, BX = read group UMI.
The tag containing the read group UMI can be modified with the `--umi-group-tag` option.

The groups flatfile contains the following columns:
- read_id
- contig
- position
- umi = raw umi
- umi_count = how many times was this umi observed at the same alignment coordinates
- final_umi = the error corrected umi
- final_umi_count = how many times was the umi observed at the same alignment coordinates, inc. error correction
- unique_id = the unique identifier for this group
    

#### Example UMI extraction: ####

In the case above the UMIs are extracted with the pattern --bc-pattern=NNXXXXNN. Below is an example of how the fastq should be formatted following extraction:

UMI is bases 3-7, bases 1-2 and 7-8 are the sample barcode and need to be removed

    @HISEQ:87:00000000 read1
    AAGGTTGCTGATTGGATGGGCTAG
    +
    DA1AEBFGGCG01DFH00B1FF0B
 
 will become:
 
     @HISEQ:87:00000000_GGTT read1
     TGATTGGATGGGCTAG
     + 
     1AFGGCG01DFH00B1

### Other options ###

See `umi_tools extract --help` and `umi_tools dedup --help` for details of
futher possibilities. 
