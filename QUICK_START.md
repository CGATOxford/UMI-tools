UMI-Tools quick start guide
=============================

- [UMI-Tools quick start guide](#umi-tools-quick-start-guide)
  * [Step 1: Install `UMI-Tools`](#step-1--install--umi-tools-)
  * [Step 2: Download the test data](#step-2--download-the-test-data)
  * [Step 3: Extract the UMIs](#step-3--extract-the-umis)
  * [Step 4: Mapping](#step-4--mapping)
  * [Step 5: Deduplication](#step-5--deduplication)
  * [Common variations](#common-variations)
    + [Paired-end sequencing](#paired-end-sequencing)
    + [Mapping to the transcriptome](#mapping-to-the-transcriptome)
    + [Other options](#other-options)

The following steps will guide you through a short example of how to
use the `UMI-tools` package to process data with UMIs added to them.

The general pipeline is:

extract UMI from raw reads -> map reads -> deduplicate reads based on UMIs

The most computationally intensive part of this is the middle part -
mapping the reads, is also the least interesting for us here. To aid
the ability of folks to follow along without having to worry if they
have the correct indexs install etc, we have provided a BAM file of
the mapped reads from this example. You can download it `here
<https://github.com/CGATOxford/UMI-tools/releases/download/v0.2.3/example.bam>`_. It
will need indexing with `samtools index` before use. You can then skip
striaght to Step 5 below.

Step 1: Install `UMI-Tools`
----------------------------

The easiest way to install `UMI-Tools` is using your favoirte python
package manager, either `pip` or `conda`. If you don't know which you
have install, we recommend trying both, starting with `conda`:

    $ conda install -c https://conda.anaconda.org/toms umi_tools

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


     Processed read: CCGGTTGCCCAATTGCCAAATTTTGGGGCCCCTATGAGCTA
                     ^^^^  

The processed reads could then be passed to a demultiplexing tool to
deal with the library part of the barcode.

Since the file we have downloaded contains only one library, here we
will treat the whole barcode as a UMI, and os the pattern will contain
only Ns.

    $ zcat SRR2057597.fastq.gz | umi_tools extract --bc-pattern=NNNNNNNNN --stdout processed.fastq.gz

Note that `extract` can output to a gziped or uncompressed file
depending on the file extension. It can also output to `stdout` if not
output is specified. A log file is output continaning the paramemeters
which which `extract` was run and the frequency of each UMI
encountered. boThis can be redirected with `--log` or supprssed with
`--supress-stats` (run parameters are still output).


Step 4: Mapping
---------------

The next step is to map the reads (in real life, you might also want
to demultiplex, trim and quality filter the reads). Below we will use
`bowtie` to map the reads to the mouse genome and `samtools` to create
a BAM file from the results. If you don't wish to spend the time doing
this, or don't have access to `bowtie` or `samtools` (or suitable
alternatives), we provide a premapped BAM file in the `example`
directory, or here.

First map the reads with your favoirte read mapper, here `bowtie`
using parameters from the paper which we stole the sample from. This
assumes the mm9 `bowtie` index and fasta are in your current
directory.
 
    $ bowtie --threads 4 -v 2 -m 10 -a mm9 <( gunzip < processed.fastq.gz ) --sam > mapped.sam

Next we need to convert the SAM file to BAM (actually `dedup` will use
SAM filesfor single ended analysis, but its much slower).

    $ samtools import mm9.fa mapped.sam mapped.bam  

The BAM now needs to be sorted and indexed:

    $ samtools sort mapped.bam -o example.bam
    $ samtools index example.bam

We have steps up to this point and made the result available if you
want to skip step 4. Get the file here (it will still need indexing):

    $ wget https://github.com/CGATOxford/UMI-tools/releases/download/v0.2.3/example.bam

Step 5: Deduplication
----------------------

Now that we have a mapped, sorted, indexed file, we can proceed to run
the deduplication proceedure on it:

    $ umi_tools dedup -I example.bam --output-stats=deduplicated -S deduplicated.bam

The `--output-stats` option is optional, but selecting output a range
of statistics about the run. One of the most interesting is the
distribution of edit distances. The content of this file after running
the above will look something like:


|directional-adjacency|directional-adjacency_null|edit_distance|unique|unique_null|
|---------------------|--------------------------|-------------|------|-----------|
|18313                |18313                     |Single_UMI   |15910 |15910      |
|0                    |7                         |0            |0     |4          |
|3                    |19                        |1            |1846  |45         |
|366                  |166                       |2            |976   |287        |
|711                  |547                       |3            |937   |1023       |
|1975                 |2186                      |4            |1906  |3498       |
|513                  |643                       |5            |306   |1114       |
|0                    |0                         |6            |0     |0          |

The first two columns show the distribution of average edit distances
between UMIs found at a single base in the genome after deduplication
with the directional-adjacency method (the default). Thus in the third
line we see that there are 3 bases in the genome where the average
edit distance between the UMIs found at that base is 1. The second
column is what we would expect to see if UMIs where randomly
distributed between mapping locations (taking into account any biases
in the overall usage of particular UMI sequences). The last two
columns the same, but for the naive `unique` deduplication method
where every UMI is assumed to represent an independent molecule in the
biological sample. Looking at the third row, we see that there are
1846 positions where we see that the average edit distance between
UMIs is 1, where as in the random null (in the final column) we would
only expect to see 45 such bases.

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
achieved by adding the `--paired` to the call to `dedup`:

    $ umi_tools dedup -I mapped.bam --paired -S deduplicated.bam

Paired deduplicating is signficantly slower and more memory intensive
than single-ended.

### Mapping to the transcriptome ###

A common practice in single cell RNA-seq is to map to the
transcriptome rather than the genome. The identify of the transcript
is usually indicated by the contig to which reads are mapped. In these
cases, the precise location of mapping is not informative (because
fragmentation happens after amplication), only the identify of the
gene. UMI-tools can be instructed to use this scheme using the
`--per-contig` option:

    $ umi_tools dedup -I transcriptome_mapped.bam --per-contig -S deduplicated.bam

### Other options ###

See `umi_tools extract --help` and `umi_tools dedup --help` for details of
futher possibilities. 
