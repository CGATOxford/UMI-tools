UMI-Tools quick start guide
=============================

The following steps will guide you through a short example of how to use the `UMI-tools` package to process data with UMIs added to them.

The general pipeline is:

extract UMI from raw reads -> map reads -> deduplicate reads based on UMIs

The most computationally intensive part of this is the middle part - mapping the reads, is also the least interesting for us here. To aid the ability of folks to follow along without having to worry if they have the correct indexs install etc, we have provided a BAM file of the mapped reads from this example. They are located in the `example` folder, part of the `UMI-Tools` download, alternatively, if you used `pip` or `conda` to install `UMI-Tools` and don't know where to find the files, you can download it `here <>`_.

Step 1: Install `UMI-Tools`
----------------------------

The easiest way to install `UMI-Tools` is using your favoirte python package manager, either `pip` or `conda`. If you don't know which you have install, we recommend trying both, starting with `conda`:


    $ conda install -c https://conda.anaconda.org/toms umi_tools

Alternatively, `pip`:

    $ pip install umi_tools

If neither of these work, see our installation guide.

Step 2: Download the test data
--------------------------------

The test data we are going to use is a control sample from a recent iCLIP experiment. Fetch it like to:

    $ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR205/007/SRR2057597/SRR2057597.fastq.gz

The file is about 200Mb, and takes a couple of minutes to download on our system. 

Step 3: Extract the UMIs
-------------------------

UMIs are strings of random nucleotides attached to the start of reads. Before the reads are mapped the random nucleotides must be removed from the reads, but the sequence must be kept. The `extract` command of `UMI-Tools` moves the UMI from the read to the read name.

Several techniques that use UMIs mix the UMI sequence in with a library barcode. In this case we want to remove the random part of the barcodes, but leave the library part so that the reads can be de-multiplexed. We specify this using the `--bc-pattern` parameter to `extract`. Ns represent the random part of the barcode and Xs the fixed part. For example, in a standard iCLIP experiment, the barcode is made of 3 random bases, followed by a 4 base library barcode, followed by 2 more random bases. Thus the `--bc-pattern` would be "NNNXXXXNN".

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

The processed reads could then be passed to a demultiplexing tool to deal with the library part of the barcode.

Since the file we have downloaded contains only one library, here we will treat the whole barcode as a UMI, and os the pattern will contain only Ns.

   $ zcat SRR2057597.fastq.gz | umi_tools extract --bc-pattern=NNNNNNNNN --stdout proccessed.fastq.gz

Note that `extract` can output to a gziped or uncompressed file depending on the file extension. It can also output to `stdout` if not output is specified. A log file is output continaning the paramemeters which which `extract` was run and the frequency of each UMI encountered. This can be redirected with `--log` or supprssed with `--supress-stats` (run parameters are still output). 

Step 4: Mapping
---------------

The next step is to map the reads (in real life, you might also want to demultiplex, trim and quality filter the reads). Below we will use `bowtie` to map the reads to the mouse genome and `samtools` to create a BAM file from the results. If you don't wish to spend the time doing this, or don't have access to `bowtie` or `samtools` (or suitable alternatives), we provide a premapped BAM file in the `example` directory, or here.

First map the reads with your favoirte read mapper, here `bowtie` using parameters from the paper which we stole the sample from. This assumes the mm9 `bowtie` index and fasta are in your current directory.
 
    $ bowtie --threads 4 -v 2 -m 10 -a mm9 <( gunzip < processed.fastq.gz ) --sam > mapped.sam

Next we need to convert the SAM file to BAM (actually `dedup` will use SAM files, but its much slower).

    $ samtools import mm9.fa mapped.sam mapped.bam  

