# FAQ
 - **Why is my `umi_tool group`/`dedup` command taking so long?**
The time taken to resolve each position is dependent on how many unique UMIs there are and how closely related they are since these factors will affect how large the network of connected UMIs is. Some of the factors which can affect run time are:
    1. UMI length (shorter => fewer possible UMIs => increased connectivity between UMIs => larger networks => longer run time)
    2. Sequencing error rate (higher => more error UMIs => larger networks => longer run time)
    3. Sequencing depth (higher = greater proportion of possible UMIs observed => larger network => longer run time)
    4. Running `umi_tools dedup --output-stats` requires a considerable amount of time and memory to generate the null distributions. If you want these stats, consider obtaining stats for just a single contig, eg. `--chrom=chr22`.
    5. If you are doing single-cell RNA-seq and you have reads from more than one cell in your BAM file, make sure you are running with the `--per-cell` swtich.
&nbsp;
- **Why is my `umi_tool group`/`dedup` command taking so much memory?**
There are a few reasons why your command could require an excessive amount of memory:
    1. Most of the above factors increasing run time can also increase memory
    2. Running `umi_tools dedup` with `chimeric-reads=use` (default) can take a large amount of memory. This is because `dedup` will wait until a whole contig has been deduplicated before re-parsing the reads from the contig and writing out the read2s which match the retained read1s. For chimeric read pairs, the read2s will not be found on the same contig and will be kept in a buffer of "orphan" read2s which may take up a lot of memory. Consider using `chimeric-reads=discard` instead to discard chimeric read pairs. Similarly, use of `--unmapped-reads=use` with `--paired` can also increase memory requirements.
&nbsp;
 - **Can I run `umi_tools` with parallel threads?**
Not yet! This is something we have been discussing for a while (see [#203](https://github.com/CGATOxford/UMI-tools/issues/203) & [#257](https://github.com/CGATOxford/UMI-tools/issues/257)) but haven't found the time to actually implement. If you'd like to help us out, get in touch!
&nbsp;
&nbsp;
- **What's the correct regex to use for technique X?**
Here is a table of the techniques we have come across that are not easily processed with the basic barcode pattern syntax, and the correct regex's to use with them:

   | Technique | regex |
   | --------- | ------ |
   | inDrop    | `(?P<cell_1>.{8,12})(?P<discard_1>GAGTGATTGCTTGTGACGCCTT)(?P<cell_2>.{8})(?P<umi_1>.{6})T{3}.*` |
   | ddSeq/Sureccell | `(?P<discard_1>.{0,5})(?P<cell_1>.{6})(?P<discard_2>TAGCCATCGCATTGC){e<=1}(?P<cell_2>.{6})(?P<discard_3>TACCTCTGAGCTGAA){e<=1}(?P<cell_3>.{6})(?P<discard_4>ACG)(?P<umi_1>.{8})(?P<discard_5>GAC).*` |
   | SPLiT-seq | `(?P<umi_1>.{10})(?<cell_1>.{8})(?P<discard_1>GTGGCCGATGTTTCGCATCGGCGTACGACT)(?P<cell_2>.{8})(?<discard_2>ATCCACGTGCTTGAGAGGCCAGAGCATTCG)(?P<cell_3>.{8})` |
   | BD-Rhapsody | `(?<cell_1>.{9})(?<discard_1>.{12})(?<cell_2>.{9})(?<discard_2>.{13})(?<cell_3>.{9})(?<umi_1>.{8})T+` |

   If you know of other, please drop us a PR with them!

- **Can I use `umi_tools` to determine consensus sequences?**
Right now, you can use `umi_tools group` to identify the duplicated read groups. From this, you can then derive consensus sequences as you wish. We have discussed adding consense sequence calling as a separate `umi_tools` command (see [#203](https://github.com/CGATOxford/UMI-tools/issues/181)). If you'd like to help us out, get in touch!
&nbsp;
&nbsp;
- **What do the `--per-gene`, `--gene-tag` and `--per-contig` options do in `umi_tools group`/`dedup`/`count`?**
These options are designed to handle samples from sequence library protocols where fragmentation occurs post amplification, e.g CEL-Seq for single cell RNA-Seq. For such sample, mapping coordinates of duplicate reads will no longer be identical and `umi_tools` can instead use the gene to which the read is mapped.
This behaviour is switched on with `--per-gene`. `umi_tools` then needs to be told if the reads are directly mapped to a transcriptome (`--per-contig`), or mapped to a genome and the transcript/gene assignment is contained in a read tag (`--gene-tag=[TAG]`). If you have assigned reads to transcripts but want to use the gene IDs to determine UMI groups, there is a further option to provide a file mapping transcript IDs to gene IDs (`--gene-transcript-map`). Finally, for single cell RNA-Seq, you must specify `--per-cell`. See `umi_tools dedup`/`group`/`count --help` for details of further related options and the [UMI-tools single cell RNA-Seq guide](https://github.com/CGATOxford/UMI-tools/blob/%7BTS%7D-AddFAQ/doc/Single_cell_tutorial.md).
&nbsp;
&nbsp;
- **Should I use `--per-gene` with my sequencing from method X**
It can be difficult to work this out sometimes! So far we have come across the following technqies that require the use of `--per-gene`: CEL-seq2, SCRB-seq, 10x Chromium, inDrop, Drop-seq and SPLiT-seq. Let us know if you know of more.
