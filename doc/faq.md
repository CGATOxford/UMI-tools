# FAQ
 - **Why is my `umi_tool group`/`dedup` command taking so long?**
The time taken to resolve each position is dependent on how many unique UMIs there are and how closely related they are since these factors will affect how large the network of connected UMIs is. Some of the factors which can affect run time are:
    1. UMI length (shorter => fewer possible UMIs => increased connectivity between UMIs => larger networks => longer run time)
    2. Sequencing error rate (higher => more error UMIs => larger networks => longer run time)
    3. Sequencing depth (higher = greater proportion of possible UMIs observed => larger network => longer run time)
&nbsp;
- **Why is my `umi_tool group`/`dedup` command taking so much memory?**
There are a few reasons why your command could require an excessive amount of memory:
    1. Most of the above factors increasing run time can also increase memory
    2. Running `umi_tools dedup --output-stats` requires a considerable amount of memory to track the UMI statistics. If you want these stats, consider obtaining stats for just a single contig, eg. `--chrom=chr22`
    3. Running `umi_tools dedup` with `chimeric-reads=use` (default) can take a large amount of memory. This is because `dedup` will wait until a whole contig has been deduplicated before re-parsing the reads from the contig and writing out the read2s which match the retain read1s. For chimeric read pairs, the read2s will not be found on the same contig and will be kept in a buffer of unfound read2s which may take up a lot of memory. Consider using `chimeric-reads=discard` instead to discard chimeric read pairs. Similarly, use of `--unmapped-reads=use` with `--paired` can also increase memory requirements.
&nbsp;
 - **Can I run `umi_tools` with parallel threads?**
Not yet! This is something we have been discussing for a while (see [#203](https://github.com/CGATOxford/UMI-tools/issues/203) & [#257](https://github.com/CGATOxford/UMI-tools/issues/257)) but haven't found the time to actually implement. If you'd like to help us out, get in touch!
&nbsp;
- **What's this regex thing all about in `umi_tools extract`/`whitelist`?**
Regexes provide a more flexible way to describe the pattern of UMI +/- cell barcode in the reads and can be used with (`--extract-method=regex`). The regex must contain groups to define how the barcodes are encoded in the read. The allowable groups in the regex are:
    - umi_n = UMI positions, where n can be any value (required)
    - cell_n = cell barcode positions, where n can be any value (optional)
    - discard_n = positions to discard, where n can be any value (optional)

  For example, the following regex can be used to extract reads from the Klein et al inDrop data:
`(?P<cell_1>.{8,12})(?P<discard_1>GAGTGATTGCTTGTGACGCCTT)(?P<cell_2>.{8})(?P<umi_1>.{6})T{3}.*`
Regexes provide a number of advantages over the simpler "string" extraction method:
    1. Reads not matching the regex will be discarded. In the above, this is used to filter reads which do not contain the adapter sequence between `cell_1` and `cell_2` groups.
    2. Variable cell barcode lengths can be encoded (see `cell_1` group above)
    3. Finally, regexes allow fuzzy matching (error-aware). For example to allow up to 2 errors to the inDrop adapter sequence:
    `(?P<discard_1>GAGTGATTGCTTGTGACGCCTT){s<=2}`
  Note that to enable fuzzy matching, `umi_tools` uses the [`regex`](https://pypi.org/project/regex/) library rather than the more standard `re` library.
&nbsp;
- **Can I use `umi_tools` to determine consensus sequences**
Right now, you can use `umi_tools group` to identify the duplicated read groups. From this, you can then derive consensus sequences as you wish. We have discussed adding consense sequence calling as a separate `umi_tools` command (see [#203](https://github.com/CGATOxford/UMI-tools/issues/181)). If you'd like to help us out, get in touch!
&nbsp;
- **What do the `--per-gene`, `--gene-tag` and `--per-contig` options do?**
  These options are designed to handle samples from sequence library protocols where amplification occurs post fragmentation, e.g CEL-Seq for single cell RNA-Seq. For such sample, mapping coordinates of duplicate reads will no longer be identical and `umi_tools` can instead use the gene to which the read is mapped.
This behaviour is switched on with `--per-gene`. `umi_tools` then needs to be told if the reads are directly mapped to a transcriptome (`--per-contig`), or mapped to a genome and the transcript/gene assignment is contained in a read tag (`--gene-tag=[TAG]`). There is a further option to provide a file mapping transcript IDs to gene IDs (`--gene-transcript-map`) if you have assigned reads to transcripts but want to use the gene IDs to determine UMI groups. Finally, for single cell RNA-Seq, you must specify `--per-cell`. See `umi_tools dedup`/`group`/`count --help` for details of further related options and the [UMI-tools single cell RNA-Seq guide](https://github.com/CGATOxford/UMI-tools/blob/%7BTS%7D-AddFAQ/doc/Single_cell_tutorial.md).
