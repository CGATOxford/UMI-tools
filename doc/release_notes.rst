Release notes
=============

1.1.2
-----
07 May 2021

**Bugfix**

 - ``whitelist --filtered-out`` with SE reads threw an unassigned error. Thanks @yech1990 for rectifying this (`#453 <https://github.com/CGATOxford/UMI-tools/issues/453>`_)

Also includes a very minor update of syntax (`#455 <https://github.com/CGATOxford/UMI-tools/issues/455>`_)


1.1.1
-----
18 Nov 2020

Updates requirements for pysam version to >0.16.0.1. Thanks @sunnymouse25 (`#444 <https://github.com/CGATOxford/UMI-tools/issues/444>`_)


1.1.0
-----
4 Nov 2020

**Additional functionality**

 - Write out reads failing regex matching with extract/whitelist (see options --filtered-out, --filtered-out2). See `#328 <https://github.com/CGATOxford/UMI-tools/issues/328>`_ for motivation
 - Ignore template length with paired-end dedup/group (see option --ignore-tlen). See `#357 <https://github.com/CGATOxford/UMI-tools/issues/357>`_ for motivation. Thanks @skitcattCRUKMI
 - Ignore read pair suffixes with extract/whitelist e.g /1 or /2. (see option --ignore-read-pair-suffixes). See (`#325 <https://github.com/CGATOxford/UMI-tools/issues/325>`_, `#391 <https://github.com/CGATOxford/UMI-tools/issues/391>`_, `#418 <https://github.com/CGATOxford/UMI-tools/issues/418>`_, `PierreBSC/Viral-Track issue 9 <https://github.com/PierreBSC/Viral-Track/issues/9>`_) for motivation

**Performance**

 - Sped up error correction mapping for cell barcodes in whitelist by using BKTree. Thanks @redst4r. Note that this adds a new python dependency (pybktree) which is available via pip and conda-forge.
 - Very slight reduction in memory usage for dedup/group via bugfix to reduce the amount of reads being retained in the buffer. Thanks to @mitrinh1 for spotting this (`#428 <https://github.com/CGATOxford/UMI-tools/issues/428>`_). The bug was equivalent to hardcoding the option -buffer-whole-contig on, which ensures all reads with the same start position are grouped together for deduplication, but at the cost of not yielding reads until the end of each contig, thus increasing memory usage. As such, the bug was not detrimental to results output.

**Bugfixes**

 - Unmapped mates were not properly discarded with dedup and group. Thanks @Daniel-Liu-c0deb0t for rectifying this.

1.0.1
-----
6 Dec 2019

Debug for KeyError when some reads are missing a cell barode tag and stats output required from umi_tools dedup. See comments from @ZHUwj0 in `#281 <https://github.com/CGATOxford/UMI-tools/issues/281>`_


1.0.0
-----
14 Feb 2019

This release is intended to be a stable release with no plans for significant updates to UMI-tools functionality in the near future. As part of this release, much of the code base has been refactored. It is possible this may have introduced bugs which have not been picked up by the regression testing. If so, please raise an issue and we'll try and rectify with a minor release update ASAP.

**Documentation**

UMI-tools documentation is now available online: https://umi-tools.readthedocs.io/en/latest/index.html

Along with the previous documentation, the readthedocs pages also include new pages:

 - FAQ
 - Making use of our Alogrithmns: The API

**New knee method for whitelist**

 - The method to detect the "knee" in whitelist has been updated (`#317 <https://github.com/CGATOxford/UMI-tools/issues/317>`_). This method should always identify a threshold and is now set as the default method. Note that this knee method appears to be slightly more conservative (fewer cells above threshold) but having identified the knee, one can always re-run whitelist and use ``--set-cell-number`` to expand the whitelist if desired
 - The old method is still available via ``--knee-method=density``
 - In addition, to run the old knee method but allow whitelist to exit without error even if a suitable knee point isn't identified, use the new ``--allow-threshold-error`` option (`#249 <https://github.com/CGATOxford/UMI-tools/issues/249>`_)
 - Putative errors in CBs above the knee can be detected using ``--ed-above-threshold`` (`#309 <https://github.com/CGATOxford/UMI-tools/issues/309>`_)

**Explicit options for handling chimeric & inproper read pairs** (`#312 <https://github.com/CGATOxford/UMI-tools/issues/312>`_)

The behaviour for chimeric read pairs, inproper read pairs and unmapped reads can now be explictly set with the ``--chimeric-pairs``, ``--unpaired-reads`` and ``--unmapped-reads`` options.

**New options**

 - ``--temp-dir``: Set the directory for temporary files (`#254 <https://github.com/CGATOxford/UMI-tools/issues/254>`_)
 - ``--either-read`` & ``--either-read-resolve``: Extract the UMI from either read (`#175 <https://github.com/CGATOxford/UMI-tools/issues/175>`_)

**Misc**

 - Updates python testing version to 3.6.7 and drops python 2 testing
 - Replace deprecated imp import (`#318 <https://github.com/CGATOxford/UMI-tools/issues/318>`_)
 - Debug error with pysam <0.14 (`#319 <https://github.com/CGATOxford/UMI-tools/issues/319>`_)
 - Refactor module files
 - Moves documentation into dedicated module



0.5.5
-----

16 Nov 2018

Mainly minor debugs and improved detection of incorrect command line options. Minor updates to documentation.

 - Resolves issues correctly skipping reads which have not been
   assigned (`#191 <https://github.com/CGATOxford/UMI-tools/issues/191>`_ & `#273 <https://github.com/CGATOxford/UMI-tools/issues/273>`_).

This involves the addition of the ``--assigned-status-tag`` option

 - Testing for OSX has been dropped due to unresolved issues with travis. We hope to resurrect this in the future!

 - In line with major python packages (e.g https://www.numpy.org/neps/nep-0014-dropping-python2.7-proposal.html), support for python 2 will be dropped from January 1st 2019.


0.5.4
-----

16 Jul 2018

 - The defualt value for ``--skip_regex`` was incorrectly
   formatted. Thanks to @ekernf01 for spotting (`#231
   <https://github.com/CGATOxford/UMI-tools/issues/231>`_ / `#256
   <https://github.com/CGATOxford/UMI-tools/issues/256>`_)


0.5.3
-----

2 Jan 2018

 - Debugs wide-format output for count (`#227 <https://github.com/CGATOxford/UMI-tools/issues/227>`_). Thanks @kevin199011

0.5.2
-----

 21 Dec 2017

 - Adds options to specify a delimiter for a cell barcode or UMI which
   should be concatenated + options to specify a string splitting the
   cell barcode or UMI into multiple parts, of which only the first
   will be used. Note, this options will only work if the barcodes are
   contained in the BAM tag - if they were appended to the read name
   using umi_tools extract there is no need for these options. See
   `#217 <https://github.com/CGATOxford/UMI-tools/issues/217>`_ for
   motivation:
    - ``--umi-tag-delimiter=[STRING]``
       remove the delimeter STRING from the UMI. Defaults to None
    - ``--umi-tag-split=[STRING]``
       split UMI by STRING and take only the first portion. Defaults to None
    - ``--cell-tag-delimiter=[STRING]``
       remove the delimeter STRING from the cell barcode. Defaults to None
    - ``--cell-tag-split=[STRING]``
       split cell barcode by STRING and take only the first
       portion. Defaults to ``-`` to deal with 10X GEMs

 - Reduced memory requirements for ``count --wide-format-cell-counts``
   (`#222 <https://github.com/CGATOxford/UMI-tools/issues/222>`_)
 - Debugs issues with --bc-pattern2 (`#201
   <https://github.com/CGATOxford/UMI-tools/issues/201>`_, `#221 <https://github.com/CGATOxford/UMI-tools/issues/221>`_)
 - Updates documentation (`#204
   <https://github.com/CGATOxford/UMI-tools/issues/204>`_,
   `#210 <https://github.com/CGATOxford/UMI-tools/issues/210>`_, `#211 <https://github.com/CGATOxford/UMI-tools/issues/211>`_). Thanks @kohlkopf, @hy09 & @cbrueffer.


0.5.1
-----

16 Oct 2017

- Minor update. Improves detection of duplicate reads with paired end
  reads, reduces run time with dedup ``--output-stats`` and a few simple
  debugs.
- Improved identification of duplicate reads from paired end reads -
  will now use the position of the FIRST splice junction in the read
  (in reference coords)
  (`#187 <https://github.com/CGATOxford/UMI-tools/issues/187>`_)
- Speeds up dedup when running with ``--output-stats`` - (`#184 <https://github.com/CGATOxford/UMI-tools/issues/184>`_)
- Fixes bugs:
    - ``whitelist --set-cell-number --plot-prefix`` -> unwanted error
    - dedup gave non-informative error when input contains zero valid
      reads/read pairs. Now raises a warning but exits with status 0
      (`#190 <https://github.com/CGATOxford/UMI-tools/issues/190>`_,
      `#195 <https://github.com/CGATOxford/UMI-tools/issues/195>`_)
    - count errored if gene identifier contained a ":" (`#198 <https://github.com/CGATOxford/UMI-tools/issues/198>`_)
    - Renames ``--whole-contig option`` to ``--buffer-whole-contig`` to
      avoid confusion with `--per-contig`` option. ``--whole-contig`` option
      will still work but will not be visible in documentation (`#196 <https://github.com/CGATOxford/UMI-tools/issues/196>`_)

0.5.0
-----

18 Aug 2017

Version 0.5.0 introduces new commands to support single-cell RNA-Seq and reduces run-time. The underlying methods have not changed hence the minor release number uptick.

**UMI-tools goes single cell**

New commands for single cell RNA-Seq (scRNA-Seq):

 - ``whitelist``
    Extract cell barcodes (CB) from droplet-based scRNA-Seq fastqs and
   estimate the number of "true" CBs. Outputs a flatfile listing the
   true cell barcodes and 'error' barcodes within a set distance. See
   `#97 <https://github.com/CGATOxford/UMI-tools/issues/97>`_ for a
   motivating example. Thanks to @Hoohm for input and patience in
   testing. Thanks to @k3yavi for input in discussions about
   implementing a 'knee' method.
 - ``count``
    Count the number of reads per cell per gene after
    de-duplication. This tool uses the same underlying methods as
    group and dedup and acts to simplify scRNA-Seq read-counting with
    umi_tools. See `#114
    <https://github.com/CGATOxford/UMI-tools/issues/114>`_, `#131
    <https://github.com/CGATOxford/UMI-tools/issues/131>`_.
 - ``count_tab``
    As per count but works from a flatfile input from e.g
    featureCounts - See `#44
    <https://github.com/CGATOxford/UMI-tools/issues/44>`_, `#121
    <https://github.com/CGATOxford/UMI-tools/issues/121>`_, `#125 <https://github.com/CGATOxford/UMI-tools/issues/125>`_

In the process of creating these commands, the options for dealing
with UMIs on a "per-gene" basis have been re-jigged to make their
purpose clearer. See e.g `#127 <https://github.com/CGATOxford/UMI-tools/issues/127>`_ for a motvating example.

To perform group, dedup or count on a per-gene, basis, the ``--per-gene`` option should be provided. This must be combined with either ``--gene-tag`` if the BAM contains gene assignments in a tag, or ``--per-contig`` if the reads have been aligned to a transcriptome. In the later case, if the reads have been aligned to a transcriptome where each contig is a transcript, the option ``--gene-transcript-map`` can be used to operate at the gene level. These options are standardised across all tools such that one can easily change e.g a ``count`` command into a ``dedup`` command.

*Additional updates*

 - ``extract`` can now accept regex patterns to describe UMI +/- CB encoding in read(s). See ``--extract-method=regex`` option.

 - We have written a guide for how to use UMI-tools for scRNA-Seq analysis including estimation of the number of true CBs, flexible extraction of cell barcodes and UMIs and ``--per-cell`` read-counting as well as common workflow variations.

 - Reduced run-time
   (`#156 <https://github.com/CGATOxford/UMI-tools/issues/156>`_)

 - Introduced a hashing step to limit the scope of the edit-distance
   comparisons required to build the networks. Big thanks to @mparker2
   for this!

 - Simplified installation (`#145 <https://github.com/CGATOxford/UMI-tools/issues/145>`_)

 - Previously extensions were cythonized and compiled on the fly using
   ``pyximport``, requiring users to have access to the install
   directory the first time the extension was required. Now the
   cythonized extension is provided, and is compiled at install-time.


0.4.4
-----

8 May 2017

 - Tweaks the way group handles paired end BAMs. To simplify the
   process and ensure all reads are written out, the paired end read
   (read 2) is now outputted without a group or UMI tag. (`#115
   <https://github.com/CGATOxford/UMI-tools/issues/115>`_).
 - Introduces the ``--skip-tags-regex`` option to enable users to skip
   descriptive gene tags, such as "Unassigned" when using
   the --gene-tag option. See `#108
   <https://github.com/CGATOxford/UMI-tools/issues/108>`_.

*Bugfixes:*
 - If the ``--transcript-gene-map`` included transcripts not observed in the BAM, this caused an error when trying to retrieve reads aligned to the transcript. This has been resolved. See `#109 <https://github.com/CGATOxford/UMI-tools/issues/109>`_
 - Allow output to zipped file with extract using python 3 `#104 <https://github.com/CGATOxford/UMI-tools/issues/104>`_
 - Improved test coverage (``--chrom`` and ``--gene-tag``
   options). Thanks @MarinusVL for kindly sharing a BAM with gene
   tags.

0.4.3
-----

28 Mar 2017

 - Improves run time for large networks (see `#94
   <https://github.com/CGATOxford/UMI-tools/issues/94>`_, `#31
   <https://github.com/CGATOxford/UMI-tools/issues/31>`_). Thanks to
   @gpratt for identifying the issue and implementing the solution



0.4.2
-----

22 Mar 2017

 - When using the directional method with the group command, the 'top' UMI within each group was not always the most abundant (see comments in `#96 <https://github.com/CGATOxford/UMI-tools/issues/96>`_). This has now been resolved

0.4.1
-----

16 Mar 2017 

 - Due to a bug in ``pysam.fetch()`` paired end files with a large number
   of contigs could take a long time to process (see `#93
   <https://github.com/CGATOxford/UMI-tools/issues/93>`_). This has
   now been resolved. Thanks to @gpratt for spotting and resolving
   this.


0.4.0
-----

9 Mar 2017

*Added functionality:*

 - Deduplicating on gene ids (`#44
   <https://github.com/CGATOxford/UMI-tools/issues/44>`_` for
   motivation)
   - The user can now group/dedup according to the gene which the read
     aligns to. This is useful for single cell RNA-Seq methods such as
     e.g CEL-Seq where the position of the read on a transcript may be
     different for reads generated from the same initial molecule. The
     following options may be used define the gene_id for each read:
      - ``--per-gene``
      - ``--gene-transcript-map``
      - ``--gene-tag``

 - Working with BAM tags (`#73
   <https://github.com/CGATOxford/UMI-tools/issues/73>`_,
   `#76 <https://github.com/CGATOxford/UMI-tools/issues/76>`_,
   `#89 <https://github.com/CGATOxford/UMI-tools/issues/89>`_):

 - UMIs can now be extracted from the BAM tags and `group` will add a
   tag to each read describing the read group and UMI. See following
   options for controlling this behaviour:
    - ``--extract-umi-method``
    - ``--umi-tag``
    - ``--umi-group-tag``

 - Ouput unmapped reads
   (`#78 <https://github.com/CGATOxford/UMI-tools/issues/78>`_)
    The group command will now output unmapped reads if
    the ``--output-unmapped`` is supplied. These reads will not be
    assigned to any group.

 - bug fixes for ``group`` command
   (`#67 <https://github.com/CGATOxford/UMI-tools/issues/67>`_, `#81
   <https://github.com/CGATOxford/UMI-tools/issues/81>`_)
 - updated documentation
   (`#77 <https://github.com/CGATOxford/UMI-tools/issues/77>`_,
   `#79 <https://github.com/CGATOxford/UMI-tools/issues/79>`_ )

0.3.6
-----

1 Feb 2017

*Improves the group command:*
 - Adds the ``--subset option`` as per the dedup command (`#74
   <https://github.com/CGATOxford/UMI-tools/issues/74>`_)
 - Corrects the flatfile output from the dedup command (`#72
   <https://github.com/CGATOxford/UMI-tools/issues/72>`_)



0.3.5
-----

27 Jan 2017

 - The code has been tweaked to improve run-time. See `#69
   <https://github.com/CGATOxford/UMI-tools/issues/69>`_ for a
   discussion about the changes implemented.


0.3.4
-----

23 Jan 2017

 - Corrects the edit distance comparison used to generate the network
   for the ``directional`` method.
  - This will only affect results generated using the directional
    method and ``--edit-distance-threshold`` >1.
  - Previously, using the ``directional`` method with the option
    ``--edit-distance-threshold`` set to > 1 did not return the
    expected set of de-duplicated reads. If you have used the
    ``directional`` method with a threshold >1, we recommend updating
    UMI-tools and re-running dedup.


0.3.3
-----

 19 Jan 2017

 - Debugs ``python 3`` compatibility issues
 - Adds ``python 3`` tests


0.3.2
-----

17 Jan 2017)

*Minor bump:*
 - Resolves setuptools-based installation issue


0.3.1
-----

1 Dec 2016

*Version bump to allow pypi update. No code changes*


0.3.0
-----

1 Dec 2016

 - Adds the new ``group`` command to group PCR duplicates and return
   the groups in a tagged BAM file and/or flat file format. This was
   motivated by multiple requests to group PCR duplicated reads for
   downstream processes, e,g `#45
   <https://github.com/CGATOxford/UMI-tools/issues/45>`_, `#54
   <https://github.com/CGATOxford/UMI-tools/issues/54>`_. Special
   thanks to Nils Koelling (@koelling) for testing the group command.


 - Adds the --umi-separator option for dedup and group for workflow
   where umi_tools extract is not used to extract the UMI. This was
   motivated by `#58 <https://github.com/CGATOxford/UMI-tools/issues/58>`_


0.2.6
-----

8 Nov 2016

 - directional-adjacency method is renamed directional

0.2.5
-----

2 Nov 2016

 - Debugs writing out paired end
 - Debugs installation

0.2.3
-----

7 Jun 2016

 - Debugs pip installation


0.2.0
-----

31 May 2016

*extract*
 - New feature: Filter out read by UMI base-call quality score
   ``--quality-threshold`` and ``--quality-encoding`` options (`#29
   <https://github.com/CGATOxford/UMI-tools/issues/29>`_, `#33  <https://github.com/CGATOxford/UMI-tools/issues/33>`_)

*dedup*
 - Improved performance for paired end files (`#31
   <https://github.com/CGATOxford/UMI-tools/issues/31>`_, `#35  <https://github.com/CGATOxford/UMI-tools/issues/35>`_)

0.0.11
------

23 May 2016

 - Debugs read extraction from 3' end

0.0.10
------

 - Improved memory performace for UMI extraction from paired end reads

0.0.9
-----
29 Apr 2016

**UMI-Tools Manuscript Release**

 - Merge pull request `#18 <https://github.com/CGATOxford/UMI-tools/issues/18>`_ from CGATOxford/TS-RefactorTools
