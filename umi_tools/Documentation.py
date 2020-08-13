# generic docstring for whitelist/extract
GENERIC_DOCSTRING_WE = '''
Barcode extraction
------------------

""""""""""""""""
``--bc-pattern``
""""""""""""""""
      Pattern for barcode(s) on read 1. See ``--extract-method``

"""""""""""""""""
``--bc-pattern2``
"""""""""""""""""
      Pattern for barcode(s) on read 2. See ``--extract-method``

""""""""""""""""""""
``--extract-method``
""""""""""""""""""""
      There are two methods enabled to extract the umi barcode (+/-
      cell barcode). For both methods, the patterns should be provided
      using the ``--bc-pattern`` and ``--bc-pattern2`` options.x

 - ``string``
       This should be used where the barcodes are always in the same
       place in the read.

       - N = UMI position (required)
       - C = cell barcode position (optional)
       - X = sample position (optional)

       Bases with Ns and Cs will be extracted and added to the read
       name. The corresponding sequence qualities will be removed from
       the read. Bases with an X will be reattached to the read.

       E.g. If the pattern is `NNNNCC`,
       Then the read::

           @HISEQ:87:00000000 read1
           AAGGTTGCTGATTGGATGGGCTAG
           +
           DA1AEBFGGCG01DFH00B1FF0B

       will become::

           @HISEQ:87:00000000_TT_AAGG read1
           GCTGATTGGATGGGCTAG
           +
           1AFGGCG01DFH00B1FF0B

       where 'TT' is the cell barcode and 'AAGG' is the UMI.

 - ``regex``
       This method allows for more flexible barcode extraction and
       should be used where the cell barcodes are variable in
       length. Alternatively, the regex option can also be used to
       filter out reads which do not contain an expected adapter
       sequence. UMI-tools uses the regex module rather than the more
       standard re module since the former also enables fuzzy matching

       The regex must contain groups to define how the barcodes are
       encoded in the read. The expected groups in the regex are:

       umi_n = UMI positions, where n can be any value (required)
       cell_n = cell barcode positions, where n can be any value (optional)
       discard_n = positions to discard, where n can be any value (optional)

       UMI positions and cell barcode positions will be extracted and
       added to the read name. The corresponding sequence qualities
       will be removed from the read.

       Discard bases and the corresponding quality scores will be
       removed from the read. All bases matched by other groups or
       components of the regex will be reattached to the read sequence

       For example, the following regex can be used to extract reads
       from the Klein et al inDrop data::

           (?P<cell_1>.{8,12})(?P<discard_1>GAGTGATTGCTTGTGACGCCTT)(?P<cell_2>.{8})(?P<umi_1>.{6})T{3}.*

       Where only reads with a 3' T-tail and `GAGTGATTGCTTGTGACGCCTT` in
       the correct position to yield two cell barcodes of 8-12 and 8bp
       respectively, and a 6bp UMI will be retained.

       You can also specify fuzzy matching to allow errors. For example if
       the discard group above was specified as below this would enable
       matches with up to 2 errors in the discard_1 group.

       ::

           (?P<discard_1>GAGTGATTGCTTGTGACGCCTT){s<=2}

       Note that all UMIs must be the same length for downstream
       processing with dedup, group or count commands


""""""""""""
``--3prime``
""""""""""""
       By default the barcode is assumed to be on the 5' end of the
       read, but use this option to sepecify that it is on the 3' end
       instead. This option only works with ``--extract-method=string``
       since 3' encoding can be specified explicitly with a regex, e.g
       ``.*(?P<umi_1>.{5})$``

""""""""""""""
``--read2-in``
""""""""""""""
        Filename for read pairs

""""""""""""""""""
``--filtered-out``
""""""""""""""""""
        Write out reads not matching regex pattern or cell barcode
        whitelist to this file

"""""""""""""""""""
``--filtered-out2``
"""""""""""""""""""
        Write out read pairs not matching regex pattern or cell barcode
        whitelist to this file

""""""""""""""
``--ignore-read-pair-suffixes``
""""""""""""""
       Ignore \1 and \2 read name suffixes. Note that this options is
       required if the suffixes are not whitespace separated from the
       rest of the read name
'''


# generic docstring for group/dedup/count/count_tab
GENERIC_DOCSTRING_GDC = '''

Extracting barcodes
-------------------

It is assumed that the FASTQ files were processed with `umi_tools
extract` before mapping and thus the UMI is the last word of the read
name. e.g::

    @HISEQ:87:00000000_AATT

where `AATT` is the UMI sequeuence.

If you have used an alternative method which does not separate the
read id and UMI with a "_", such as bcl2fastq which uses ":", you can
specify the separator with the option ``--umi-separator=<sep>``,
replacing <sep> with e.g ":".

Alternatively, if your UMIs are encoded in a tag, you can specify this
by setting the option --extract-umi-method=tag and set the tag name
with the --umi-tag option. For example, if your UMIs are encoded in
the 'UM' tag, provide the following options:
``--extract-umi-method=tag`` ``--umi-tag=UM``

Finally, if you have used umis to extract the UMI +/- cell barcode,
you can specify ``--extract-umi-method=umis``

The start position of a read is considered to be the start of its alignment
minus any soft clipped bases. A read aligned at position 500 with
cigar 2S98M will be assumed to start at position 498.

""""""""""""""""""""""""
``--extract-umi-method``
""""""""""""""""""""""""
      How are the barcodes encoded in the read?

      Options are:

      - read_id (default)
            Barcodes are contained at the end of the read separated as
            specified with ``--umi-separator`` option

      - tag
            Barcodes contained in a tag(s), see ``--umi-tag``/``--cell-tag``
            options

      - umis
            Barcodes were extracted using umis (https://github.com/vals/umis)

"""""""""""""""""""""""""""""""
``--umi-separator=[SEPARATOR]``
"""""""""""""""""""""""""""""""
      Separator between read id and UMI. See ``--extract-umi-method``
      above. Default=``_``

"""""""""""""""""""
``--umi-tag=[TAG]``
"""""""""""""""""""
      Tag which contains UMI. See ``--extract-umi-method`` above

"""""""""""""""""""""""""""
``--umi-tag-split=[SPLIT]``
"""""""""""""""""""""""""""
      Separate the UMI in tag by SPLIT and take the first element

"""""""""""""""""""""""""""""""""""
``--umi-tag-delimiter=[DELIMITER]``
"""""""""""""""""""""""""""""""""""
      Separate the UMI in by DELIMITER and concatenate the elements

""""""""""""""""""""
``--cell-tag=[TAG]``
""""""""""""""""""""
      Tag which contains cell barcode. See `--extract-umi-method` above

""""""""""""""""""""""""""""
``--cell-tag-split=[SPLIT]``
""""""""""""""""""""""""""""
      Separate the cell barcode in tag by SPLIT and take the first element

""""""""""""""""""""""""""""""""""""
``--cell-tag-delimiter=[DELIMITER]``
""""""""""""""""""""""""""""""""""""
      Separate the cell barcode in by DELIMITER and concatenate the elements


UMI grouping options
---------------------------

""""""""""""
``--method``
""""""""""""
    What method to use to identify group of reads with the same (or
    similar) UMI(s)?

    All methods start by identifying the reads with the same mapping position.

    The simplest methods, unique and percentile, group reads with
    the exact same UMI. The network-based methods, cluster, adjacency and
    directional, build networks where nodes are UMIs and edges connect UMIs
    with an edit distance <= threshold (usually 1). The groups of reads
    are then defined from the network in a method-specific manner. For all
    the network-based methods, each read group is equivalent to one read
    count for the gene.

      - unique
          Reads group share the exact same UMI

      - percentile
          Reads group share the exact same UMI. UMIs with counts < 1% of the
          median counts for UMIs at the same position are ignored.

      - cluster
          Identify clusters of connected UMIs (based on hamming distance
          threshold). Each network is a read group

      - adjacency
          Cluster UMIs as above. For each cluster, select the node (UMI)
          with the highest counts. Visit all nodes one edge away. If all
          nodes have been visited, stop. Otherwise, repeat with remaining
          nodes until all nodes have been visted. Each step
          defines a read group.

      - directional (default)
          Identify clusters of connected UMIs (based on hamming distance
          threshold) and umi A counts >= (2* umi B counts) - 1. Each
          network is a read group.

"""""""""""""""""""""""""""""
``--edit-distance-threshold``
"""""""""""""""""""""""""""""
       For the adjacency and cluster methods the threshold for the
       edit distance to connect two UMIs in the network can be
       increased. The default value of 1 works best unless the UMI is
       very long (>14bp).

"""""""""""""""""""""""
``--spliced-is-unique``
"""""""""""""""""""""""
       Causes two reads that start in the same position on the same
       strand and having the same UMI to be considered unique if one is spliced
       and the other is not. (Uses the 'N' cigar operation to test for
       splicing).

"""""""""""""""""""""""""
``--soft-clip-threshold``
"""""""""""""""""""""""""
       Mappers that soft clip will sometimes do so rather than mapping a
       spliced read if there is only a small overhang over the exon
       junction. By setting this option, you can treat reads with at least
       this many bases soft-clipped at the 3' end as spliced. Default=4.

""""""""""""""""""""""""""""""""""""""""""""""
``--multimapping-detection-method=[NH/X0/XT]``
""""""""""""""""""""""""""""""""""""""""""""""
      If the sam/bam contains tags to identify multimapping reads, you can
      specify for use when selecting the best read at a given loci.
      Supported tags are "NH", "X0" and "XT". If not specified, the read
      with the highest mapping quality will be selected.

"""""""""""""""""
``--read-length``
"""""""""""""""""
      Use the read length as a criteria when deduping, for e.g sRNA-Seq.


Single-cell RNA-Seq options
---------------------------

""""""""""""""
``--per-gene``
""""""""""""""
      Reads will be grouped together if they have the same gene.  This
      is useful if your library prep generates PCR duplicates with non
      identical alignment positions such as CEL-Seq. Note this option
      is hardcoded to be on with the count command. I.e counting is
      always performed per-gene. Must be combined with either
      ``--gene-tag`` or ``--per-contig`` option.

""""""""""""""
``--gene-tag``
""""""""""""""
      Deduplicate per gene. The gene information is encoded in the bam
      read tag specified

"""""""""""""""""""""""""
``--assigned-status-tag``
"""""""""""""""""""""""""
      BAM tag which describes whether a read is assigned to a
      gene. Defaults to the same value as given for ``--gene-tag``

"""""""""""""""""""""
``--skip-tags-regex``
"""""""""""""""""""""
      Use in conjunction with the ``--assigned-status-tag`` option to
      skip any reads where the tag matches this regex.  Default
      (``"^[__|Unassigned]"``) matches anything which starts with "__"
      or "Unassigned":

""""""""""""""""
``--per-contig``
""""""""""""""""
      Deduplicate per contig (field 3 in BAM; RNAME).
      All reads with the same contig will be considered to have the
      same alignment position. This is useful if you have aligned to a
      reference transcriptome with one transcript per gene. If you
      have aligned to a transcriptome with more than one transcript
      per gene, you can supply a map between transcripts and gene
      using the ``--gene-transcript-map`` option

"""""""""""""""""""""""""
``--gene-transcript-map``
"""""""""""""""""""""""""
      File mapping genes to transcripts (tab separated), e.g::

          gene1   transcript1
          gene1   transcript2
          gene2   transcript3

""""""""""""""
``--per-cell``
""""""""""""""
      Reads will only be grouped together if they have the same cell
      barcode. Can be combined with ``--per-gene``.

SAM/BAM Options
---------------

"""""""""""""""""""""
``--mapping-quality``
"""""""""""""""""""""
      Minimium mapping quality (MAPQ) for a read to be retained. Default is 0.

""""""""""""""""""""
``--unmapped-reads``
""""""""""""""""""""
     How should unmapped reads be handled. Options are:
      - discard (default)
          Discard all unmapped reads
      - use
          If read2 is unmapped, deduplicate using read1 only. Requires
          ``--paired``
      - output
          Output unmapped reads/read pairs without UMI
          grouping/deduplication. Only available in umi_tools group

""""""""""""""""""""
``--chimeric-pairs``
""""""""""""""""""""
     How should chimeric read pairs be handled. Options are:
      - discard
          Discard all chimeric read pairs
      - use (default)
          Deduplicate using read1 only
      - output
          Output chimeric read pairs without UMI
          grouping/deduplication.  Only available in umi_tools group

""""""""""""""""""""
``--unpaired-reads``
""""""""""""""""""""
     How should unpaired reads be handled. Options are:
      - discard
          Discard all unpaired reads
      - use (default)
          Deduplicate using read1 only
      - output
          Output unpaired reads without UMI
          grouping/deduplication. Only available in umi_tools group

""""""""""""""""
``--ignore-umi``
""""""""""""""""
      Ignore the UMI and group reads using mapping coordinates only

""""""""""""
``--subset``
""""""""""""
      Only consider a fraction of the reads, chosen at random. This is useful
      for doing saturation analyses.

"""""""""""
``--chrom``
"""""""""""
      Only consider a single chromosome. This is useful for
      debugging/testing purposes


Input/Output Options
---------------------

"""""""""""""""""""""""
``--in-sam, --out-sam``
"""""""""""""""""""""""
      By default, inputs are assumed to be in BAM format and outputs are written
      in BAM format. Use these options to specify the use of SAM format for
      input or output.

""""""""""""
``--paired``
""""""""""""
       BAM is paired end - output both read pairs. This will also
       force the use of the template length to determine reads with
       the same mapping coordinates.

'''


# generic docstring for dedup/group
GROUP_DEDUP_GENERIC_OPTIONS = '''
Group/Dedup options
-------------------

""""""""""""""""""""
``--no-sort-output``
""""""""""""""""""""
       By default, output is sorted. This involves the
       use of a temporary unsorted file since reads are considered in
       the order of their start position which may not be the same
       as their alignment coordinate due to soft-clipping and reverse
       alignments. The temp file will be saved (in ``--temp-dir``) and deleted
       when it has been sorted to the outfile. Use this option to turn
       off sorting.


"""""""""""""""""""""""""
``--buffer-whole-contig``
"""""""""""""""""""""""""
      forces dedup to parse an entire contig before yielding any reads
      for deduplication. This is the only way to absolutely guarantee
      that all reads with the same start position are grouped together
      for deduplication since dedup uses the start position of the
      read, not the alignment coordinate on which the reads are
      sorted. However, by default, dedup reads for another 1000bp
      before outputting read groups which will avoid any reads being
      missed with short read sequencing (<1000bp).

'''
