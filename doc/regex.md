# Specifiying cell barcode and UMI

There are two ways to specify the location of cell barcodes and UMI bases when using `whitelist` and `extract`: string/basic (default) or regex. Sting or basic mode is simple and easy to understand and can be used for most techniques, like iCLIP, 10x chromium or Drop-seq. The regex mode is a little harder to explain, but allows a great deal of flexibility in specifying the layout of your read.

## Basic string mode

In the default basic mode we specify the location of barcodes and UMIs in a read using a simple string. `N`s specify the location of bases to be treated as UMIs, `C`s as bases to treated as cell barcode and `X`s as bases that are neither and that should be retained on the read. By default this pattern is applied to the 5' end of the read, but we can tell `extract` to look on the 3' end of the read using `--3-prime`. 

So for example in 10x Chromium, read one consists of 16 bases of cell barcode followed by 10 bases of UMI, so the correct pattern to use is `CCCCCCCCCCCCCCCCNNNNNNNNNN`:

```
@ST-K00126:308:HFLYFBBXX:1:1101:25834:1173 1:N:0:NACCACCA
NGGGTCAGTCTAGTGTGGCGATTCAC
+
#AAFFJJJJJJJJJJJJJJJJJJJJJ

-------CB-------|---UMI---
```

To understand how `X`s are treated, consider the read layout for an iCLIP experiment following the protocol of [Mueller-Mcnicoll et al,  _Genes Dev_  (2016) 30: 553](http://genesdev.cshlp.org/content/30/5/553). Here at the start of each read 1 there are three UMI bases. This is then followed by 4 bases representing a library barcode, and then there are two more UMI bases. The correct pattern here is `NNNXXXXNN`.

```
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
```

## Regex (regular expression) mode

Regexes provide a more flexible way to describe the pattern of UMI +/- cell barcode in the reads and can be used with (`--extract-method=regex`). Its needed for techniques such as inDrop, SLiT-seq and ddSeq among others.  If you know nothing about regular expressions, see [regular what-nows?](#regular-what-nows) below. 

Regexes provide a number of advantages over the simpler "string" extraction method: 

   1. Reads not matching the regex will be discarded. E.g. this can be used to filter reads which do not contain adapter sequences.
   2. Variable cell barcode lengths can be encoded.
   3. Finally, regexes allow fuzzy matching (error-aware) in, for example, adaptors.
  Note that to enable fuzzy matching, `umi_tools` uses the [`regex`](https://pypi.org/project/regex/) library rather than the more standard `re` library.
  
The regex must contain named capture groups to define how the barcodes are encoded in the read. Named capture groups are a non-standard regex feature available in the python regex dialect. A capture group is a sub part of a pattern enclosed in brackets. Matches to sub-pattern are "captured" and can be extracted for reuse. So the pattern `(.{4})TTTTT` will match any four characters followed by 5 Ts, and return what those 4 characters were. In most cases we refer to a capture group by its positions (1st group, 2nd group etc). *Named* capture groups allow use to give names to each group using a `(?P<name>` syntax. Thus, `(?P<prefix>.{4})TTTTT` matches the same as the pattern above, but the captured group is given the name "prefix". 

When passing a regex to `whitelist`/`extract`, the allowable groups in the regex are:

   * `umi_n` = UMI positions, where n can be any value (required)
   * `cell_n` = cell barcode positions, where n can be any value (optional)
   * `discard_n` = positions to discard, where n can be any value (optional)
 
 We specify fuzzy matching by adding something like `{s<=X}` after a group. This specifies that the group should be matched with up to X **s**ubstitutions. The allowed error types are `s`: substitutions, `i`: insertions, `d`:deletions, `e`: any error (or the Levenshtein distance). See the  [`regex`](https://pypi.org/project/regex/) package documentation for more details.

 ### Example: extracting UMIs

Suppose we have a FASTQ file with reads with a 4 base UMI at their 5' end and a 4 base UMI at their 3' end. For example:

```
@EWSim-1.1-umi5-reada-umix
AAAAATGGCATCCACCGATTTCTCCAAGATTGAACGTA
+
IHHHIIIHHHHIHHIIHIIHIIHHIHHHIIIIIIHIII
@EWSim-2.1-umi5-reada-umiy
AAAAATGGCATCCACCGATTTCTCCAAGATTGAAATAT
+
GABIFCD@ABEAA?EAHH?AFA@IEDGGA@CCDGFI@C
@EWSim-3.1-umi5-readb-umix
AAAATCTAGATTAGAAAGATTGACCTCATTAACGTA
+
E?IHIIH@DAC?CIAGDCIIEF@BIEDG?EDH@I?A
```

To extract these UMIs, we can use the regular expression:

```
^(?P<umi_1>.{4}).+(?P<umi_2>.{4})$
```

To break this down:

* `^(?P<umi_1>.{4})` matches the 4 bases at the start of the read and extracts this as a UMI, named `<umi_1>`.
* `.+` matches all the bases between the first 4 bases and last 4 bases (`+` indicates that there must be at least one).
* `(?P<umi_2>.{4})$` matches the 4 bases at the end of the read and extracts this as a UMI, named `<umi_2>`.

Applying this to each example read in our FASTQ file would give the following values:

```
@EWSim-1.1-umi5-reada-umix
AAAAATGGCATCCACCGATTTCTCCAAGATTGAACGTA # Read
AAAA                                   # <umi_1>
    ATGGCATCCACCGATTTCTCCAAGATTGAA     # Read-post UMI extraction
                                  CGTA # <umi_2>
@EWSim-2.1-umi5-reada-umiy
AAAAATGGCATCCACCGATTTCTCCAAGATTGAAATAT # Read
AAAA                                   # <umi_1>
    ATGGCATCCACCGATTTCTCCAAGATTGAA     # Read-post UMI extraction
                                  ATAT # <umi_2>
@EWSim-3.1-umi5-readb-umix
AAAATCTAGATTAGAAAGATTGACCTCATTAACGTA   # Read
AAAA                                   # <umi_1>
    TCTAGATTAGAAAGATTGACCTCATTAA       # Read-post UMI extraction
                                CGTA   # <umi_2>
```

If we were to use `umi_tools extract` to extract these UMIs, using this regex, our resulting FASTQ file would look like the following:

```
@EWSim-1.1-umi5-reada-umix_AAAACGTA
ATGGCATCCACCGATTTCTCCAAGATTGAA
+
IIIHHHHIHHIIHIIHIIHHIHHHIIIIII
@EWSim-2.1-umi5-reada-umiy_AAAAATAT
ATGGCATCCACCGATTTCTCCAAGATTGAA
+
FCD@ABEAA?EAHH?AFA@IEDGGA@CCDG
@EWSim-3.1-umi5-readb-umix_AAAACGTA
TCTAGATTAGAAAGATTGACCTCATTAA
+
IIH@DAC?CIAGDCIIEF@BIEDG?EDH
```

The extracted UMIs are appended together (`<umi_1><umi2>`) and are appended to the end of the head header with a delimiter, `_`.

 ### Example: extracting barcodes and UMIs

Suppose we have a multiplexed FASTQ file with reads with a 4 base UMI at their 5' end and a 4 base UMI at their 3' end, and with a 3 base barcode at the 3' end. For example:

```
@EWSim-1.1-umi5-reada-umix-bar0.0
AAAAATGGCATCCACCGATTTCTCCAAGATTGAACGTAACG
+
IHIHIIIHHIIHHIHIIHHHIIHIHIIIHIIIHHIIIIHIH
@EWSim-2.1-umi5-reada-umiy-bar1.0
AAAAATGGCATCCACCGATTTCTCCAAGATTGAAATATGAC
+
?I?DDIECDEIA?GBGEGGB?EG@GDFEG?GDB?A?IDIDA
@EWSim-3.1-umi5-readb-umix-bar2.0
AAAATCTAGATTAGAAAGATTGACCTCATTAACGTACGA
+
?I?CGGIFCICHHCFFDFB@DDGA??BFDDHABBIF@GC
```

To extract these UMIs, and also the barcode, we can use the regular expression:

```
^(?P<umi_1>.{4}).+(?P<umi_2>.{4})(?P<cell_1>.{3})$
```

To break this down:

* `^(?P<umi_1>.{4})` matches the 4 bases at the start of the read and extracts this as a UMI, named `<umi_1>`.
* `.+` matches all the bases between the first 4 bases and last 7 bases (`+` indicates that there must be at least one).
* `(?P<umi_2>.{4})$` matches the first 4 of the 7 bases at the end of the read and extracts this as a UMI, named `<umi_2>`.
* `(?P<cell_1>.{4})$` matches the 3 bases at the end of the read and extracts this as a barcode, named `<cell_1>`.

Applying this to each example read in our multiplexed FASTQ file would give the following values:

```
@EWSim-1.1-umi5-reada-umix-bar0.0
AAAAATGGCATCCACCGATTTCTCCAAGATTGAACGTAACG # Read
AAAA                                      # <umi_1>
    ATGGCATCCACCGATTTCTCCAAGATTGAA        # Read-post UMI extraction
                                  CGTA    # <umi_2>
				      ACG # <cell_1>
@EWSim-2.1-umi5-reada-umiy-bar1.0
AAAAATGGCATCCACCGATTTCTCCAAGATTGAAATATGAC # Read
AAAA                                      # <umi_1>
    ATGGCATCCACCGATTTCTCCAAGATTGAA        # Read-post UMI extraction
                                  ATAT    # <umi_2>
                                      GAC # <cell_1>
@EWSim-3.1-umi5-readb-umix-bar2.0
AAAATCTAGATTAGAAAGATTGACCTCATTAACGTACGA   # Read
AAAA                                      # <umi_1>
    TCTAGATTAGAAAGATTGACCTCATTAA        # Read-post UMI extraction
                                CGTA      # <umi_2>
                                    CGA   # <cell_1>
```

If we were to use `umi_tools extract` to extract these UMIs and barcode, using this regex, our resulting FASTQ file would look like the following:

```
@EWSim-1.1-umi5-reada-umix-bar0.0_ACG_AAAACGTA
ATGGCATCCACCGATTTCTCCAAGATTGAA
+
IIIHHIIHHIHIIHHHIIHIHIIIHIIIHH
@EWSim-2.1-umi5-reada-umiy-bar1.0_GAC_AAAAATAT
ATGGCATCCACCGATTTCTCCAAGATTGAA
+
DIECDEIA?GBGEGGB?EG@GDFEG?GDB?
@EWSim-3.1-umi5-readb-umix-bar2.0_CGA_AAAACGTA
TCTAGATTAGAAAGATTGACCTCATTAA
+
GGIFCICHHCFFDFB@DDGA??BFDDHA
```

The extracted UMIs are appended together (`<umi_1><umi2>`) and are then appended to the barcode, with a delimiter, `_`. Together the barcode and UMIs are then appended to the end of the read header, again with a delimiter `_`.


 ### Example: the inDrop barcode read
Read 1 from the inDrop technique of Klein *et al* consists of a two part cell barcode separated by an a 22 base adapter sequence. The first part of the cell barcode can be between 8 and 12 bases in length, and the second part is always 8 bases long. A 6 base UMI follows the second part of the cell barcode, and this is then followed by at least three T bases. 

This gives us the regex:
`(?P<cell_1>.{8,12})(?P<discard_1>GAGTGATTGCTTGTGACGCCTT){s<=2}(?P<cell_2>.{8})(?P<umi_1>.{6})T{3}.*`

The first part of the cell barcode is matched by `.{8,12}`. We wish to capture this and use it as the first part of the cell barcode, so we name this group `cell_1`. This is followed by the adapter sequence `GAGTGATTGCTTGTGACGCCTT` which is captured in a group `discard_1` so that it is discarded. We wish to allow up to two mismatches when matching the adapter sequence, so we follow the group with `{s<=2}` (this is a syntax specific to the  [`regex`](https://pypi.org/project/regex/) library). Next up is the second part of the cell barcode, 8 anythings in a group called `cell_2`: `(?<cell_2>.{8})`, the final cell barcode attached to the read header will consist of the concatenation of these two parts. `(?<umi_1?.{6})` captures the six base UMI. Finally three Ts and then any number of any base is matched with `T{3}.*`. Because these are not in a capture group they are left on the read and not discarded or moved to the read header. 

### Regular what-nows?
Regular expressions (abbreviated regex) are a language used for flexible string matching. We will not go into the full details of this here, there are many tutorials available on the web. A good reference manual is the python regular expression page. 

The basics you will need are that patterns consist of characters that must be matched. For for example `A` matches the character A. If we put sets of characters in square brackets, we match any of the characters, so the pattern `[AT]` matches an A or a T. We also have wildcards. There are many wildcards, but the most useful is `.` which matches anything. 

We can also repeat things. `*` and `+` mean zero-or-more and one-or-more respectively, so `A*` matches zero-or-more A characters and `.+` matches or or more of any character. We can be more restrictive using curly braces: `A{3}` matches exactly 3 A characters and `A{1,3}` matches 1 to 3 A characters.

Sometimes we want to "capture" part of the match to the pattern so that we can tell what it was. For example the pattern `.{4}TTTT` matches any four characters and then four Ts. After testing if a string matches a pattern we might want to do what the four characters were. We do this by enclosing the bit we are interested in round brackets - this is a so called "capture group". So if we test the pattern `(.{4})TTTT` against the string `ATCGTTTT` it will match and the content of the group will be `ATCG`. Sometimes we want to capture more than one group, so if we search `(.{4})TTTT(.{4})` we have two groups (numbered from the left). If we search against `ATCGTTTTAAAA`, when we fetch group 1, if will contain `ATCG` and when we fetch group 2 it will contain `AAAA`.
 
&nbsp;

