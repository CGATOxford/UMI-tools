Tools for dealing with Unique Molecular Identifiers
====================================================

This repository contains a number of tools for dealing with Unique Molecular Identifiers (UMIs)/Random Molecular Tags (RMTs). Currently there are two tools:

extract_umi.py:   Flexible removal of UMI sequences from fastq reads. UMIs are removed and appended
                  to the read name. Any other barcode, for example a library barcode, is left on the
                  read.

dedup_umi.py:     Implements a number of different UMI deduplication schemes. The recommended methods
                  are `directional_adjecency` and `adjecency`. In general `directional_adjecency` seems
                  to be less sensitive to starting conditions, but there are situations where `adjecency`
                  might out perform. See simulation results at the  `CGAT <https://cgatoxford.wordpress.com/2015/08/14/unique-molecular-identifiers-the-problem-the-solution-and-the-proof/>`_.

Instalation
------------

Both tools are just python scripts. Type

```
python dedup_umi.py --help
```

or


```
python extract_umi.py --help
```

for help. `dedup_umi.py` is dependent on `numpy`, `pandas` and both are dependent, at the moment, on `CGAT <https://www.cgat.org/downloads/public/cgat/documentation/cgat.html#cgat>`_.
