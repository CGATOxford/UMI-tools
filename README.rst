Tools for dealing with Unique Molecular Identifiers
====================================================

This repository contains tools for dealing with Unique Molecular Identifiers (UMIs)/Random Molecular Tags (RMTs). Currently there are two tools:

* extract: Flexible removal of UMI sequences from fastq reads.
    UMIs are removed and appended to the read name. Any other barcode, for example a library barcode, is left on the read.

* dedup: Implements a number of different UMI deduplication schemes.
    The recommended method is `directional_adjecency`.

See simulation results at the `CGAT blog <https://cgatoxford.wordpress.com/2015/08/14/unique-molecular-identifiers-the-problem-the-solution-and-the-proof/>`_.

`Genome Science 2015 poster <http://f1000research.com/posters/4-728>`_.

`Biorxiv Preprint <http://dx.doi.org/10.1101/051755>`_.



Installation
------------

If you're using Conda, you can use:

.. code:: bash

   conda install -c https://conda.anaconda.org/toms umi_tools

Or pip:

.. code:: bash

   pip install umi_tools


Or if you'd like to work directly from the git repository:

.. code:: bash

   git clone git@github.com:CGATOxford/UMI-tools.git

Enter repository and run:

.. code:: bash

   python setup.py install


Help
----- 

To get help on umi_tools run

.. code:: bash

   `umi_tools --help`

To get help on umi_tools extract run

.. code:: bash

   `umi_tools extract --help`

To get help on umi_tools dedup run

.. code:: bash

   `umi_tools dedup --help`


Dependencies
------------
umi_tools is dependent on `numpy`, `pandas`, `cython`, `pysam` and `future`
