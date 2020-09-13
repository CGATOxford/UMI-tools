.. image:: https://user-images.githubusercontent.com/6096414/93030687-c7cf7300-f61c-11ea-92b8-102ec17ef6aa.png

UMI-tools was published in `Genome Research <http://genome.cshlp.org/content/early/2017/01/18/gr.209601.116.abstract>`_ on 18 Jan '17 (open access)

For full documentation see https://umi-tools.readthedocs.io/en/latest/

Tools for dealing with Unique Molecular Identifiers
====================================================

This repository contains tools for dealing with Unique Molecular
Identifiers (UMIs)/Random Molecular Tags (RMTs) and single cell
RNA-Seq cell barcodes. Currently there are 6
commands. 

The ``extract`` and ``whitelist`` commands are used to prepare a
fastq containg UMIs +/- cell barcodes for alignment. 

* whitelist:
   **Builds a whitelist of the 'real' cell barcodes**
      This is useful for droplet-based single cell RNA-Seq where the
      identity of the true cell barcodes is unknown. Whitelist can
      then be used to filter with extract (see below)

* extract:
   **Flexible removal of UMI sequences from fastq reads.**
      UMIs are removed and appended to the read name. Any other
      barcode, for example a library barcode, is left on the read. Can
      also filter reads by quality or against a whitelist (see above)

The remaining commands, ``group``, ``dedup`` and ``count``/``count_tab``, are used to
identify PCR duplicates using the UMIs and perform different levels of
analysis depending on the needs of the user. A number of different UMI
deduplication schemes are enabled - The recommended method is
*directional*.

* dedup:
   **Groups PCR duplicates and deduplicates reads to yield one read per group**
      Use this when you want to remove the PCR duplicates prior to any
      downstream analysis

* group: 
   **Groups PCR duplicates using the same methods available through `dedup`.**
      This is useful when you want to manually interrogate the PCR duplicates
   
* count:
   **Groups and deduplicates PCR duplicates and counts the unique molecules per gene**
      Use this when you want to obtain a matrix with unique molecules
      per gene, per cell, for scRNA-Seq.

* count_tab:
   **As per count except input is a flatfile**

See `QUICK_START.md <./doc/QUICK_START.md>`_ for a quick tutorial on
the most common usage pattern.

If you want to use UMI-tools in single-cell RNA-Seq data processing,
see `Single_cell_tutorial.md <./doc/Single_cell_tutorial.md>`_

**Important update**: We now recommend the use of `alevin` for droplet-based
scRNA-Seq (e.g 10X, inDrop etc). `alevin` is an accurate, fast and convenient end-to-end tool to go from fastq -> count matrix and  extends the UMI error correction in `UMI-tools` within a framework that also enables quantification of droplet scRNA-Seq without discarding multi-mapped reads.  See `alevin documentation <https://salmon.readthedocs.io/en/latest/alevin.html>`_ and `alevin pre-print <https://www.biorxiv.org/content/10.1101/335000v2>`_ for more information

The ``dedup``, ``group``, and ``count`` / ``count_tab`` commands make use of network-based methods to resolve similar UMIs with the same alignment coordinates. For a background regarding these methods see:

`Genome Research Publication <http://genome.cshlp.org/content/early/2017/01/18/gr.209601.116.abstract>`_

`Blog post discussing network-based methods <https://cgatoxford.wordpress.com/2015/08/14/unique-molecular-identifiers-the-problem-the-solution-and-the-proof/>`_.


Installation
------------

If you're using Conda, you can use:

.. code:: bash

   $ conda install -c bioconda -c conda-forge umi_tools

Or pip:

.. code:: bash

   $ pip install umi_tools


Or if you'd like to work directly from the git repository:

.. code:: bash

   $ git clone https://github.com/CGATOxford/UMI-tools.git

Enter repository and run:

.. code:: bash

   $ python setup.py install

For more detail see `INSTALL.rst <./doc/INSTALL.rst>`_

Help
----- 

For full documentation see https://umi-tools.readthedocs.io/en/latest/

See `QUICK_START.md <./doc/QUICK_START.md>`_ and
`Single_cell_tutorial.md <./doc/Single_cell_tutorial.md>`_ for tutorials on the most common usage patterns.

To get help on umi_tools run

.. code:: bash

   $ umi_tools --help

To get help on the options for a specific [COMMAND], run

.. code:: bash

   $ umi_tools [COMMAND] --help


Dependencies
------------
umi_tools is dependent on `python>=3.5`, `numpy`, `pandas`, `scipy`, `cython`, `pysam`,
`future`, `regex` and `matplotlib`

