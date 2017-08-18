.. image:: https://images1-focus-opensocial.googleusercontent.com/gadgets/proxy?url=https://cloud.githubusercontent.com/assets/6096414/19521726/a4dea98e-960c-11e6-806a-a18ff04a391e.png&container=focus&resize_w=550

UMI-tools was published in `Genome Research <http://genome.cshlp.org/content/early/2017/01/18/gr.209601.116.abstract>`_ on 18 Jan '17 (open access)

Tools for dealing with Unique Molecular Identifiers
====================================================

This repository contains tools for dealing with Unique Molecular
Identifiers (UMIs)/Random Molecular Tags (RMTs) and single Cell
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

The remaining commands, group, dedup and count/count_tab, are used to
identify PCR duplicates using the UMIs and perform different levels of
analysis depending on the needs of the user. A number of different UMI
deduplication schemes are enabled - The recommended method is
_directional_.

* group: 
   **Groups PCR duplicates using the same methods available through `dedup`.**
      This is useful when you want to manually interrogate the PCR duplicates

* dedup:
   **Groups PCR duplicates and deduplicates reads to yield one read per group**
      Use this when you want to remove the PCR duplicates prior to any
      downstream analysis
    
* count:
   **Groups and deduplicates PCR duplicates and counts the unique molecules per gene**
      Use this when you want to obtain a matrix with unique molecules
      per gene. Can also perform per-cell counting for scRNA-Seq.

* count_tab:
   **As per count except input is a flatfile**

See `QUICK_START.md <./docs/QUICK_START.md>`_ for a quick tutorial on
the most common usage pattern.

If you want to use UMI-tools in single-cell RNA-Seq data processing,
see `Single_cell_tutorial.md <./docs/Single_cell_tutorial.md>`_


The ``dedup``, ``group``, and ``count`` / ``count_tab`` commands make use of network-based methods to resolve similar UMIs with the same alignment coordinates. For a background regarding these methods see:

`Genome Research Publication <http://genome.cshlp.org/content/early/2017/01/18/gr.209601.116.abstract>`_

`Blog post discussing network-based methods <https://cgatoxford.wordpress.com/2015/08/14/unique-molecular-identifiers-the-problem-the-solution-and-the-proof/>`_.


Installation
------------

If you're using Conda, you can use:

.. code:: bash

   $ conda install -c https://conda.anaconda.org/toms umi_tools

Or pip:

.. code:: bash

   $ pip install umi_tools


Or if you'd like to work directly from the git repository:

.. code:: bash

   $ git clone https://github.com/CGATOxford/UMI-tools.git

Enter repository and run:

.. code:: bash

   $ python setup.py install

For more detail see `INSTALL.rst <./docs/INSTALL.rst>`_

Help
----- 

See `QUICK_START.md <./docs/QUICK_START.md>`_ and
`Single_cell_tutorial.md <./docs/Single_cell_tutorial.md>`_ for tutorials on the most common usage patterns.

To get detailed help on umi_tools run

.. code:: bash

   $ umi_tools --help

To get help on a specific [COMMAND] run

.. code:: bash

   $ umi_tools [COMMAND] --help


Dependencies
------------
umi_tools is dependent on `numpy`, `pandas`, `scipy`, `cython`, `pysam`,
`future`, `regex` and `matplotlib`
