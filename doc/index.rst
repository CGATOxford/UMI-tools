.. UMI-tools documentation master file, created by
   sphinx-quickstart on Sun Feb  3 18:27:13 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: https://images1-focus-opensocial.googleusercontent.com/gadgets/proxy?url=https://cloud.githubusercontent.com/assets/6096414/19521726/a4dea98e-960c-11e6-806a-a18ff04a391e.png&container=focus&resize_w=550

Tools for dealing with Unique Molecular Identifiers
====================================================

.. note::

   **Important update**: We now recommend the use of `alevin`
   for 10x Chromium and Drop-seq droplet-based scRNA-Seq. `alevin` is an
   accurate, fast and convenient end-to-end tool to go from fastq ->
   count matrix and extends the UMI error correction in `UMI-tools`
   within a framework that also enables quantification of droplet
   scRNA-Seq without discarding multi-mapped reads.  See `alevin
   documentation
   <https://salmon.readthedocs.io/en/latest/alevin.html>`_ and `alevin
   pre-print <https://www.biorxiv.org/content/10.1101/335000v2>`_ for
   more information

Welcome to the UMI-tools documentation. UMI-tools contains tools for
dealing with Unique Molecular Identifiers (UMIs)/Random Molecular Tags
(RMTs) and single cell RNA-Seq cell barcodes.

UMI-tools was published in `Genome Research <http://genome.cshlp.org/content/early/2017/01/18/gr.209601.116.abstract>`_ on 18 Jan '17 (open access)


.. toctree::
   :maxdepth: 1
   :caption: Contents:

   INSTALL
   QUICK_START
   the_methods
   Single_cell_tutorial
   regex
   faq
   API
   common_options
   release_notes

Tools
-----
   
Currently there are 6 commands. The ``extract`` and ``whitelist``
commands are used to prepare a fastq containg UMIs +/- cell barcodes
for alignment.

* whitelist:
   **Builds a whitelist of the 'real' cell barcodes**
      This is useful for droplet-based single cell RNA-Seq where the
      identity of the true cell barcodes is unknown. The whitelist can
      then be used to filter cell barcodes with extract (see below)

* extract:
   **Flexible removal of UMI sequences from fastq reads.**
      UMIs are removed and appended to the read name. Any other
      barcode, for example a library barcode, is left on the read. Can
      also filter reads by quality or against a whitelist (see above)

The remaining commands, ``group``, ``dedup`` and
``count``/``count_tab``, are used to identify PCR duplicates using the
UMIs and perform different levels of analysis depending on the needs
of the user. A number of different UMI deduplication schemes are
enabled - The recommended method is *directional*. For more deails about the
deduplication schemes see :doc:`the_methods`

* dedup:
   **Groups PCR duplicates and deduplicates reads to yield one read per group**
      Use this when you want to remove the PCR duplicates prior to any
      downstream analysis

* group: 
   **Groups PCR duplicates using the same methods available through `dedup`.**
      This is useful when you want to manually interrogate the PCR
      duplicates or perform bespoke downstream processing such as
      generating consensus sequences
   
* count:
   **Groups and deduplicates PCR duplicates and counts the unique molecules per gene**
      Use this when you want to obtain a matrix with unique molecules
      per gene, per cell, for scRNA-Seq

* count_tab:
   **As per count except input is a flatfile**

.. toctree::
   :maxdepth: 1
   :caption: Commands:

   whitelist <reference/whitelist>
   extract <reference/extract>
   group <reference/group>
   dedup <reference/dedup>
   count <reference/count>
   count_tab <reference/count_tab>

Each tool has a set of :doc:`common_options` for input/output,
profiling and debugging.

See :doc:`QUICK_START` for a quick tutorial on
the most common usage pattern.

If you want to use UMI-tools in single-cell RNA-Seq data processing,
see :doc:`Single_cell_tutorial`

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

For more detail see :doc:`INSTALL`


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

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
