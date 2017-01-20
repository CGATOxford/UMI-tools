.. image:: https://images1-focus-opensocial.googleusercontent.com/gadgets/proxy?url=https://cloud.githubusercontent.com/assets/6096414/19521726/a4dea98e-960c-11e6-806a-a18ff04a391e.png&container=focus&resize_w=550

UMI-tools was published in `Genome Research <http://genome.cshlp.org/content/early/2017/01/18/gr.209601.116.abstract>`_ on 18 Jan '17 (early access)

Tools for dealing with Unique Molecular Identifiers
====================================================

This repository contains tools for dealing with Unique Molecular Identifiers (UMIs)/Random Molecular Tags (RMTs). Currently there are three tools:

* extract: Flexible removal of UMI sequences from fastq reads.
    UMIs are removed and appended to the read name. Any other barcode, for example a library barcode, is left on the read.

* dedup: Removes PCR duplicates. Implements a number of different UMI deduplication schemes.
    The recommended method is `directional`.
    
* group: Groups PCR duplicates using the same methods available through `dedup`.
    This is useful when you want to interrogate the PCR duplicates

See `QUICK_START.md <QUICK_START.md>`_ for a quick tutorial on the most common usage pattern.

The `dedup` and `group` commands make use of network-based methods to resolve similar UMIs with the same alignment coordinates. For a background regarding these methods see:

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

For more detail see `INSTALL.rst <INSTALL.rst>`_

Help
----- 

See `QUICK_START.md <QUICK_START.md>`_ for a quick tutorial on the most common usage pattern.

To get detailed help on umi_tools run

.. code:: bash

   $ umi_tools --help

To get help on umi_tools extract run

.. code:: bash

   $ umi_tools extract --help

To get help on umi_tools dedup run

.. code:: bash

   $ umi_tools dedup --help


Dependencies
------------
umi_tools is dependent on `numpy`, `pandas`, `cython`, `pysam` and `future`
