Installation Guide
===================

.. contents::

We currently support installation on Linux (tested on Ubuntu and
CentOS) and Mac OSX. We do not currently support Windows.

There are three possible ways to install UMI-tools: conda, pip or from
source, in decending order of ease.

.. note::
   As of version 1.0.0 UMI-tools requires python 3.5 or better. If you are
   still using python 2.7, we recommend you switch to python 3. If you
   use `conda` it is possible to have both python 2 and python 3 environments.
   If you can't do this and really need 2.7, we recommend you use UMI-tools
   0.5.5, but note that this will not be updated. 
   
   
Quick Start
-------------

Try one of the following::

    $ conda install -c bioconda -c conda-forge umi_tools

or::

    $ pip install umi_tools

or grab a zip of the `latest release`_ from github and unpack
(replace `0.5.5` with version number; replace ``wget`` with ``curl -O`` for OS X)::

    $ unzip 1.0.0.zip
    $ cd UMI-tools-1.0.0
    $ python setup.py install --user

If these options don't work, see below.

Conda package manager
----------------------

This is the easiest way to install ``UMI-tools`` if you are already using
either anaconda python or miniconda: all depedencies, whether they be
python libraries or system libraries are automatically installed. You
can also do all your installations in seperate isolated "environments"
where installing new software will not affect packages in other
environments. The downside of ``conda`` installation is that if you do not
already use anaconda or miniconda, you will be building a completely
new python environment, would have to reinstall all of your libraries
etc. You can read more `about conda here`_.

1. Install miniconda (if not already installed), see `conda
   installation instructions here`_.

2. Type::

    $ conda install -c bioconda -c conda-forge umi_tools

That's it, simple as that.


Installation from PyPI using the pip package manager
-----------------------------------------------------

The pip python pacakge manager is the standard package manager. The
advantage over conda is that it is probably already installed on your
system, will use your existing python environments, and plays nicely
with the virtualenv system. On the downside, the installation of
dependencies is not handled as cleanly as in conda. You will need

* python version greater than 3.5
* gcc or compatible c compiler 
* zlib with development headers
* the pip python package manager version at least 1.4

Linux
++++++

Most systems will already have ``gcc``, ``pip`` and ``zlib`` installed, so its
worth just trying::

    $ pip install umi_tools

If you also have python 2 on your system, you may need to use ``pip3``
rather than ``pip``. If you get a permissions error try adding ``--user``
to the ``pip`` command. Note that ``umi_tools`` will now only be 
installed for the current user.

If that doesn't work, then you need to find what is missing. You can
check for gcc and pip by typing gcc or pip at a terminal
prompt. Installing GCC and zlib is much easier if you have root access
on your machine. In the unlikely event that you don't have these
installed AND you don't have root access, please speak to someone who
does have root access. Python and pip can be installed without root. 

1.  **Install gcc**: the easiest way is using your package manager. In
    ubuntu or debian::

        $ sudo apt-get install gcc

    or in CentOS/Redhat/Fedora::

        $ sudo yum install gcc

2.  **Install zlib**: again, use the package manager. Ubuntu::

        $ sudo apt-get zlib1g-dev

    or in CentOS/Redhat/Fedora::

        $ sudo yum install zlib-devel

3.  **Install pip**: pip is also probably available from your package
    manager. In ubuntu, Centos, RHEL and fedora the package is called
    `python-pip`. In CentOS/RHEL the package is located in the EPEL
    repository which needs to be installed first. You could also
    install pip from the web::
    
        $ wget https://bootstrap.pypa.io/get-pip.py
        $ python get-pip.py --user

    but in this case you'll need to make sure that the ``python-dev``
    (Ubuntu) or ``python-devel`` (CentOS/RHEL/fedora) packages are
    installed.

The pip command at the top should now work. 


Apple OS X
+++++++++++

The good news is that `zlib` is installed by default of OS X. The
bad news is that `gcc` and `pip` are generally not included (although
many users may have installed them already). Furthermore, it's generally
not advisable to use the default python since installation of third party
python libraries leads to difficulties with permissions, especially since the
introduction of System Integrity Protection (SIP) from OS X El Capitan onwards.
For this reason, we recommend using a non-default python. 

If you only have the default python (e.g /usr/local/bin/python) there are a number of ways
to install another instance of python. Many OS X users recommend using the ``homebrew``
package manager to manage command line packages on OS X. You can find `instructions here`_
for installation python via ``homebrew``. This will also install setuptools and pip.
You can install gcc via homebrew by following `these instructions`_::
    
    $ brew install gcc48

You may also need to install ``freetype``::

    $ brew install freetype


**Install UMI-tools**: You should now have everything you need to
install ``UMI-tools``::

        $ pip install umi_tools

We have had reports that the current version of one of the
``UMI-tools`` dependencies, ``pysam``, is causing problems on the latest
versions of OS X. If your installation is failing on the
installation of pysam, try forcing an older version with::

        $ pip install pysam==0.8.4

before installing ``umi_tools``.

If you don't want to do use homebrew, here are non-homebrew instructions for installing gcc and pip as needed:

1.  **Install gcc**: Apples XCode suite includes ``gcc``. Installation depends
    on which version of OS X you are using

    - *Mac OS X 10.9* or higher: Open a terminal and run::

        $ xcode-select --install

    - *Mac OS X 10.8* or lower: go to Apple's `developer download
      page`_ and download Command Line Tools for XCode. You'll need a
      developer account.

2.  **Install pip**: In a terminal type::

        $ curl -O https://bootstrap.pypa.io/get-pip.py
        $ python get-pip.py


Installing from source
-----------------------

There are several reaons you might want to install from source. If for
example you need to install the most up-to-date version, or if you
can't or don't want to use one of the package managers above. There
are two levels of installing from source. The first is to install the
dependencies using one of the pacakge managers above, and then just
install ``umi_tools`` from source. The second is to install everything
from source without the help of pip or conda.


Depedencies from conda/PyPI manager
++++++++++++++++++++++++++++++++++++

1.  Download the UMI-tools code, either the `latest release`_ or the
    `master branch`_ (which should contain the lastest development
    version) and unpack the zip or tar and enter the directory::

        $ unzip 1.0.0.zip
        $ cd UMI-tools-1.0.0

    or clone the repository::

        $ git clone https://github.com/CGATOxford/UMI-tools.git

3.  Use your python package manager to install the
    dependencies. e.g. for ``pip``

        $ pip install -r requirements.txt

    or with ``conda``::

        $ conda install setuptools
        $ conda install pandas
        $ conda install future
        $ conda install scipy
        $ conda install matplotlib
        $ conda config --add channels bioconda
        $ conda install regex
        $ conda install pysam

4.  Install UMI-tools using the ``setup.py`` script::

        $ python setup.py install --user

Completely from source
+++++++++++++++++++++++

.. WARNING::
    **This section is deprecated and no longer updateed**. Once upon a time it
    was possible for us to provide complete instructions for installing completely 
    from source without a package manager. Unfortunately, our dependencies have 
    multiplied and the dependencies of our dependencies have also multiplied. 
    You can try the below and it may work as the system libraries required are not
    particularly rare, especially if you are already doing bioinformatics. However, 
    if one of the dependencies fails to install, I'm afraid you are on your own. 

This method will allow you to install without installing pip or
conda. It is in theory possible to install completely without root by
installing gcc, zlib and python-dev in your home directory, but that
is beyond the scope of this document. You are also going to need a ``g++``
compatiable compiler. On OS X ``XCode`` has one of these by default. On
Linux install the ``build-essential`` or ``g++`` packages.

1.  Download and install `Cython`. For OS X replace ``wget`` with ``curl
    -O``::

       $  wget https://pypi.python.org/packages/c6/fe/97319581905de40f1be7015a0ea1bd336a756f6249914b148a17eefa75dc/Cython-0.24.1.tar.gz
        $ tar -xzf Cython-0.24.1.tar.gz
        $ cd Cython-0.24.1.tar.gz
        $ python setup.py install --user

2.  Download and install ``UMI-tools``::

        $ wget https://github.com/CGATOxford/UMI-tools/archive/master.zip
        $ unzip master.zip
        $ cd UMI-tools-master
        $ python setup.py install --user

    running this is probably going to take quite a long time. You will
    probably see quite a lot of warning messages that look like
    errors. 

    The most likely fail point is installing ``pysam``. Due to a bug in 
    pysam, when it is installed from source, the recorded install version
    is wrong. Thus, if you get the error::

        $ pysam 0.2.3 is installed by 0.8.4 is required by umi_tools

    try just running setup again. 

    In addition, as we pointed out above, we have had reports that 
    installation of the lastest ``pysam`` fails on the latest OS X. If
    this is the case, try installing an older version of ``pysam``::

        $ curl -O https://pypi.python.org/packages/27/89/bf8c44d0bfe9d0cadab062893806994c168c9f490f67370fc56d6e8ba224/pysam-0.8.4.tar.gz
        $ tar -xzf pysam-0.8.4.tar.gz
        $ cd pysam-0.8.4
        $ python setup.py install --user

Running tests
+++++++++++++

After installing from source you can run the test suite to make sure everything is working. To do this you'll need to install `nose` and `pyyaml` using your favourite package manager and then run::

    $ nosetests tests/test_umi_tools.py
    

Getting further help
---------------------

If you are still having trouble with installation, contact us by by
creating an issue on our `github issues page`_.

.. _about conda here: http://conda.pydata.org/docs/intro.html
.. _conda installation instructions here: http://conda.pydata.org/docs/installation.html
.. _developer download page: https://developer.apple.com/downloads/index.action#
.. _latest release: https://github.com/CGATOxford/UMI-tools/releases/latest
.. _master branch: https://github.com/CGATOxford/UMI-tools/archive/master.zip
.. _github issues page: https://github.com/CGATOxford/UMI-tools/issues/new
.. _instructions here: http://docs.python-guide.org/en/latest/starting/install/osx/
.. _these instructions: http://www-scf.usc.edu/~csci104/installation/gccmac.html
