Installation Guide
===================

We currently support installation on Linux (tested on Ubuntu and
CentOS) and Mac OSX. We do not currently support Windows.

There are three possible ways to install UMI-tools: conda, pip or from
source, in decending order of ease.

Conda package manager
----------------------

This is the easiest way to install UMI-tools if you are already using
either anaconda python or miniconda: all depednencies, whether they be
python libraries or system libraries are automatically installed. You
can also do all your installations in seperate isolated "environments"
where installing new software will not affect packages in other
environments. The downside of conda installation is that if you do not
already use anaconda or miniconda, you will be building a completely
new python environment, would have to reinstall all of your libraries
etc. You can read more `about conda here`_.

1. Install miniconda (if not already installed), see `conda
   installation instructions here`_.

2. Type::

   conda install -c https://conda.anaconda.org/toms umi_tools

Thats it, simple as that.


Installation from PyPI using the pip package manager
-----------------------------------------------------

The pip python pacakge manager is the standard package manager. The
advantage over conda is that it is probably already installed on your
system, will use your existing python environments, and plays nicely
with the virtualenv system. On the downside, the installation of
dependencies is not handled as cleanly as in conda. You will need

* python version greater than 2.7
* gcc or compatible c compiler 
* zlib with development headers
* the pip python package manager version at least 1.4

Linux
------

Most systems will already have gcc, pip and zlib installed, so its
worth just trying

    pip install umi_tools

If you get a permissions error try adding `--user` to the pip
command. Note that `umi_tools` will now only be install for the
current user.

If that doesn't work, then you need to find what is missing. You can
check for gcc and pip by typing gcc or pip at a terminal
prompt. Installing GCC and zlib is much easier if you have root access
on your machine. In the unlikely event that you don't have these
installed AND you don't have root access, please speak to someone who
does have root access. Python and pip can be installed without root. 

1.  **Install gcc**: the easiest way is using your package manager. In
    ubuntu or debian::

        sudo apt-get install gcc

    or in CentOS/Redhat/Fedora::

        sudo yum install gcc

2.  **Install zlib**: again, use the package manager. Ubuntu::

        sudo apt-get zlib1g-dev

    or in CentOS/Redhat/Fedora::

        sudo yum install

3.  **Install pip**: pip is also probably available from your package
    manager, but can also be installed without this if you don't have
    root access::

        wget https://bootstrap.pypa.io/get-pip.py
        python get-pip.py --user

The pip command at the top should now work. 


Apple OS-X
-----------

The good news is that `zlib` is installed by default of OS X and
modern versions include an upto date `python`. The bad news is that
`gcc` and `pip` are generally not included (although many users many
have installed them already).
