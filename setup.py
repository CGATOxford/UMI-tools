from __future__ import absolute_import

import sys
import os
import glob

from ez_setup import use_setuptools
use_setuptools("10.0")
import setuptools

from umi_tools import __version__
from setuptools import setup, find_packages, Extension

from distutils.version import LooseVersion
if LooseVersion(setuptools.__version__) < LooseVersion('1.1'):
    print ("Version detected:", LooseVersion(setuptools.__version__))
    raise ImportError(
        "umi_tools requires setuptools 1.1 higher")

###############################################################
###############################################################
# Define dependencies 
# Perform a umi_tools Installation

major, minor1, minor2, s, tmp = sys.version_info

if (major == 2 and minor1 < 7) or major < 2:
    raise SystemExit("""umi_tools requires Python 2.7 or later.""")

umi_tools_packages = ["umi_tools"]
umi_tools_package_dirs = {'umi_tools': 'umi_tools'}

# debugging pip installation
#install_requires = []
#for requirement in (
#        l.strip() for l in open('requirements.txt') if not l.startswith("#")):
#    install_requires.append(requirement)

install_requires = [
    "setuptools>=1.1",
    "numpy>=1.7",
    "pandas>=0.12.0",
    "future",
    "regex",
    "scipy",
    "matplotlib",
    "pybktree"]

# This is a hack. When Pysam is installed from source, the recorded
# version is 0.2.3, even though a more recent version is actaully
# installed. In the following, if pysam is not detected, pysam will be
# install, presumably this will be the lastest version. If pysam is
# present detect its version with pysam.__version__.  The only problem
# with this is that if pysam is present, but out of date, the system
# will not recognise the update

try:
    import pysam
    if LooseVersion(pysam.__version__) < LooseVersion('0.8.4'):
        print("""
        
    ######################################################################
    #
    # WARNING:
    # Pysam is installed, but not recent enough. We will update pysam, but
    # the system may fail to detect that pysam has been updated. If this
    # happens please run setup again"
    #
    ######################################################################

        """)

        install_requires.append("pysam>=0.8.4")

except ImportError:
    install_requires.append("pysam")

##########################################################
##########################################################
# Classifiers
classifiers = """
Development Status :: 3 - Alpha
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved
Programming Language :: Python
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
"""

setup(
    # package information
    name='umi_tools',
    version=__version__,
    description='umi_tools: Tools for UMI analyses',
    author='Ian Sudbery',
    author_email='i.sudbery@sheffield.ac.uk',
    license="MIT",
    platforms=["any"],
    keywords="computational genomics",
    long_description='umi_tools: Tools for UMI analyses',
    classifiers=list(filter(None, classifiers.split("\n"))),
    url="https://github.com/CGATOxford/UMI-tools",
    download_url="https://github.com/CGATOxford/UMI-tools/tarball/%s" % __version__,
    # package contents
    packages=umi_tools_packages,
    package_dir=umi_tools_package_dirs,
    include_package_data=True,
    # dependencies
    #setup_requires=['cython'],
    install_requires=install_requires,
    # extension modules
    ext_modules=[Extension("umi_tools._dedup_umi", ["umi_tools/_dedup_umi.c"])],
    entry_points={
        'console_scripts': ['umi_tools = umi_tools.umi_tools:main']
    },
    # other options
    zip_safe=False,
)
