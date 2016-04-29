import sys
import os
import glob

try:
    import Cython
except ImportError:
    raise ImportError(
        "UMI-tools requires cython to "
        "be installed before running setup.py (pip install cython)")

try:
    import pysam
except ImportError:
    raise ImportError(
        "UMI-tools requires pysam to "
        "be installed before running setup.py (pip install pysam)")

########################################################################
########################################################################
# Import setuptools
# Use existing setuptools, otherwise try ez_setup.

try:
    import setuptools
except ImportError:
    raise ImportError(
        "UMI-tools requires setuptools"
        "be installed before running setup.py (pip install setuptools)")

from setuptools import setup, find_packages, Extension

from distutils.version import LooseVersion
if LooseVersion(setuptools.__version__) < LooseVersion('1.1'):
    print ("Version detected:", LooseVersion(setuptools.__version__))
    raise ImportError(
        "UMI-tools requires setuptools 1.1 higher")

from Cython.Build import cythonize
########################################################################
########################################################################
# collect umi_tools version
sys.path.insert(0, "umi_tools")
import version

version = version.__version__

###############################################################
###############################################################
# Define dependencies 
# Perform a umi_tools Installation

major, minor1, minor2, s, tmp = sys.version_info

if (major == 2 and minor1 < 7) or major < 2:
    raise SystemExit("""UMI-tools requires Python 2.7 or later.""")

umi_tools_packages = ["umi_tools"]
umi_tools_package_dirs = {'umi_tools': 'umi_tools'}

install_requires = []

for requirement in (
        l.strip() for l in open('requirements.txt') if not l.startswith("#")):
    install_requires.append(requirement)


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
    version=version,
    description='umi-tools: Tools for UMI analyses',
    author='Ian Sudbery',
    author_email='i.sudbery@sheffield.ac.uk',
    license="MIT",
    platforms=["any"],
    keywords="computational genomics",
    long_description='umi-tools: Tools for UMI analyses',
    classifiers=list(filter(None, classifiers.split("\n"))),
    url="https://github.com/CGATOxford/UMI-tools",
    download_url="https://github.com/CGATOxford/UMI-tools/tarball/%s" % version,
    # package contents
    packages=umi_tools_packages,
    package_dir=umi_tools_package_dirs,
    include_package_data=True,
    # dependencies
    install_requires=install_requires,
    # extension modules
    ext_modules=cythonize("umi_tools/_dedup_umi.pyx"),
    entry_points={
        'console_scripts': ['umi_tools = umi_tools.umi_tools:main']
    },
    # other options
    zip_safe=False,
)
