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

from Cython.Distutils import build_ext

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

# Need to check Python 3 compatibility
if major == 3:
    raise SystemExit("""UMI-tools is not fully python3 compatible""")

if (major == 2 and minor1 < 7) or major < 2:
    raise SystemExit("""UMI-tools requires Python 2.7 or later.""")

# use requires.txt to identify requirements
#requires = []
#for requirement in (
#        l.strip() for l in open('requires.txt') if not l.startswith("#")):
#    requires.append(requirement)

umi_tools_packages = ["umi_tools"]
umi_tools_package_dirs = {'umi_tools': 'umi_tools'}

# automatically build pyximport script extensions
pyx_files = glob.glob("umi_tools/*.pyx")
tool_extensions = []
pysam_dirname = os.path.dirname(pysam.__file__)
include_dirs = [pysam.get_include()]

extra_link_args = [os.path.join(pysam_dirname, x) for x in 
                   pysam.get_libraries()]

for pyx_file in pyx_files:
    script_name = os.path.basename(pyx_file)
    script_prefix = script_name[:-4]
    tool_extensions.append(
        Extension("umi_tools.%s" % (script_prefix),
                  sources=[pyx_file],
                  extra_link_args=extra_link_args,
                  include_dirs=include_dirs,
                  define_macros=pysam.get_defines())
    )

ext_modules = tool_extensions

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
    author='Ian Sudbury',
    author_email='i.sudbery@sheffield.ac.uk',
    license="BSD",
    platforms=["any"],
    keywords="computational genomics",
    long_description='umi-tools: Tools for UMI analyses',
    classifiers=filter(None, classifiers.split("\n")),
    url="https://github.com/CGATOxford/UMI-tools",
    # package contents
    packages=umi_tools_packages,
    package_dir=umi_tools_package_dirs,
    include_package_data=True,
    # dependencies
    install_requires=["setuptools>=1.1",
                      "cython>=0.19",
                      "numpy>=1.7",
                      "pandas>=0.12.0",
                      "pysam>=0.8.4",
                      "future"],
    #cmdclass={'build_ext': build_ext},
    entry_points={
        'console_scripts': ['umi_tools = umi_tools.umi_tools:main']
    },
    # other options
    zip_safe=False,
)
