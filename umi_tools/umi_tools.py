'''
umi_tools.py - Tools for UMI analyses
===============================================

:Author: Tom Smith & Ian Sudbury, CGAT
:Release: $Id$
:Date: |today|
:Tags: Genomics UMI

There are 6 tools:

  - whitelist
  - extract
  - group
  - dedup
  - count
  - count_tab

To get help on a specific tool, type:

    umi_tools <tool> --help

To use a specific tool, type::

    umi_tools <tool> [tool options] [tool arguments]
'''

from __future__ import absolute_import
import os
import sys
import importlib
from umi_tools import __version__


def main(argv=None):

    argv = sys.argv

    path = os.path.abspath(os.path.dirname(__file__))

    if len(argv) == 1 or argv[1] == "--help" or argv[1] == "-h":
        print("For full UMI-tools documentation, see: "
              "https://umi-tools.readthedocs.io/en/latest/\n")
        print(globals()["__doc__"])

        return

    elif len(argv) == 1 or argv[1] == "--version" or argv[1] == "-v":
        print("UMI-tools version: %s" % __version__)

        return

    elif argv[2] in ["--help", "-h", "--help-extended"]:
        print("UMI-Tools: Version %s" % __version__)

    command = argv[1]

    module = importlib.import_module("umi_tools." + command, "umi_tools")
    # remove 'umi-tools' from sys.argv
    del sys.argv[0]
    module.main(sys.argv)

if __name__ == "__main__":
    sys.exit(main())
