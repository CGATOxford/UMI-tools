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


def main():
    argv = sys.argv

    if len(argv) <2 or argv[1] == "--help" or argv[1] == "-h":
        print("For full UMI-tools documentation, see: "
              "https://umi-tools.readthedocs.io/en/latest/\n")
        print(globals()["__doc__"])
        return 0

    if argv[1] == "--version" or argv[1] == "-v":
        print("UMI-tools version: %s" % __version__)
        return 0

    command = argv[1]
    
    try:
        module = importlib.import_module("umi_tools." + command, "umi_tools")
    except ModuleNotFoundError:
        print("'%s' is not a UMI-tools command. See 'umi_tools -h'." % command)
        return 1

    del sys.argv[0]  # remove 'umi-tools' from sys.argv
    return module.main(sys.argv)

if __name__ == "__main__":
    sys.exit(main())
