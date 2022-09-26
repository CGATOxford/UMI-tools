'''
umi_tools.py - Tools for UMI analyses
=====================================

:Author: Tom Smith & Ian Sudbury, CGAT
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

    elif len(argv) > 2 and  argv[2] in ["--help", "-h", "--help-extended"]:
        print("UMI-Tools: Version %s" % __version__)

        return 0

    command = argv[1]

    try:
        module = importlib.import_module("umi_tools." + command, "umi_tools")
    except ImportError:
        print("'%s' is not a UMI-tools command. See 'umi_tools -h'.\n" % command)
        print("For full UMI-tools documentation, see: "
              "https://umi-tools.readthedocs.io/en/latest/\n")
        print(globals()["__doc__"])

        return 1

    # remove 'umi-tools' from sys.argv
    del sys.argv[0]
    module.main(sys.argv)


if __name__ == "__main__":
    sys.exit(main())
