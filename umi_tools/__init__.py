from __future__ import absolute_import

from umi_tools.version import __version__

# Import happens on setup, but on setup networks dependencies won't have
# been installed. So if we can't import, don't worry about it.]
try:
    from umi_tools.network import UMIClusterer
except ImportError:
    pass

