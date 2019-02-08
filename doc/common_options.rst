.. _common_options:

==============
Common options
==============

Each of the ``umi_tools`` commands has a set of common options to deal
with input and output files, logging, profiling and debugging


Input/Output options
--------------------

By default, each tool reads from ``stdin`` and outputs to ``stdout``,
with the exception of ``dedup``, ``group`` and ``count``, which cannot
work from ``stdin`` since in some cases they need to parse the input
multiple times. 

By default, logging is sent to ``stdout``. In most cases, you will
want to direct this to a dedicated logfile (``-L``/ ``--log``) or send
the logging to ``stderr`` (``--log2stderr``).

"""""""""""""""
``-I, --stdin``
"""""""""""""""

    File to read stdin from [default = stdin].

""""""""""""""""
``-S, --stdout``
""""""""""""""""

    File where output is to go [default = stdout].

"""""""""""""
``-L, --log``
"""""""""""""

    File with logging information [default = stdout].

""""""""""""""""
``--log2stderr``
""""""""""""""""

    Send logging information to stderr [default = False].

""""""""""""""""""""""""""
``-v, --verbose``
""""""""""""""""""""""""""

    Log level. The higher, the more output [default = 1].

"""""""""""""""
``-E, --error``
"""""""""""""""

    File with error information [default = stderr].

""""""""""""""
``--temp-dir``
""""""""""""""

    Directory for temporary files. If not set, the bash environmental
    variable TMPDIR is used[default = None].

"""""""""""""""""""
``--compresslevel``
"""""""""""""""""""

    Level of Gzip compression to use. Default=6 matches GNU gzip
    rather than python gzip default (which is 9)


profiling and debugging options
-------------------------------

""""""""""""
``--timeit``
""""""""""""

    Store timeing information in file [default=none].

"""""""""""""""""
``--timeit-name``
"""""""""""""""""

    Name in timing file for this class of jobs [default=all].

"""""""""""""""""""
``--timeit-header``
"""""""""""""""""""

    Add header for timing information [default=none].

"""""""""""""""""
``--random-seed``
"""""""""""""""""

    Random seed to initialize number generator with [default=none].
