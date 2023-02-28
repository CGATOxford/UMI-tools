"""test_scripts.py
==================
nose test script for umi_tools
Test data for scripts are contained in the directory umi_tools
described in the :file:`tests.yaml`

This code is modified with permission from:
https://github.com/CGATOxford/cgat/blob/master/tests/test_scripts.py
"""

from __future__ import print_function

import subprocess
import tempfile
import os
import shutil
import re
import glob
import gzip
import yaml
import time
import hashlib
import sys
import platform

from nose.tools import ok_

PYTHON_VERSION = platform.python_version()
IS_PY3 = sys.version_info.major >= 3


SUBDIRS = ("gpipe", "optic")

# Setup logging
LOGFILE = open("test_scripts.log", "a")
DEBUG = os.environ.get("CGAT_DEBUG", False)


def check_main(script):
    '''test is if a script can be imported and has a main function.
    '''

    # The following tries importing, but I ran into problems - thus simply
    # do a textual check for now
    # path, basename = os.path.split(script)
    # pyxfile = os.path.join( path, "_") + basename + "x"
    # ignore script with pyximport for now, something does not work
    # if not os.path.exists( pyxfile ):
    #     with warnings.catch_warnings() as w:
    #         warnings.simplefilter("ignore")
    #         (file, pathname, description) =
    #                imp.find_module( basename[:-3], [path,])
    #         module = imp.load_module( basename, file, pathname, description)
    #     ok_( "main" in dir(module), "no main function" )

    # subsitute gpipe and other subdirectories.
    for s in SUBDIRS:
        script = re.sub("%s_" % s, "%s/" % s, script)

    # check for text match
    ok_([x for x in open(script) if x.startswith("def main(")],
        "no main function")


def compute_checksum(filename):
    '''return md5 checksum of file.'''
    return hashlib.md5(open(filename, 'rb').read()).hexdigest()

#########################################
# List of tests to perform.
#########################################
# The fields are:


def check_script(test_name,
                 stdin,
                 options, outputs,
                 references,
                 current_dir,
                 sort=False):
    '''check script.
    # 1. Name of the script
    # 2. Filename to use as stdin
    # 3. Option string
    # 4. List of output files to collect
    # 5. List of reference files

    '''
    working_dir = "tests"

    tmpdir = tempfile.mkdtemp()

    t1 = time.time()

    stdout = os.path.join(tmpdir, 'stdout')

    if stdin:
        stdin = '--stdin=%s' % (os.path.join(os.path.abspath(working_dir), stdin))
    else:
        stdin = ""

    if options:
        options = re.sub("%TMP%", tmpdir, options)
        options = re.sub("<TMP>", tmpdir, options)
        options = re.sub("%DIR%", os.path.abspath(working_dir), options)
        options = re.sub("<DIR>", os.path.abspath(working_dir), options)
    else:
        options = ""

    options = re.sub("\n", "", options)

    # use /bin/bash in order to enable "<( )" syntax in shells
    statement = ("/bin/bash -c "
                 "'umi_tools %(options)s %(stdin)s > %(stdout)s'") % locals()

    process = subprocess.Popen(statement,
                               shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               cwd=tmpdir)

    if DEBUG:
        print("tmpdir={}".format(tmpdir), end=" ")

    process_stdout, process_stderr = process.communicate()

    fail = False
    msg = ""

    if process.returncode != 0:
        fail = True
        msg = "error in statement: %s; stderr=%s" %\
              (statement, process_stderr)

    if not fail:
        # compare line by line, ignoring comments
        for output, reference in zip(outputs, references):
            if output == "stdout":
                output = stdout
            elif output.startswith("<DIR>/") or \
                    output.startswith("%DIR%/"):
                output = os.path.join(working_dir, output[6:])
            else:
                output = os.path.join(tmpdir, output)

            if not os.path.exists(output):
                fail = True
                msg = "output file '%s'  does not exist: %s" %\
                      (output, statement)

            reference = os.path.join(working_dir, reference)
            if not fail and not os.path.exists(reference):
                fail = True
                msg = "reference file '%s' does not exist (%s): %s" %\
                      (reference, tmpdir, statement)

            if not fail:

                a = _read(output)
                b = _read(reference)

                if sort:
                    a = sorted(a)
                    b = sorted(b)

                if a != b:
                    fail = True
                    msg = ("files %s and %s are not the same\n"
                           "%s\nmd5: output=%i, %s reference=%i, %s") %\
                        (output, reference, statement,
                         len(a),
                         compute_checksum(output),
                         len(b),
                         compute_checksum(reference))

                    diffs = []
                    for aa, bb in zip(a, b):
                        if aa != bb:
                            diffs.append((aa, bb))
                            if len(diffs) > 10:
                                break

                    if sort:
                        msg += ('\n\nNote that files were sorted prior to diff so line numbers '
                                'in diff are with respect to sorted reference '
                                'and output files\n\n')

                    msg += "first 10 differences: {}".format(
                        "\n--\n".join(
                            ["\n".join(map(str, (x)))
                             for x in diffs]))
                    break

    t2 = time.time()
    LOGFILE.write("t%s\t%f\n" % (test_name, t2-t1))
    LOGFILE.flush()

    # preserve coverage information, this gets stored it tmpdir, but
    # needs to be moved to the current directory to be merged.
    coverage_files = glob.glob(os.path.join(tmpdir, ".coverage*"))
    for f in coverage_files:
        shutil.move(os.path.abspath(f),
                    os.path.join(current_dir, os.path.basename(f)))

    if not DEBUG:
        shutil.rmtree(tmpdir)
    ok_(not fail, msg)


def get_tests_directory():
    """discover and return the absolute path of the root directory"""
    testdir = os.getcwd()
    if not testdir.endswith("tests"):
        testdir = os.path.join(testdir, "tests")
        if not os.path.exists(testdir):
            raise ValueError("can not find test directory")

    return testdir


def test_tool():
    '''yield list of scripts to test.'''
    # the current directory
    current_dir = os.getcwd()

    # directory location of tests
    testing_dir = get_tests_directory()

    # directory location of scripts
    tool_dir = os.path.join(os.path.dirname(testing_dir), "umi_tools")

    #for test_script in test_dirs:

    check_main.description = os.path.join(tool_dir, "def_main")

    tool = "tests/umi_tools.py"
    tool_name = os.path.basename(tool)

    yield (check_main,
           os.path.abspath(os.path.join(tool_dir, tool_name)))

    fn = 'tests/tests.yaml'
    assert os.path.exists(fn), "tests.yaml does not exist!"

    tool_tests = yaml.safe_load(open(fn))

    for test, values in sorted(list(tool_tests.items())):
        check_script.description = os.path.join(tool_name, test)
        if "skip_python" in values:
            versions = [x.strip() for x in
                        str(values["skip_python"]).split(",")]
            versions = [x for x in versions
                        if PYTHON_VERSION.startswith(x)]
            if len(versions) > 0:
                continue

        yield(check_script,
              test,
              values.get('stdin', None),
              values['options'],
              values['outputs'],
              values['references'],
              current_dir,
              values.get('sort', False))


def _read(fn):
    if fn.endswith(".gz"):
        with gzip.open(fn) as inf:
            data = inf.read()
    else:
        with open(fn, "rb") as inf:
            data = inf.read()

    if IS_PY3:
        try:
            data = data.decode("ascii")
        except UnicodeDecodeError:
            return data

    data = [x for x in data.splitlines()
            if not x.startswith("#")]

    return data
