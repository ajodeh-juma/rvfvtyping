#!/usr/bin/env python

"""
utility functions
"""

import os
import errno
import shutil
import argparse
import subprocess
from datetime import datetime
from itertools import groupby

# time format
f = '%Y-%m-%d %H:%M:%S'


def start_time():
    """
    determine the start time of a process
    :return:
    """
    # set start time format
    st = datetime.now()
    s = st.strftime(f)
    print("Started:\t{}".format(s))
    return st


def end_time(stime):
    """
    :param stime:
    :return:
    """
    # set end time format
    et = datetime.now()
    e = et.strftime(f)
    print("Completed:\t{}".format(e))
    return et


def run_time(stime, etime):
    """
    determine the total time to run a process
    :param stime:
    :param etime:
    :return:
    """

    tdelta = etime - stime

    # format the time delta object to human readable form
    d = dict(days=tdelta.days)
    d['hrs'], rem = divmod(tdelta.seconds, 3600)
    d['min'], d['sec'] = divmod(rem, 60)

    if d['min'] == 0:
        fmt = '{sec} sec'
    elif d['hrs'] == 0:
        fmt = '{min} min {sec} sec'
    elif d['days'] == 0:
        fmt = '{hrs} hr(s) {min} min {sec} sec'
    else:
        fmt = '{days} day(s) {hrs} hr(s) {min} min {sec} sec'
    print("[ALL done] Runtime: " + '\t' + fmt.format(**d))


def find_executable(names, default=None):
    """
    find an executable PATH from the given list of names.
    Raises an error if no executable is found.  Provide a
    value for the "default" parameter to instead return a value.

    :param names: list of given executable names
    :param default:
    :return: <str> path to the first executable found in PATH from the given list of names.
    """
    exe = next(filter(shutil.which, names), default)

    if exe is None:
        print("Unable to find any of {} in PATH={}".format(names, os.environ["PATH"]))
        print("\nHint: You can install the missing program using conda or homebrew or apt-get.\n")
        raise Exception
    return exe


def run_shell_command(cmd, raise_errors=False, extra_env=None):
    """
    run the given command string via Bash with error checking

    :param cmd: command given to the bash shell for executing
    :param raise_errors: bool to raise error if running command fails/succeeds
    :param extra_env: mapping that provides keys and values which are overlayed onto the default subprocess environment.
    :return:
    """

    env = os.environ.copy()

    if extra_env:
        env.update(extra_env)

    try:
        p = subprocess.check_call("set -euo pipefail; " + cmd,
                                  shell=True,
                                  universal_newlines=True,
                                  executable="/bin/bash"
                                  )
    except (subprocess.CalledProcessError, OSError) as error:
        rc = error.returncode
        if rc == 127:
            extra = "Are you sure this program is installed?"
        else:
            extra = " "
        print("Error occurred: shell exited with return code: {}\ncommand running: {}".format(
            error.returncode, cmd, extra)
        )
        if raise_errors:
            raise
        else:
            return False
    except FileNotFoundError as error:
        print("Unable to run shell command using {}! tool requires {} to be installed.".format(
            error.filename, error.filename)
        )
        if raise_errors:
            raise
        else:
            return False
    else:
        return True


def mkdir(directory):
    """
    recursivley create a directory if it does not exist

    :param directory: path to the directory to be created
    :return: directory
    """
    directory = os.path.abspath(directory)
    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise  # raises the error again
    return directory


def fasta_iterator(fasta_name):
    """
    Given a fasta file. yield tuples of header, sequence
    Author: Brent Pedersen
    https://www.biostars.org/p/710/
    """
    fa_fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fa_fh, lambda line: line[0] == ">"))

    for header in faiter:
        # drop the ">"
        header = next(header)[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq
    fa_fh.close()


def str2bool(v):
    """
    convert string to a boolean

    :param v: string value either true or false
    :return: False or True
    """
    if isinstance(v, bool):
        return v
    if v.lower() in ['true', 't']:
        return True
    elif v.lower() in ['false', 'f']:
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')
