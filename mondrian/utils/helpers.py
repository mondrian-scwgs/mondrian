'''
Created on Feb 19, 2018

@author: dgrewal
'''
import errno
import gzip
import logging
import os

import pandas as pd
from subprocess import Popen, PIPE


def run_cmd(cmd, output=None):
    stdout = PIPE
    if output:
        stdout = open(output, "w")

    p = Popen(cmd, stdout=stdout, stderr=PIPE)

    cmdout, cmderr = p.communicate()
    retc = p.returncode

    if retc:
        raise Exception(
            "command failed. stderr:{}, stdout:{}".format(
                cmdout,
                cmderr))

    if output:
        stdout.close()

    print(cmdout)
    print(cmderr)


class getFileHandle(object):
    def __init__(self, filename, mode='rt'):
        self.filename = filename
        self.mode = mode

    def __enter__(self):
        if self.get_file_format(self.filename) in ["csv", 'plain-text']:
            self.handle = open(self.filename, self.mode)
        elif self.get_file_format(self.filename) == "gzip":
            self.handle = gzip.open(self.filename, self.mode)
        elif self.get_file_format(self.filename) == "h5":
            self.handle = pd.HDFStore(self.filename, self.mode)
        return self.handle

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.handle.close()

    def get_file_format(self, filepath):
        if filepath.endswith('.tmp'):
            filepath = filepath[:-4]

        _, ext = os.path.splitext(filepath)

        if ext == ".csv":
            return "csv"
        elif ext == ".gz":
            return "gzip"
        elif ext == ".h5" or ext == ".hdf5":
            return "h5"
        elif ext == '.yaml':
            return 'plain-text'
        else:
            logging.getLogger("single_cell.helpers").warning(
                "Couldn't detect output format. extension {}".format(ext)
            )
            return "plain-text"


def makedirs(directory, isfile=False):
    if isfile:
        directory = os.path.dirname(directory)
        if not directory:
            return

    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
