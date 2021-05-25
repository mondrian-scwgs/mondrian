'''
Created on Feb 19, 2018

@author: dgrewal
'''
import errno
import gzip
import logging
import os
import shutil
import tarfile

import pandas as pd
import yaml
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



class InputException(Exception):
    pass

def copyfile(source, dest):
    shutil.copyfile(source, dest)


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


def get_compression_type_pandas(filepath):
    if get_file_format(filepath) == 'gzip':
        return 'gzip'
    else:
        return None


def get_file_format(filepath):
    if filepath.endswith('.tmp'):
        filepath = filepath[:-4]

    _, ext = os.path.splitext(filepath)

    if ext == ".csv":
        return "csv"
    elif ext == ".gz":
        return "gzip"
    elif ext == ".h5" or ext == ".hdf5":
        return "h5"
    else:
        logging.getLogger("single_cell.plot_metrics").warning(
            "Couldn't detect output format. extension {}".format(ext)
        )
        return "csv"


def write_to_yaml(outfile, data):
    with open(outfile, 'w') as output:
        yaml.safe_dump(data, output, default_flow_style=False)


def eval_expr(val, operation, threshold):
    if operation == "gt":
        if val > threshold:
            return True
    elif operation == 'ge':
        if val >= threshold:
            return True
    elif operation == 'lt':
        if val < threshold:
            return True
    elif operation == 'le':
        if val <= threshold:
            return True
    elif operation == 'eq':
        if val == threshold:
            return True
    elif operation == 'ne':
        if not val == threshold:
            return True
    elif operation == 'in':
        if val in threshold:
            return True
    elif operation == 'notin':
        if not val in threshold:
            return True
    else:
        raise Exception("unknown operator type: {}".format(operation))

    return False


def filter_metrics(metrics, cell_filters):
    # cells to keep
    for metric_col, operation, threshold in cell_filters:
        if metrics.empty:
            continue

        rows_to_keep = metrics[metric_col].apply(eval_expr, args=(operation, threshold))

        metrics = metrics[rows_to_keep]

    return metrics


def get_incrementing_filename(path):
    """
    avoid overwriting files. if path exists then return path
    otherwise generate a path that doesnt exist.
    """

    if not os.path.exists(path):
        return path

    i = 0
    while os.path.exists("{}.{}".format(path, i)):
        i += 1

    return "{}.{}".format(path, i)


def build_shell_script(command, tag, tempdir):
    outfile = os.path.join(tempdir, "{}.sh".format(tag))
    with open(outfile, 'w') as scriptfile:
        scriptfile.write("#!/bin/bash\n")
        if isinstance(command, list) or isinstance(command, tuple):
            command = ' '.join(command) + '\n'
        scriptfile.write(command)
    return outfile



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


def make_tarfile(output_filename, source_dir):
    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))


def extract_tar(input_tar, outdir):
    with tarfile.open(input_tar) as tar:
        tar.extractall(path=outdir)


def gunzip_file(infile, outfile):
    with gzip.open(infile, 'rb') as f_in:
        with open(outfile, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
