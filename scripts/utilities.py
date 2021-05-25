import pickle
import json
import csv
import re
import tarfile

def load(filename):
    """Load contents of parsed results file."""
    D = pickle.load(open(filename, "rb"))
    R = D['setup']['radius']
    H = D['setup']['height']
    T = D['setup']['time']
    N = D['agents'].shape[0]
    C = D['setup']['coords']
    POPS = D['setup']['pops']
    TYPES = D['setup']['types']
    return D, R, H, T, N, C, POPS, TYPES

def load_json(json_file, tar=None):
    """Load .json file."""
    if tar:
        file = tar.extractfile(json_file)
        contents = [line.decode("utf-8") for line in file.readlines()]
        return json.loads("".join(contents))
    else:
        return json.load(open(json_file, "r"))

def is_tar(file):
    """Check if file has .tar.xz extension."""
    return file[-7:] == ".tar.xz"
