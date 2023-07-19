import numpy as np
import pathlib
import os
import sys

datadir = pathlib.Path(__file__).parent
epoch_example_dir = datadir / "epoch_examples"


def run_in_dir(dir):
    def decorator(func):
        def function_in_dir(*args, **kwargs):
            cwd = os.getcwd()
            try:
                os.chdir(dir)
                res = func(*args, **kwargs)
            finally:
                os.chdir(cwd)
            return res
        return function_in_dir
    return decorator


def silence_function(func):
    def silent_function(*args, **kwargs):
        with open(os.devnull, 'w') as devnull:
            try:
                sys.stdout = devnull
                res = func(*args, **kwargs)
            finally:
                sys.stdout = sys.__stdout__
        return res
    return silent_function


def unflatten_dict(dic):
    new_dict = {}
    for key in dic:
        parts = key.split("/")
        parent_dict = new_dict
        for p in parts[:-1]:
            if (p not in parent_dict):
                parent_dict[p] = {}
            parent_dict = parent_dict[p]

        parent_dict[parts[-1]] = dic[key]
    return new_dict


def flatten_dict(dic):
    newdict = {}
    for key in dic:
        if (type(dic[key]) == dict):
            flat_child = flatten_dict(dic[key])
            for key2 in flat_child:
                newdict["%s/%s" % (key, key2)] = flat_child[key2]
        else:
            newdict[key] = dic[key]

    return newdict


def _load_data_npz(filename):
    with np.load(datadir / filename) as npz:
        data = unflatten_dict(npz)
    return data


def make_whole_data():
    return _load_data_npz("make_whole.npz")


@run_in_dir(epoch_example_dir)
@silence_function
def load_cfg():
    from functional_sampling_tool.inout import load_options
    return load_options(epoch_example_dir / "config.py")
