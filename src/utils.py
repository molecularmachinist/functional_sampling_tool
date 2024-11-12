# -*- coding: utf-8 -*-
from numpy.typing import ArrayLike, NDArray
from MDAnalysis.coordinates.base import Timestep
from MDAnalysis.core.groups import AtomGroup
import MDAnalysis as mda
from typing import Any, Callable, Union, Optional, List, Dict
from . import __version__ as fst_version
import pathlib
import subprocess as subp
import numpy as np
import math
import sys
import re
import inspect
import hashlib

from .exceptions import (NoSbatchLaunchError,
                         NonzeroReturnError,
                         ExternalProgramMissingError)


# Type hinting
# Type aliases
transform_type = Callable[[Timestep], Timestep]


def rolling_mean(data: ArrayLike, window: int = 10,
                 center: bool = True, fill: float = np.nan) -> NDArray[np.float_]:
    if (window < 1):
        raise ValueError(
            "rolling mean window smaller than 1 "
            f"is nonsensical (is {window})")
    if (window == 1):
        return data
    if (center):
        start_offset = math.floor(window/2)
        end_offset = -math.ceil(window/2)+1
        if (end_offset == 0):
            end_offset = None
    else:
        start_offset = window-1
        end_offset = None

    window -= 1

    cumsum = np.nancumsum(data, axis=0, dtype=float)
    mean = np.full_like(cumsum, fill)
    mean[start_offset:end_offset] = (cumsum[window:]-cumsum[:-window])/window

    return mean


def check_num(prefix: pathlib.Path) -> List[int]:
    """
    Checks filenames prefix01, prefix02, etc and returns a list of integers that were found.
    """
    prog = re.compile("%s(?P<num>[0-9]+)$" % prefix.name)
    nums = []
    for f in prefix.parent.iterdir():
        m = prog.match(f.name)
        if (m):
            nums.append(int(m.group("num")))

    nums.sort()
    return nums


def read_ndx(ndx: pathlib.Path) -> Dict[str, List[int]]:
    # Return empty dictionary if ndx is None
    if (ndx is None):
        return {}
    print("\nReading index groups from %s" % ndx)
    indexes = {}
    groups = []
    current = None

    with ndx.open() as f:
        for line in f:
            line = line.strip()
            if (line.startswith("[") and line.endswith("]")):
                current = line[1:-1].strip()
                indexes[current] = []
                groups.append(current)
                continue
            parts = line.split()
            for i in parts:
                indexes[current].append(int(i))

    print("Found %d groups" % len(groups))
    print("    %-20s%10s\n" % ("Group name", "atoms"))
    for i, g in enumerate(groups):
        print("%2d. %-20s%10d" % (i+1, g, len(indexes[g])))

    print()
    return indexes


def gromacs_command(gmx: str, cmd: str, *args: Any, directory: pathlib.Path = pathlib.Path("."),
                    input: Optional[bytes] = None, nobackup=False, **kwargs: Any) -> int:
    """ Call the gromacs subcommand cmd in directory. Both args and keys of kwargs should be without the leading dash.
        output is redirected to output_<cmd>.txt. Returns the return code of the command.
    """

    command = [gmx]
    if (nobackup):
        command.append("-nobackup")
    command.append(cmd)
    command += ["-"+str(a) for a in args]
    for k in kwargs:
        command += ["-"+k, str(kwargs[k])]

    print(f"Running: {' '.join(command)}")
    with (directory / ("output_%s.txt" % cmd)).open("w") as fout:
        try:
            compProc = subp.run(command, stdout=fout,
                                stderr=subp.STDOUT, input=input,
                                cwd=directory)
        except FileNotFoundError as e:
            if (e.filename == gmx):
                raise ExternalProgramMissingError(
                    f"GROMACS command at '{gmx}' not found.\n"
                    "Original message:\n"
                    f"{e.__class__.__name__}: {e}"
                )
            raise

    return compProc.returncode


def rsync_command(send_from: Union[str, list], send_to: str, excludes: List[str] = []) -> int:
    """ A utility function for keeping the remote dir in sync.
        Basically run rsync with given source and target, and given excludes.
        Uses -ravP by default.
        Prints stdout and sterr when program finishes and returns the return code.
    """
    if (type(send_from) == list):
        command = ["rsync", "-avP", "--partial"] +\
            ["--exclude=%s" % e for e in excludes] + \
            send_from + \
            [send_to]
    else:
        command = ["rsync", "-avP", "--partial"] + \
            ["--exclude=%s" % e for e in excludes] + \
            [send_from, send_to]

    print(f"Running: {' '.join(command)}")
    try:
        compProc = subp.run(command)
    except FileNotFoundError as e:
        raise ExternalProgramMissingError(
            f"rsync command not found.\n"
            "Original message:\n"
            f"{e.__class__.__name__}: {e}"
        )
    return compProc.returncode


def rsync_down(cfg: Any, all=False) -> None:
    """ Wrapper function for rsync_command, to sync the local dir to the remote,
        ie. "pull down"
    """
    if (all):
        epoch = "epoch*"
        local_target = "."
    else:
        epoch_nums = check_num(pathlib.Path("epoch"))
        epoch = "epoch%02d/" % epoch_nums[-1]
        local_target = epoch
    rc = rsync_command(f"{cfg.remote_name}:{cfg.remote_dir}/{epoch}",
                       local_target,
                       excludes=cfg.rsync_excludes)
    if (rc):
        raise NonzeroReturnError("rsync process returned %d" % rc, code=rc)


def rsync_up(cfg: Any, all=False) -> None:
    """ Wrapper function for rsync_command, to sync the remote dir to the local,
        ie. "push up"
    """
    epoch_nums = check_num(pathlib.Path("epoch"))
    if (all):
        send_dirs = ["epoch%02d" % e for e in epoch_nums]
    else:
        send_dirs = ["epoch%02d" % epoch_nums[-1]]
    rc = rsync_command(send_dirs,
                       f"{cfg.remote_name}:{cfg.remote_dir}",
                       excludes=cfg.rsync_excludes)
    if (rc):
        raise NonzeroReturnError("rsync process returned %d" % rc, code=rc)


def copy_sbatch_template(fin: pathlib.Path, fout: pathlib.Path, enum: int, cfg: Any = None) -> None:
    """
    Copies the sbatch template from fin to fout.

    To be depracated: if the template has {account} {email} and {i}, it will be filled as before, from cfg
    """
    if (not fin.exists()):
        raise NoSbatchLaunchError(
            f"No file found at {fin}, make sure the launch script exists and the path is correct.")
    with fin.open() as fi:
        with fout.open("w") as fo:
            for line in fi:
                if (line.startswith("##")):
                    continue

                fo.write(line.replace("{epoch_num}", str(enum)))


def copy_config(fin: pathlib.Path, fout: pathlib.Path, default_values: Dict[str, Any] = {}):
    """
    Copy the config file to the specified location.
    """
    with fout.open("w")as fo:
        fo.write(f"# Copy of {fout}\n")
        fo.write("# Made with functional_sampling_tool "
                 f"version {fst_version}\n")
        fo.write("# Command line call:\n")
        fo.write("# "+" ".join(sys.argv)+"\n\n")

        with fin.open() as fi:
            fo.write(fi.read())

        fo.write("\n#  Default values: \n")
        for key, val in default_values.items():
            fo.write(f"{key} = {val!r}\n")


def hash_func(f: Callable) -> str:
    """
    Reads function as string and calculates the md5sum as a hex string
    """
    fstr = inspect.getsource(f)
    return hashlib.md5(fstr.encode('utf-8')).hexdigest()


def load_sel(sel_str: str,
             struct: Union[mda.Universe, AtomGroup],
             ndx: Dict[str, List[int]]) -> AtomGroup:
    if (sel_str in ndx):
        sel = struct.atoms[np.array(ndx[sel_str])-1]
    else:
        sel = struct.select_atoms(sel_str)
    return sel


def load_origin_data(filename: pathlib.Path, e: int) -> Dict[str, Union[int, float]]:
    """
    Load the reps origin data from the given filename. If the file is not found, but the
    epoch is the first one, the info is guessed. If the epoch is not the first one a message
    is printed on stderr and an empty dict is returned.
    If the file is found, but parsing fails, a message is printed in stderr and empty dict returned.
    """
    data = {}
    try:
        with (filename).open() as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.split()
                break

        if (len(parts) != 4):
            raise ValueError()
        data = {"epc":  int(parts[0]),
                "rep":  int(parts[1]),
                "frm":  int(parts[2]),
                "fval": float(parts[3])}

    except FileNotFoundError:
        if (e == 1):
            data = {"epc":  0,
                    "rep":  0,
                    "frm":  0,
                    "fval": float("nan")}
        else:
            print(f"origin file in {filename} could not be found. "
                  "This might be a problem if you're making an ancestry graph",
                  file=sys.stderr)
    except ValueError:
        print(f"origin file in {filename} could not be parsed. "
              "This might be a problem if you're making an ancestry graph",
              file=sys.stderr)

    return data
