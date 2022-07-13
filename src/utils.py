# -*- coding: utf-8 -*-
import subprocess as subp
import MDAnalysis as mda
import numpy as np
from typing import Any, Callable, TypeAlias,Union, Optional
import math, os
import inspect, hashlib
import warnings

# Types for typing
ts_type:TypeAlias = mda.coordinates.base.Timestep
transform_type:TypeAlias = Callable[[ts_type],ts_type]
ag_type:TypeAlias = mda.core.groups.AtomGroup
a_type:TypeAlias = mda.core.groups.Atom


class DeprecatedUsageWarning(UserWarning):
    pass


def rolling_mean(data: np.ndarray, window:int=10, center:bool = True, fill:float=np.nan) -> np.ndarray:
    if(center):
        start_offset = math.floor(window/2)
        end_offset   = -math.ceil(window/2)+1
        if(end_offset==0): end_offset=None
    else:
        start_offset = window-1
        end_offset   = None

    window-=1
    if(window==0): window=None

    cumsum = np.nancumsum(data, axis=0, dtype=float)
    mean   = np.full_like(cumsum, fill)
    mean[start_offset:end_offset] = (cumsum[window:]-cumsum[:-window])/window

    return mean


def read_ndx(ndx:str) -> dict[str,list[int]]:
    # Return empty dictionary if ndx is None
    if(ndx is None):
        return {}
    print("\nReading index groups from %s"%ndx)
    indexes = {}
    groups  = []
    current=None

    with open(ndx) as f:
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

    print("Found %d groups"%len(groups))
    print("    %-20s%10s\n"%("Group name", "atoms"))
    for i, g in enumerate(groups):
        print("%2d. %-20s%10d"%(i+1, g, len(indexes[g])))

    print()
    return indexes


def gromacs_command(gmx:str, cmd:str, *args:str, directory:str=".", input:Optional[bytes]=None, **kwargs:str):
    """ Call the gromacs subcommand cmd in directory. Both args and keys of kwargs should be without the leading dash.
        output is redirected to output_<cmd>.txt. Returns the return code of the command.
    """

    # Save original working dir to come back to
    prevdir = os.getcwd()
    try:
        os.chdir(directory)
        command = [gmx, cmd]+["-"+a for a in args]
        for k in kwargs:
            command += ["-"+k, kwargs[k]]

        print(f"Running: {' '.join(command)}")
        with open("output_%s.txt"%cmd, "w") as fout:
            compProc = subp.run(command, stdout=fout, stderr=subp.STDOUT, input=input)

    finally:
        # Whatever happens, we go back to the original working dir
        os.chdir(prevdir)
        
    return compProc.returncode



def rsync_command(send_from:str, send_to:str, excludes:list[str]=[]):
    """ A utility function for keeping the remote dir in sync.
        Basically run rsync with given source and target, and given excludes.
        Uses -ravP by default.
        Prints stdout and sterr when program finishes and returns the return code.
    """
    command = ["rsync", "-avP", "--partial"]+["--exclude=%s"%e for e in excludes] + [send_from, send_to]
    print(f"Running: {' '.join(command)}")
    compProc = subp.run(command)
    return compProc.returncode



def rsync_down(cfg:Any) -> None:
    """ Wrapper function for rsync_command, to sync the local dir to the remote,
        ie. "pull down"
    """
    rc = rsync_command("%s:%s/"%(cfg.remote_name,cfg.remote_dir), "./", excludes=cfg.rsync_excludes)
    print("Process returned %d"%rc)


def rsync_up(cfg:Any) -> None:
    """ Wrapper function for rsync_command, to sync the remote dir to the local,
        ie. "push up"
    """
    rc = rsync_command( "./", "%s:%s/"%(cfg.remote_name,cfg.remote_dir), excludes=cfg.rsync_excludes)
    print("Process returned %d"%rc)

def make_pdb(cfg:Any) -> None:
    """
    Makes the pdb to load mdtraj
    """
    print("Making initial/start.pdb")

    # Save original working dir to come back to
    prevdir = os.getcwd()
    try:
        os.chdir("initial")

        rc=gromacs_command(cfg.gmx, "trjconv", f="start.gro",
        s="../epoch01/rep01/mdrun.tpr", o="start.pdb", input=b"System")

        print("Process returned %d"%rc)
        assert rc == 0, "Nonzero returncode from trjconv, see initial/output_trjconv.txt for more detail."
    finally:
        # Whatever happens, we go back to the original working dir
        os.chdir(prevdir)

def __depr_copy_sbatch_template(fin:str,fout:str,enum:int,cfg:Any) -> None:
    """
    To be deprecated. Gets email and account from config.
    """
    warnings.warn("email and account for 'sbatch_launch.sh' should no longer come from config. " +
                    "It will now be filled in as before, but this will be " +
                    "deprecated in the future.", DeprecatedUsageWarning)
    with open(fin) as f:
        with open(fout, "w") as fo:
            fo.write(f.read().format(i=enum, account=cfg.account, email=cfg.email))

def copy_sbatch_template(fin:str,fout:str,enum:int,cfg:Any=None) -> None:
    """
    Copies the sbatch template from fin to fout.

    To be depracated: if the template has {account} {email} and {i}, it will be filled as before, from cfg
    """
    use_depr = False
    with open(fin) as fi:
        templ = fi.read()
        if("{i}" in templ and "{account}" in templ and "{email}" in templ):
            use_depr = True
    if(use_depr):
        __depr_copy_sbatch_template(fin,fout,enum,cfg)
        return
    with open(fin) as fi:
        with open(fout,"w") as fo:
            for line in fi:
                if(line.startswith("##")):
                    continue

                fo.write(line.replace("{epoch_num}", str(enum)))

def hash_func(f:Callable) -> str:
    """
    Reads function as string and calulates the md5sum as a hex string
    """
    fstr = inspect.getsource(f)
    return hashlib.md5(fstr.encode('utf-8')).hexdigest()




def load_sel(sel_str:str, struct:Union[mda.Universe,ag_type], ndx:dict[str,list[int]]) -> ag_type:
    if(sel_str in ndx):
        sel = struct.atoms[np.array(ndx[sel_str])-1]
    else:
        sel = struct.select_atoms(sel_str)
    return sel