# -*- coding: utf-8 -*-
import subprocess as subp
import numpy as np
import math
import inspect, hashlib



def rolling_mean(data, window=10, center = True, fill=np.nan):
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

def gromacs_command(gmx, cmd, *args, input=None, **kwargs):
    """ Call the gromacs subcommand cmd. Both args and keys of kwargs should be without the leading dash.
        output is redirected to output_<cmd>.txt. Returns the return code of the command.
    """
    command = [gmx, cmd]+["-"+a for a in args]
    for k in kwargs:
        command += ["-"+k, kwargs[k]]

    print(f"Running: {' '.join(command)}")
    with open("output_%s.txt"%cmd, "w") as fout:
        compProc = subp.run(command, stdout=fout, stderr=subp.STDOUT, input=input)

    return compProc.returncode



def rsync_command(send_from, send_to, excludes=[]):
    """ A utility function for keeping the remote dir in sync.
        Basically run rsync with given source and target, and given excludes.
        Uses -ravP by default.
        Prints stdout and sterr when program finishes and returns the return code.
    """
    command = ["rsync", "-avP", "--partial"]+["--exclude=%s"%e for e in excludes] + [send_from, send_to]
    print(f"Running: {' '.join(command)}")
    compProc = subp.run(command)
    return compProc.returncode



def rsync_down(cfg):
    """ Wrapper function for rsync_command, to sync the local dir to the remote,
        ie. "pull down"
    """
    rc = rsync_command("%s:%s/"%(cfg.remote_name,cfg.remote_dir), "./", excludes=cfg.rsync_excludes)
    print("Process returned %d"%rc)


def rsync_up(cfg):
    """ Wrapper function for rsync_command, to sync the remote dir to the local,
        ie. "push up"
    """
    rc = rsync_command( "./", "%s:%s/"%(cfg.remote_name,cfg.remote_dir), excludes=cfg.rsync_excludes)
    print("Process returned %d"%rc)

def make_pdb(cfg):
    """
    Makes the pdb to load mdtraj
    """
    print("Making initial/start.pdb")

    # Save original working dir to come back to
    prevdir = os.getcwd()
    try:
        os.chdir("initial")

        rc=gromacs_command(cfg.gmx, "trjconv", f="start.gro",
        s="../epoch01/rep01/mdrun.tpr", o="start.pdb")

        print("Process returned %d"%rc)
        assert rc == 0, "Nonzero returncode from trjconv, see initial/output_trjconv.txt for more detail."
    finally:
        # Whatever happens, we go back to the original working dir
        os.chdir(prevdir)


def hash_func(f):
    """
    Reads function as string and calulates the md5sum as a hex string
    """
    fstr = inspect.getsource(f)
    return hashlib.md5(fstr.encode('utf-8')).hexdigest()
