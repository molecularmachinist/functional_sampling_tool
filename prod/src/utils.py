# -*- coding: utf-8 -*-
import subprocess as subp

from . import config as cfg

def gromacs_command(cmd, *args, input=None, **kwargs):
    """ Call the gromacs subcommand cmd. Both args and keys of kwargs should be without the leading dash.
        output is redirected to output_<cmd>.txt. Returns the return code of the command.
    """
    command = [gmx, cmd]+["-"+a for a in args]
    for k in kwargs:
        command += ["-"+k, kwargs[k]]

    with open("output_%s.txt"%cmd, "w") as fout:
        compProc = subp.run(command, stdout=fout, stderr=subp.STDOUT, input=input)

    return compProc.returncode



def rsync_command(send_from, send_to, excludes=[]):
    """ A utility function for keeping the remote dir in sync.
        Basically run rsync with given source and target, and given excludes.
        Uses -ravP by default.
        Prints stdout and sterr when program finishes and returns the return code.
    """
    command = ["rsync", "-ravP"]+["--exclude=%s"%e for e in excludes] + [send_from, send_to]
    compProc = subp.run(command,capture_output=True,text=True)
    print(compProc.stdout)
    print(compProc.stderr)
    return compProc.returncode



def rsync_down():
    """ Wrapper function for rsync_command, to sync the local dir to the remote,
        ie. "pull down"
    """
    rc = rsync_command("%s:%s/"%(cfg.remote_name,cfg.remote_dir), "./", excludes=cfg.rsync_excludes)
    print("Process returned %d"%rc)


def rsync_up():
    """ Wrapper function for rsync_command, to sync the remote dir to the local,
        ie. "push up"
    """
    rc = rsync_command( "./", "%s:%s/"%(cfg.remote_name,cfg.remote_dir), excludes=cfg.rsync_excludes)
    print("Process returned %d"%rc)
