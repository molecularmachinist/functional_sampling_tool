# -*- coding: utf-8 -*-
import os, shutil

from . import utils


def init_rep(i,cfg,d="epoch01"):
    """ Initializes rep i
    """
    try:
        os.makedirs("%s/rep%02d"%(d,i))
    except FileExistsError:
        print("Directory %s/rep%02d exists, skipping rep %d"%(d,i,i))
        return

    # If not found, just load the default
    shutil.copyfile("initial/start.gro", "%s/rep%02d/start.gro"%(d,i))
    print("Copied initial/start.gro to rep%02d, starting to grompp..."%(i))

    # Save original working dir to come back to
    prevdir = os.getcwd()
    try:
        os.chdir("%s/rep%02d"%(d,i))

        rc=utils.gromacs_command(cfg.gmx, "grompp", c="start.gro", f="../../"+cfg.mdp, n="../../"+cfg.ndx,
                           p="../../"+cfg.topol, o="mdrun.tpr", maxwarn="1")

        print("Process returned %d"%rc)

    finally:
        # Whatever happens, we go back to the original working dir
        os.chdir(prevdir)



def next_rep(i,cfg,newepoch,oldepoch,rep,frm, val):
    """ Initializes rep i of newepoch, taking the frame frm from rep of oldepoch
    """
    os.makedirs("epoch%02d/rep%02d"%(newepoch, i))

    print("Reading frame %d of epoch %d, rep %d"%(frm, oldepoch, rep))
    init_struct = mdtraj.load_frame("epoch%02d/rep%02d/mdrun.xtc"%(oldepoch, rep), frm, top=cfg.struct)

    init_struct.save_gro("epoch%02d/rep%02d/start.gro"%(newepoch, i))
    print("Wrote structure to epoch%02d/rep%02d/start.gro, starting to grompp..."%(newepoch, i))

    # Make a note of the origin
    with open("epoch%02d/rep%02d/origin.txt"%(newepoch, i), "w") as f:
        f.write("# Origin of this trajectory:\n")
        f.write("# epoch rep frame val\n")
        f.write("%d %d %d %f\n"%(oldepoch,rep,frm, val))

    # Save original working dir to come back to
    prevdir = os.getcwd()
    try:
        os.chdir("epoch%02d/rep%02d"%(newepoch, i))

        # The pdb structure seems to change the atom names, so maxwarn is 2
        rc=utils.gromacs_command("grompp", c="start.gro", f="../../"+cfg.mdp, n="../../"+cfg.ndx,
                           p="../../"+cfg.topol, o="mdrun.tpr", maxwarn="2")

        print("Process returned %d"%rc)

    finally:
        # Whatever happens, we go back to the original working dir
        os.chdir(prevdir)

def start_epoch(nextepoch, cfg, val=None, epc=None, rep=None, frm=None):
    """ Start epoch nextepoch. If it is one, the initial epoch is started,
          otherwise the following parameters should be supplied. cfg is the config.
       "Optional" input parameters:
            - val : Length N array of the fvals for each new rep
            - epc : Length N array of the epochs each new rep comes from
            - rep : Length N array of the rep within the epoch each new rep comes from
            - frm : Length N array of the frame within the rep each new rep comes from
    """
    if(nextepoch==1):
        # Initial structures and first epoch
        os.makedirs("epoch01", exist_ok=True)

        for i in range(1,cfg.N+1):
            init_rep(i,cfg)
    elif(None in [val,epc,rep,frm]):
        raise ValueError("val, epc, rep or frm cannot be None if nextepoch!=1")
    else:

        os.makedirs("epoch%02d"%(nextepoch))

        for i, (v, e, r, f) in enumerate(zip(val,epc,rep,frm)):
            next_rep(i+1,cfg, nextepoch, epoch, rep, frm, val)

    # copy sbatch template
    with open("templates/sbatch_launch.sh") as fin:
        with open("epoch%02d/sbatch_launch.sh"%nextepoch, "w") as fout:
            fout.write(fin.read().format(i=nextepoch, account=cfg.account, email=cfg.email))
