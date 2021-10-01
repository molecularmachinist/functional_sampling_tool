# -*- coding: utf-8 -*-
import os, shutil




def init_rep(i,d="epoch01"):
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

        rc=gromacs_command("grompp", c="start.gro", f="../../"+mdp, n="../../"+ndx,
                           p="../../"+topol, o="mdrun.tpr", maxwarn="1")

        print("Process returned %d"%rc)

    finally:
        # Whatever happens, we go back to the original working dir
        os.chdir(prevdir)


def start_epoch(nextepoch, cfg, chosen_bins=None):
    """ Start epoch nextepoch. If it is one, the initial epoch is started,
        otherwise chosen_bins should be supplied.
    """
    if(nextepoch==1):
        # Initial structures and first epoch
        os.makedirs("epoch01", exist_ok=True)

        for i in range(1,N+1):
            init_rep(i)
    elif(chosen_bins is None):
        raise ValueError("chosen_bins cannot be None if nextepoch!=1")
    else:
        v, e, r, f = choose_frames(chosen_bins)
        print(v, e, r, f)

        os.makedirs("epoch%02d"%(nextepoch))

        for i, (val, epoch, rep, frm) in enumerate(zip(v,e,r,f)):
            next_rep(i+1, nextepoch, epoch, rep, frm, val)

    # copy sbatch template
    with open("templates/sbatch_launch.sh") as fin:
        with open("epoch%02d/sbatch_launch.sh"%nextepoch, "w") as fout:
            fout.write(fin.read().format(i=nextepoch, account=cfg.account, email=cfg.email))
