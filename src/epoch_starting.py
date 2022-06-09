# -*- coding: utf-8 -*-
import os, shutil
import numpy as np
from multiprocessing import Pool

from . import utils


def init_rep(i,cfg,pool,d="epoch01"):
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

    kwargs = {"c": "start.gro", "f": "../../"+cfg.mdp,
                                "n": "../../"+cfg.ndx, "p": "../../"+cfg.topol,
                                "o": "mdrun.tpr", "maxwarn": str(cfg.maxwarn),
                                "directory": "%s/rep%02d"%(d, i)
                                }

    if(cfg.restraint_file=="initial"):
        kwargs["r"] = "../../initial/start.gro"
    elif(cfg.restraint_file=="start"):
        kwargs["r"] = "start.gro"
    elif(cfg.restraint_file):
        kwargs["r"] = cfg.restraint_file




    # Do the gromacs job asynchronously in the worker pool
    rc=pool.apply_async(utils.gromacs_command,
                        args=(cfg.gmx, "grompp"),
                        kwds=kwargs
                        )


    return ("%s/rep%02d"%(d,i),rc)


def next_rep(i,cfg,newepoch,oldepoch,rep,frm, val, pool):
    """ Initializes rep i of newepoch, taking the frame frm from rep of oldepoch
    """
    os.makedirs("epoch%02d/rep%02d"%(newepoch, i))

    print("Reading frame %d of epoch %d, rep %d, fval %f"%(frm, oldepoch, rep, val))
    cfg.struct.load_new("epoch%02d/rep%02d/mdrun.xtc"%(oldepoch, rep))
    cfg.struct.trajectory[frm]

    cfg.struct.atoms.write("epoch%02d/rep%02d/start.gro"%(newepoch, i))
    print("Wrote structure to epoch%02d/rep%02d/start.gro, starting to grompp..."%(newepoch, i))

    # Make a note of the origin
    with open("epoch%02d/rep%02d/origin.txt"%(newepoch, i), "w") as f:
        f.write("# Origin of this trajectory:\n")
        f.write("# epoch rep frame val\n")
        f.write("%d %d %d %f\n"%(oldepoch,rep,frm, val))

    # Do the gromacs job asynchronously in the worker pool
    rc=pool.apply_async(utils.gromacs_command,
                        args=(cfg.gmx, "grompp"),
                        kwds={"c": "start.gro", "f": "../../"+cfg.mdp,
                                "n": "../../"+cfg.ndx, "p": "../../"+cfg.topol,
                                "o": "mdrun.tpr", "maxwarn": str(cfg.maxwarn+cfg.maxwarn_add),
                                "directory": "epoch%02d/rep%02d"%(newepoch, i)
                                }
                        )


    return ("epoch%02d/rep%02d"%(newepoch,i),rc)


def start_epoch(nextepoch, cfg, val=None, epc=None, rep=None, frm=None):
    """ Start epoch nextepoch. If it is one, the initial epoch is started,
          otherwise the following parameters should be supplied. cfg is the config.
       "Optional" input parameters:
            - val : Length N array of the fvals for each new rep
            - epc : Length N array of the epochs each new rep comes from
            - rep : Length N array of the rep within the epoch each new rep comes from
            - frm : Length N array of the frame within the rep each new rep comes from
    """

    # The grompping itself will be done asynchronously in the worker pool
    # No point in having more than 4 workers, since reading the frame from file is done
    # in the main process. Only exception is initialization, as there is no need to read a trajectory. 
    numproc = min(os.cpu_count(),4,cfg.N) if(nextepoch!=1) else min(os.cpu_count(),cfg.N)

    with Pool(numproc) as p:
        res = []
        if(nextepoch==1):
            # Initial structures and first epoch
            os.makedirs("epoch01", exist_ok=True)

            for i in range(cfg.N):
                dostop=False
                for d,rc in res:
                    if(rc.ready() and rc.get()!=0):
                        dostop=True
                        break

                if(dostop): break                        
                res.append(init_rep(cfg.first_rep_num + i,cfg,p))
        elif(np.any([d is None for d in (val,epc,rep,frm)])):
            raise ValueError("val, epc, rep or frm cannot be None if nextepoch!=1")
        else:

            os.makedirs("epoch%02d"%(nextepoch),exist_ok=True)

            for i, (v, e, r, f) in enumerate(zip(val,epc,rep,frm)):
                dostop=False
                for d,rc in res:
                    if(rc.ready() and rc.get()!=0):
                        dostop=True
                        break

                if(dostop): break
                res.append(next_rep(cfg.first_rep_num + i,cfg, nextepoch, e, r, f, v, p))
        
        print("Waiting for grompping to finish")
        p.close()
        p.join()

    for d,rc in res:
        assert rc.get() == 0, "Nonzero returncode from grompp, see %s/output_grompp.txt for more detail."%(d)

    # copy sbatch template
    with open("templates/sbatch_launch.sh") as fin:
        with open("epoch%02d/sbatch_launch.sh"%nextepoch, "w") as fout:
            fout.write(fin.read().format(i=nextepoch, account=cfg.account, email=cfg.email))
