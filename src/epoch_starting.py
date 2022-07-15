# -*- coding: utf-8 -*-
import os, time
import numpy as np
import warnings
from multiprocessing import Pool
import pathlib

from . import utils
from . import inout

# Type hints
from typing import Any, Tuple, Optional
from MDAnalysis.core.groups import AtomGroup
from multiprocessing.pool import AsyncResult
from numpy.typing import NDArray


def init_rep(i: int,cfg: Any, atoms: AtomGroup, pool: Pool, d: str="epoch01") -> Tuple[pathlib.Path,AsyncResult]:
    """ Initializes rep i from atom group atoms
    """
    d = pathlib.Path(d) / ("rep%02d"%i)
    d.mkdir()

    # If not found, just load the default
    atoms.write(str(d / "start.gro"))

    print("Wrote start.gro to rep%02d, starting to grompp..."%(i))

    kwargs = {"c": pathlib.Path("start.gro"),
              "f": pathlib.Path("..", "..") / cfg.mdp,
              "n": pathlib.Path("..", "..") / cfg.ndx,
              "p": pathlib.Path("..", "..") / cfg.topol,
              "o": pathlib.Path("mdrun.tpr"),
              "directory": d }

    if(cfg.maxwarn):
        kwargs["maxwarn"] = str(cfg.maxwarn)
    if(cfg.restraint_file=="initial"):
        kwargs["r"] = pathlib.Path("..", "..") / cfg.initial_struct
    elif(cfg.restraint_file=="start"):
        kwargs["r"] = "start.gro"
    elif(cfg.restraint_file):
        kwargs["r"] = cfg.restraint_file




    # Do the gromacs job asynchronously in the worker pool
    rc=pool.apply_async(utils.gromacs_command,
                        args=(cfg.gmx, "grompp"),
                        kwds=kwargs
                        )


    return (d,rc)


def next_rep(i: int, cfg: Any, newepoch: int, oldepoch: int, rep: int, frm: int, val: float, pool: Pool) -> Tuple[pathlib.Path,AsyncResult]:
    """ Initializes rep i of newepoch, taking the frame frm from rep of oldepoch
    """
    d = pathlib.Path("epoch%02d"%newepoch) / ("rep%02d"%i)
    d.mkdir()

    print("Reading frame %d of epoch %d, rep %d, fval %f"%(frm, oldepoch, rep, val))
    d_old = pathlib.Path("epoch%02d"%oldepoch) / ("rep%02d"%rep)
    next_xtc = d_old / "mdrun.xtc"
    if(next_xtc != cfg.struct.filename):
        cfg.struct.load_new(str(next_xtc))
    cfg.struct.trajectory[frm]

    cfg.struct.atoms.write(str(d / "start.gro"))
    print("Wrote structure to %s, starting to grompp..."%(d / "start.gro"))

    # Make a note of the origin
    with (d / "origin.txt").open("w") as f:
        f.write("# Origin of this trajectory:\n")
        f.write("# epoch rep frame val\n")
        f.write("%d %d %d %f\n"%(oldepoch,rep,frm, val))

    kwargs = {"c": pathlib.Path("start.gro"),
              "f": pathlib.Path("..", "..") / cfg.mdp,
              "n": pathlib.Path("..", "..") / cfg.ndx,
              "p": pathlib.Path("..", "..") / cfg.topol,
              "o": pathlib.Path("mdrun.tpr"),
              "directory": d }

    if(cfg.maxwarn):
        kwargs["maxwarn"] = cfg.maxwarn
    if(cfg.restraint_file=="initial"):
        kwargs["r"] = pathlib.Path("..", "..") / cfg.initial_struct
    elif(cfg.restraint_file=="start"):
        kwargs["r"] = pathlib.Path("start.gro")
    elif(cfg.restraint_file):
        kwargs["r"] = pathlib.Path("..", "..") / cfg.restraint_file

    # Do the gromacs job asynchronously in the worker pool
    rc=pool.apply_async(utils.gromacs_command,
                        args=(cfg.gmx, "grompp"),
                        kwds=kwargs
                        )


    return (d,rc)


def start_epoch(nextepoch: int, cfg: Any,
                val: Optional[NDArray[np.float_]] = None,
                epc: Optional[NDArray[np.int_]]   = None,
                rep: Optional[NDArray[np.int_]]   = None,
                frm: Optional[NDArray[np.int_]]   = None) -> None:
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
            struct = inout.load_starter_structures(cfg.initial_struct)
            num_frames = len(struct.trajectory)
            print("Found %d starting structures"%num_frames)
            if(num_frames>cfg.N):
                warnings.warn(f"{num_frames} starting structures found, but only {cfg.N} " \
                              f"repetitions will be started. Only the first {cfg.N} struc" \
                                "tures will be used. If this was intentional, feel free t" \
                                "o ignore this warning.\nContinuing in 2 seconds.",
                                UserWarning)
                time.sleep(2)
                
            os.makedirs("epoch01", exist_ok=True)

            for i in range(cfg.N):
                dostop=False
                for d,rc in res:
                    if(rc.ready() and rc.get()!=0):
                        dostop=True
                        break

                if(dostop): break

                struct.trajectory[i%num_frames]
                res.append(init_rep(cfg.first_rep_num + i,cfg,struct.atoms,p))
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
    utils.copy_sbatch_template(cfg.sbatch, pathlib.Path("epoch%02d"%nextepoch) / cfg.sbatch.name, nextepoch, cfg)
