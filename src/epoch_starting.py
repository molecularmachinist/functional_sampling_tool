# -*- coding: utf-8 -*-
import os
import time
import numpy as np
import warnings
from multiprocessing import Pool
import pathlib
import MDAnalysis as mda

from .exceptions import RepExistsError, NonzeroReturnError

from . import utils

# Type hints
from typing import Any, Tuple, Optional
from MDAnalysis.core.groups import AtomGroup
from multiprocessing.pool import AsyncResult
from numpy.typing import NDArray


def init_rep(i: int,
             cfg: Any,
             atoms: AtomGroup,
             origin: Tuple[int, int],
             pool: Pool,
             d: str = "epoch01") -> Tuple[pathlib.Path, AsyncResult]:
    """ Initializes rep i from atom group atoms
    """
    d = pathlib.Path(d) / ("rep%02d" % i)
    if (d.exists()):
        raise RepExistsError(
            f"{d} already exists, will not try to overwrite. Try running \"fst clean\" to remove the latest epoch directory.")
    d.mkdir()

    atoms.write(d / "start.gro")

    # Make a note of the origin
    with (d / "origin.txt").open("w") as f:
        f.write("# Origin of this trajectory:\n")
        f.write(f"# Frame is from {cfg.initial_struct[origin[0]]}\n")
        f.write("# Initial value not yet calculated, set to NaN\n")
        f.write("# Rep is the index of the origin file in 'initial_struct'\n")
        f.write("# epoch rep frame val\n")
        f.write("%d %d %d nan\n" % (0, origin[0], origin[1]))

    print("Wrote start.gro to rep%02d, starting to grompp..." % (i))

    kwargs = {"c": pathlib.Path("start.gro"),
              "f": pathlib.Path("..", "..") / cfg.mdp,
              "n": pathlib.Path("..", "..") / cfg.ndx,
              "p": pathlib.Path("..", "..") / cfg.topol,
              "o": pathlib.Path("mdrun.tpr"),
              "directory": d}

    if (cfg.maxwarn):
        kwargs["maxwarn"] = str(cfg.maxwarn)
    if (cfg.restraint_file == "initial"):
        kwargs["r"] = pathlib.Path("..", "..") / cfg.initial_struct[0]
    elif (cfg.restraint_file == "start"):
        kwargs["r"] = "start.gro"
    elif (cfg.restraint_file):
        kwargs["r"] = cfg.restraint_file

    # Do the gromacs job asynchronously in the worker pool
    rc = pool.apply_async(utils.gromacs_command,
                          args=(cfg.gmx, "grompp"),
                          kwds=kwargs
                          )

    return (d, rc)


def next_rep(i: int,
             cfg: Any,
             newepoch: int,
             oldepoch: int,
             rep: int,
             frm: int,
             val: float,
             pool: Pool) -> Tuple[pathlib.Path, AsyncResult]:
    """ Initializes rep i of newepoch, taking the frame frm from rep of oldepoch
    """
    d = pathlib.Path("epoch%02d" % newepoch) / ("rep%02d" % i)
    if (d.exists()):
        raise RepExistsError(
            f"{d} already exists, will not try to overwrite. Try running \"fst clean\" to remove the latest epoch directory.")
    d.mkdir()

    print("Reading frame %d of epoch %d, rep %d, fval %f" %
          (frm, oldepoch, rep, val))
    d_old = pathlib.Path("epoch%02d" % oldepoch) / ("rep%02d" % rep)
    next_xtc = d_old / "mdrun.xtc"
    if (next_xtc != cfg.struct.filename):
        cfg.struct.load_new(next_xtc)
    cfg.struct.trajectory[frm]

    cfg.struct.atoms.write(d / "start.gro")
    print("Wrote structure to %s, starting to grompp..." % (d / "start.gro"))

    # Make a note of the origin
    with (d / "origin.txt").open("w") as f:
        f.write("# Origin of this trajectory:\n")
        f.write("# epoch rep frame val\n")
        f.write("%d %d %d %f\n" % (oldepoch, rep, frm, val))

    kwargs = {"c": pathlib.Path("start.gro"),
              "f": pathlib.Path("..", "..") / cfg.mdp,
              "n": pathlib.Path("..", "..") / cfg.ndx,
              "p": pathlib.Path("..", "..") / cfg.topol,
              "o": pathlib.Path("mdrun.tpr"),
              "directory": d}

    if (cfg.maxwarn):
        kwargs["maxwarn"] = cfg.maxwarn
    if (cfg.restraint_file == "initial"):
        kwargs["r"] = pathlib.Path("..", "..") / cfg.initial_struct[0]
    elif (cfg.restraint_file == "start"):
        kwargs["r"] = pathlib.Path("start.gro")
    elif (cfg.restraint_file):
        kwargs["r"] = pathlib.Path("..", "..") / cfg.restraint_file

    # Do the gromacs job asynchronously in the worker pool
    rc = pool.apply_async(utils.gromacs_command,
                          args=(cfg.gmx, "grompp"),
                          kwds=kwargs
                          )

    return (d, rc)


def start_epoch(nextepoch: int, cfg: Any,
                val: Optional[NDArray[np.float_]] = None,
                epc: Optional[NDArray[np.int_]] = None,
                rep: Optional[NDArray[np.int_]] = None,
                frm: Optional[NDArray[np.int_]] = None,
                numproc: Optional[int] = None) -> None:
    """ Start epoch nextepoch. If it is one, the initial epoch is started,
          otherwise the following parameters should be supplied. cfg is the config.
       "Optional" input parameters:
            - val : Length N array of the fvals for each new rep
            - epc : Length N array of the epochs each new rep comes from
            - rep : Length N array of the rep within the epoch each new rep comes from
            - frm : Length N array of the frame within the rep each new rep comes from
    """
    # Make epoch dir
    edir = pathlib.Path("epoch%02d" % nextepoch)
    edir.mkdir(exist_ok=True)

    # The grompping itself will be done asynchronously in the worker pool
    if (numproc is None):
        # Default max number of threads is half the cpu count. It is assumed that due to hyperthreading the number of logical
        # cores is twice the number of physical cores
        numproc = min(os.cpu_count()//2, cfg.N)
    else:
        # If numproc is specified, it is still limited by N
        numproc = min(numproc, cfg.N)

    with Pool(numproc) as p:
        res = []
        if (nextepoch == 1):
            # Initial structures and first epoch
            struct = mda.Universe(*cfg.initial_struct)
            num_frames = len(struct.trajectory)
            print("Found %d starting structures" % num_frames)
            if (num_frames > cfg.N):
                warnings.warn(f"{num_frames} starting structures found, but only {cfg.N} "
                              f"repetitions will be started. Only the first {cfg.N} struc"
                              "tures will be used. If this was intentional, feel free t"
                              "o ignore this warning.\nContinuing in 2 seconds.",
                              UserWarning)
                time.sleep(2)

            for i in range(cfg.N):
                dostop = False
                for d, rc in res:
                    if (rc.ready() and rc.get() != 0):
                        dostop = True
                        break

                if (dostop):
                    break

                struct.trajectory[i % num_frames]

                if (len(cfg.initial_struct) == 1):
                    origin = (0, 0)
                elif (len(cfg.initial_struct) == 2):
                    origin = (1, i % num_frames)
                else:
                    ofile, oframe = struct.trajectory._get_local_frame(
                        struct.trajectory.frame)
                    origin = (ofile+1, oframe)
                res.append(
                    init_rep(cfg.first_rep_num + i,
                             cfg,
                             struct.atoms,
                             origin,
                             p)
                )
        elif (np.any([d is None for d in (val, epc, rep, frm)])):
            raise ValueError(
                "val, epc, rep or frm cannot be None if nextepoch!=1")
        else:

            for i, (v, e, r, f) in enumerate(zip(val, epc, rep, frm)):
                dostop = False
                for d, rc in res:
                    if (rc.ready() and rc.get() != 0):
                        dostop = True
                        break

                if (dostop):
                    break
                res.append(next_rep(cfg.first_rep_num + i,
                           cfg, nextepoch, e, r, f, v, p))

        print("Waiting for grompping to finish")
        p.close()
        p.join()

    for d, rc in res:
        if (rc.get()):
            raise NonzeroReturnError(
                "Nonzero returncode from grompp, "
                "see %s/output_grompp.txt for more detail." % (d),
                code=rc.get())

    # copy sbatch template
    utils.copy_sbatch_template(cfg.sbatch, pathlib.Path(
        "epoch%02d" % nextepoch) / cfg.sbatch.name, nextepoch, cfg)

    utils.copy_config(pathlib.Path(cfg.__file__),
                      pathlib.Path("epoch%02d" % nextepoch) / "config.py",
                      cfg.default_items)
