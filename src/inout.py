# -*- coding: utf-8 -*-
import sys
import shutil
import time
import numpy as np
import re
import MDAnalysis as mda
import importlib
import pathlib

from . import utils
from . import transformations
from . import default_config

from .exceptions import NoConfigError

# Type hints
from numpy.typing import NDArray
from typing import Any, Tuple, Union, List, Dict


class LoadError(ValueError):
    pass


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


def get_data_from_archive(d: pathlib.Path, cfg: Any) -> Tuple[NDArray[np.float_], NDArray[np.float_]]:
    with np.load(d/cfg.npz_file_name) as npz:
        dat = dict(npz)

    if (dat["xtc_mod_t"] != (d/"mdrun.xtc").stat().st_mtime):
        raise LoadError(
            "Modification time of %s does not match, reloading" % str(d/"mdrun.xtc"))

    if ((not np.array_equal(cfg.sel.indices, dat["sel"])) or (not np.array_equal(cfg.sel_clust.indices, dat["sel_clust"]))):
        raise LoadError("Selections in %s do not match, reloading" %
                        str(d/cfg.npz_file_name))

    transform_opt = [cfg.unwrap_mols, cfg.unwrap_mols and cfg.mols_in_box,
                     cfg.clust_centre and not cfg.clust_superpos, cfg.clust_superpos]
    if (not np.array_equal(transform_opt, dat["transform_opt"])):
        raise LoadError("Trajectory transformations changed from %s, reloading" % str(
            d/cfg.npz_file_name))

    # if "stride" is found in dat, it was made with v0.0.1, so data might be missing
    if ("stride" in dat):
        raise LoadError("%s made with tool version 0.0.1 meaning data might be missing, reloading" % str(
            d/cfg.npz_file_name))

    unwrap_sel = cfg.traj_transforms[0].sel if cfg.unwrap_mols else np.zeros(
        0, dtype=int)
    if (not np.array_equal(unwrap_sel, dat["unwrap_sel"])):
        raise LoadError("Unwrapping selection in %s does not match, reloading" % str(
            d/cfg.npz_file_name))

    if (dat["func_hash"] != utils.hash_func(cfg.function_val)):
        print("Function hash changed, recalculating fval")
        fval = cfg.function_val(dat["fval_crd"])
        dat["fval"] = fval
        print("Saving modified fval to %s" % cfg.npz_file_name)
        dat["func_hash"] = utils.hash_func(cfg.function_val)
        np.savez_compressed(d/cfg.npz_file_name, **dat)
    else:
        fval = dat["fval"]

    crd = dat["crd"]
    return fval, crd


def get_data_from_xtc(d: pathlib.Path, cfg: Any) -> Tuple[NDArray[np.float_], NDArray[np.float_]]:
    cfg.struct.load_new(str(d / "mdrun.xtc"))
    cfg.struct.trajectory.add_transformations(*cfg.traj_transforms)
    trjlen = len(cfg.struct.trajectory)
    fval_crd = np.empty([trjlen, len(cfg.sel), 3])
    crd = np.empty([trjlen, len(cfg.sel_clust), 3])
    print("Copying crd")
    for j, ts in enumerate(cfg.struct.trajectory):
        fval_crd[j] = cfg.sel.positions
        ts = cfg.clust_transform(ts)
        crd[j] = cfg.sel_clust.positions

    print("Calculating fval")
    fval = cfg.function_val(fval_crd)
    print("Saving %s" % cfg.npz_file_name)
    xtc_mod_t = (d/"mdrun.xtc").stat().st_mtime
    func_hash = utils.hash_func(cfg.function_val)
    unwrap_sel = cfg.traj_transforms[0].sel if cfg.unwrap_mols else np.zeros(
        0, dtype=int)
    transform_opt = [cfg.unwrap_mols, cfg.unwrap_mols and cfg.mols_in_box,
                     cfg.clust_centre and not cfg.clust_superpos, cfg.clust_superpos]
    np.savez_compressed(d/cfg.npz_file_name,
                        fval=fval,
                        fval_crd=fval_crd,
                        crd=crd,
                        xtc_mod_t=xtc_mod_t,
                        func_hash=func_hash,
                        unwrap_sel=unwrap_sel,
                        transform_opt=transform_opt,
                        sel=cfg.sel.indices,
                        sel_clust=cfg.sel_clust.indices)

    return fval, crd


def load_from_dir(d: pathlib.Path, cfg: Any, load_fval: bool) -> Tuple[NDArray[np.float_], NDArray[np.float_]]:
    if (not load_fval):
        try:
            fval, crd = get_data_from_archive(d, cfg)
            print("Loaded data from", d / cfg.npz_file_name)
            return fval, crd
        except FileNotFoundError:
            print("Did not find %s, loading from xtc" %
                  str(d/cfg.npz_file_name))
        except LoadError as e:
            print(e)
        except KeyError as e:
            print(d/cfg.npz_file_name, "missing data:")
            print(e)
            print("reloading from xtc")

    print("Loading trajectory", d / "mdrun.xtc")

    return get_data_from_xtc(d, cfg)


def load_epoch_data(epoch: int, cfg: Any, load_fval: bool) -> Tuple[
        NDArray[np.int_], NDArray[np.float_], NDArray[np.int_], NDArray[np.float_]]:
    edir = pathlib.Path("epoch%02d" % epoch)
    rep_nums = check_num(edir / "rep")
    fval = []
    crd = []
    reps = []
    for i in rep_nums:
        if ((epoch, i) in cfg.ignore_reps):
            continue
        d = edir / ("rep%02d" % i)
        if ((d/"mdrun.xtc").exists()):
            fv, cr = load_from_dir(d, cfg, load_fval)
            if (cfg.ignore_from_start != 0 or cfg.stride != 1):
                # If we are ignoring data, we'll make a copy so as to not drag the ignored data in memory
                # Otherwise numpy will just make a view of the original array and keep it in memory
                fv = fv[cfg.ignore_from_start::cfg.stride].copy()
                cr = cr[cfg.ignore_from_start::cfg.stride].copy()

            fval.append(fv)
            crd.append(cr)
            reps.append(np.full(fval[-1].shape, i))
        else:
            print(f"No data loaded from {d}")

    # if nothing was found, we continue with no error
    if (not fval):
        print(f"No data loaded from {edir}")
        return np.zeros(0, dtype=int), np.zeros(0, dtype=float), np.zeros(0, dtype=int), np.zeros((0, len(cfg.sel_clust), 3), dtype=float)

    frms = [np.arange(len(f))*cfg.stride+cfg.ignore_from_start for f in fval]

    return np.concatenate(reps), np.concatenate(fval), np.concatenate(frms), np.concatenate(crd)


def load_extract_data(cfg: Any, doignore: bool = True) -> Dict[str, Dict[int, Dict[int, Union[str, NDArray[np.float_]]]]]:
    epochs = check_num(pathlib.Path("epoch"))
    data = {"fval": {}, "fnames": {}}
    for e in epochs:
        if (doignore and (e in cfg.ignore_epcs)):
            continue

        for key in data:
            data[key][e] = {}

        edir = pathlib.Path("epoch%02d" % e)
        rep_nums = check_num(edir / "rep")

        for r in rep_nums:
            if (doignore and ((e, r) in cfg.ignore_reps)):
                continue
            d = edir / ("rep%02d" % r)
            if (not (d/cfg.npz_file_name).is_file()):
                continue

            # mdrun filename
            data["fnames"][e][r] = str(d / "mdrun.xtc")
            with np.load(d/cfg.npz_file_name) as npz:
                data["fval"][e][r] = npz["fval"]

    return data


def load_flat_extract_data(cfg: Any, doignore: bool = True) -> Tuple[
        Dict[str, NDArray[Union[np.float_, np.int_]]], Dict[int, Dict[int, str]]]:
    data = load_extract_data(cfg, doignore)

    flat_data = {"fval": [], "frms": [], "reps": [], "epcs": []}
    for e in data["fval"]:
        for r in data["fval"][e]:
            flat_data["fval"].append(data["fval"][e][r])
            ndat = len(data["fval"][e][r])
            flat_data["frms"].append(np.arange(ndat))
            flat_data["reps"].append(np.full(ndat, r))
            flat_data["epcs"].append(np.full(ndat, e))

    for key in flat_data:
        flat_data[key] = np.concatenate(flat_data[key])

    return flat_data, data["fnames"]


def load_data(cfg: Any, load_fval: bool):
    epochs = check_num(pathlib.Path("epoch"))
    if (not epochs):
        raise FileNotFoundError(
            "No epochs found, cannot load data. Are you sure you are in the correct folder?")
    fval = []
    reps = []
    frms = []
    crds = []
    epcs = []
    for i in epochs:
        if (i in cfg.ignore_epcs):
            continue
        r, f, fr, crd = load_epoch_data(i, cfg, load_fval)
        reps.append(r)
        fval.append(f)
        frms.append(fr)
        crds.append(crd)
        epcs.append(np.full(f.shape, i))

    if (not fval):
        raise FileNotFoundError(
            "No data loaded, cannot continue. Run some more simulations.")

    reps = np.concatenate(reps)
    fval = np.concatenate(fval)
    epcs = np.concatenate(epcs)
    frms = np.concatenate(frms)
    crds = np.concatenate(crds)
    return fval, crds, frms, reps, epcs


def import_cfg(cfgpath: pathlib.Path) -> Any:
    """
    Only import the config, does not load structs or do anything with it.
    """
    if (not cfgpath.exists()):
        raise NoConfigError("Config file does not exist at %s." % cfgpath)
    print("Loading config from %s" % cfgpath)
    spec = importlib.util.spec_from_file_location("config", cfgpath)
    writebytecode = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    cfg = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cfg)
    sys.dont_write_bytecode = writebytecode

    # Load defaults
    cfg.default_items = {}
    for item in dir(default_config):
        if (item.startswith("__")):
            continue
        if (not hasattr(cfg, item)):
            value = getattr(default_config, item)
            cfg.default_items[item] = value
            setattr(cfg, item, value)

    cfg.ignore_epcs = set(cfg.ignore_epcs)
    cfg.ignore_reps = set(cfg.ignore_reps)

    if (type(cfg.initial_struct) != list):
        cfg.initial_struct = [cfg.initial_struct]

    cfg.initial_struct = [pathlib.Path(item) for item in cfg.initial_struct]

    # Make path variables into Path objects
    for pathvar in ("npz_file_name", "fig_output_dir", "mdp", "topol", "ndx", "sbatch"):
        aspath = pathlib.Path(getattr(cfg, pathvar))
        setattr(cfg, pathvar, aspath)

    return cfg


def load_options(cfgpath: pathlib.Path) -> Any:
    """
    Imports the config, and loads the structs and selections
    """
    cfg = import_cfg(cfgpath)
    print("Loading structure")
    cfg.struct = mda.Universe(str(cfg.initial_struct[0]))
    cfg.indexes = utils.read_ndx(cfg.index_file)
    cfg.sel = utils.load_sel(cfg.select_str, cfg.struct, cfg.indexes)
    cfg.sel_clust = utils.load_sel(
        cfg.select_str_clust, cfg.struct, cfg.indexes)

    print("Selected %d atoms" % len(cfg.sel))
    print("Selected %d atoms for clustering" % len(cfg.sel_clust))

    if (cfg.unwrap_mols):
        # Preparing molecule unwrapper
        mdrunpath = pathlib.Path("epoch01")/"rep01"/"mdrun.tpr"
        bonded_struct = mda.Universe(str(mdrunpath),
                                     str(cfg.initial_struct[0]))
        unwrap_sel = utils.load_sel(cfg.unwrap_sel, cfg.struct, cfg.indexes)
        unwrap_sel = bonded_struct.atoms[unwrap_sel.indices]
        if (cfg.unwrap_starters is None):
            unwrap_starters = []
        else:
            unwrap_starters = utils.load_sel(
                cfg.unwrap_starters, unwrap_sel, cfg.indexes)
        cfg.traj_transforms = [
            transformations.Unwrapper(unwrap_sel, unwrap_starters)]
        print("Selected %d atoms for unwrapping" % len(unwrap_sel))
        if (cfg.mols_in_box):
            print("Also putting mol COMs back in box")
            cfg.traj_transforms.append(transformations.MolWrapper(unwrap_sel))

    else:
        cfg.traj_transforms = []
    cfg.startval = cfg.function_val(np.array([cfg.sel.positions]))[0]
    print("Initial function value %g" % cfg.startval)

    if (cfg.clust_superpos):
        print("Clustering coordinates will be superpositioned")
    elif (cfg.clust_centre):
        print("Clustering coordinates will be centred")

    cfg.clust_transform = transformations.Superpos(cfg.sel_clust,
                                                   cfg.clust_centre,
                                                   cfg.clust_superpos)

    # Min and maxvals
    if (cfg.minval == "start"):
        cfg.minval = cfg.startval
    elif (cfg.minval is None):
        cfg.minval = -float("inf")

    if (cfg.maxval == "start"):
        cfg.maxval = cfg.startval
    elif (cfg.maxval is None):
        cfg.maxval = float("inf")

    # Make the global rng
    cfg.rng = np.random.default_rng(cfg.rng_seed)

    return cfg


def clean_latest_epoch(force: bool = False) -> None:
    epochs = check_num(pathlib.Path("epoch"))
    if (not epochs):
        raise FileNotFoundError(
            "No epochs found, are you sure you are in the correct folder?")
    epoch = epochs[-1]

    # First, make sure no rep is already simulated
    edir = pathlib.Path("epoch%02d" % epoch)
    for r in check_num(edir / "rep"):
        mdrun_path = edir / ("rep%02d" % r) / "mdrun.xtc"
        if (mdrun_path.is_file()):
            if (not force):
                print(f"Found {mdrun_path}, will not clean up epoch {epoch}")
                print("If you do want it cleaned, use the --force option")
                return

            print(f"Found {mdrun_path}, but will still force clean")
            print("Data in epoch%02d will be lost forever" % epoch)
            print("You have 5 seconds to cancel with ctrl-C")
            time.sleep(5)
            print("Continuing, data in epoch%02d will now be lost forever" % epoch)
            break

    print(f"Cleaning up epoch {epoch}")
    # If nothing was found (or force=True), we can just do it
    shutil.rmtree(edir)
