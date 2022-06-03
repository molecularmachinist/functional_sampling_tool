# -*- coding: utf-8 -*-
import os
import numpy as np

from . import utils

class LoadError(ValueError):
    pass

def check_num(prefix):
    """
    Checks filenames prefix01, prefix02, etc and returns the highest found (0 if none found)
    """
    parts = prefix.split("/")
    dir = "/".join(prefix.split("/")[:-1])
    # if dir is not empty, add trailing slash
    if(dir):
        dir+="/"
    else:
        dir="./"
    fprefix = parts[-1]
    files = os.listdir(dir)
    i=0
    while(True):
        i+=1
        if("%s%02d"%(fprefix,i) not in files):
            return i-1

def get_data_from_archive(d, cfg):
    with np.load(d+"fval_data.npz") as npz:
        dat = dict(npz)

    if(dat["xtc_mod_t"]!=os.path.getmtime(d+"mdrun.xtc")):
        raise LoadError("Modification time of %s does not match, reloading"%(d+"mdrun.xtc"))

    if((not np.array_equal(cfg.sel.indices,dat["sel"])) or (not np.array_equal(cfg.sel_clust.indices,dat["sel_clust"]))):
        raise LoadError("Selections in %sfval_data.npz do not match, reloading"%d)

    transform_opt = [cfg.unwrap_mols, cfg.unwrap_mols and cfg.mols_in_box,
                     cfg.clust_centre and not cfg.clust_superpos, cfg.clust_superpos]
    if(not np.array_equal(transform_opt,dat["transform_opt"])):
        raise LoadError("Trajetory transformations changed from %sfval_data.npz, reloading"%d)

    unwrap_sel = cfg.traj_transforms[0].sel if cfg.unwrap_mols else np.zeros(0,dtype=int)
    if(not np.array_equal(unwrap_sel,dat["unwrap_sel"])):
        raise LoadError("Unwrapping selection in %sfval_data.npz does not match, reloading"%d)

    if(dat["func_hash"]!=utils.hash_func(cfg.function_val)):
        print("Function hash changed, recalculating fval")
        fval= cfg.function_val(dat["fval_crd"])
        dat["fval"] = fval
        print("Saving modified fval to fval_data.npz")
        dat["func_hash"]=utils.hash_func(cfg.function_val)
        np.savez_compressed(d+"fval_data.npz",**dat)
    else:
        fval = dat["fval"]

    crd = dat["crd"]
    return fval, crd

def get_data_from_xtc(d, cfg):
    cfg.struct.load_new(d+"mdrun.xtc")
    cfg.struct.trajectory.add_transformations(*cfg.traj_transforms)
    fval_crd = np.empty([len(cfg.struct.trajectory), len(cfg.sel), 3])
    crd      = np.empty([len(cfg.struct.trajectory), len(cfg.sel_clust), 3])
    print("Copying crd")
    for j,ts in enumerate(cfg.struct.trajectory):
        fval_crd[j] = cfg.sel.positions
        ts = cfg.clust_transform(ts)
        crd[j]      = cfg.sel_clust.positions

    print("Calculating fval")
    fval = cfg.function_val(fval_crd)
    print("Saving fval_data.npz")
    xtc_mod_t = os.path.getmtime(d+"mdrun.xtc")
    func_hash = utils.hash_func(cfg.function_val)
    unwrap_sel = cfg.traj_transforms[0].sel if cfg.unwrap_mols else np.zeros(0,dtype=int)
    transform_opt = [cfg.unwrap_mols, cfg.unwrap_mols and cfg.mols_in_box,
                     cfg.clust_centre and not cfg.clust_superpos, cfg.clust_superpos]
    np.savez_compressed(d+"fval_data.npz",
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

def load_from_dir(d, cfg, load_fval):
    if(not load_fval):
        try:
            fval,crd = get_data_from_archive(d, cfg)
            print("Loaded data from %sfval_data.npz"%d)
            return fval, crd
        except FileNotFoundError:
            print("Did not find %s, loading from xtc"%(d+"fval_data.npz"))
        except LoadError as e:
            print(e)
        except KeyError as e:
            print(d+"fval_data.npz"+" missing data:")
            print(e)
            print("reloading from xtc")

    print("Loading trajectory %s"%(d+"mdrun.xtc"))

    return get_data_from_xtc(d, cfg)

def load_epoch_data(epoch, cfg, load_fval):
    N = check_num("epoch%02d/rep"%epoch)
    fval = []
    crd  = []
    reps = []
    for i in range(1,N+1):
        if((epoch,i) in cfg.ignore_reps):
            continue
        d="epoch%02d/rep%02d/"%(epoch, i)
        fv,cr = load_from_dir(d, cfg, load_fval)
        fval.append(fv)
        crd.append(cr)
        reps.append(np.full(fv.shape,i))

    frms = [np.arange(len(f)) for f in fval]

    return np.concatenate(reps), np.concatenate(fval), np.concatenate(frms), np.concatenate(crd)


def load_data(cfg, load_fval):
    epochs = check_num("epoch")
    fval = []
    reps = []
    frms = []
    crds = []
    epcs = []
    for i in range(1,epochs+1):
        if (i in cfg.ignore_epcs):
            continue
        r,f,fr,crd = load_epoch_data(i, cfg, load_fval)
        reps.append(r)
        fval.append(f)
        frms.append(fr)
        crds.append(crd)
        epcs.append(np.full(f.shape,i))


    reps = np.concatenate(reps)
    fval = np.concatenate(fval)
    epcs = np.concatenate(epcs)
    frms = np.concatenate(frms)
    crds = np.concatenate(crds)
    return fval, crds, frms, reps, epcs
