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

def get_data_from_archive(d, sel, sel_clust, function_val):
    with np.load(d+"fval_data.npz") as dat:
        if(dat["xtc_mod_t"]!=os.path.getmtime(d+"mdrun.xtc")):
            raise LoadError("Modification time of %s does not match, reloading"%(d+"mdrun.xtc"))
        if((not np.array_equal(sel.indices,dat["sel"])) or (not np.array_equal(sel_clust.indices,dat["sel_clust"]))):
            raise LoadError("Selections in %sfval_data.npz do not match, reloading"%d)
        if(dat["func_hash"]!=utils.hash_func(function_val)):
            print("Function hash changed, recalculating fval")
            fval= function_val(dat["fval_crd"])
            dat["fval"] = fval[-1]
            print("Saving modified fval to fval_data.npz")
            dat["func_hash"]=utils.hash_func(function_val)
            np.savez_compressed(d+"fval_data.npz",**dat)
        else:
            fval = dat["fval"]

        crd = dat["crd"]
    return fval, crd

def get_data_from_xtc(d, struct, sel, sel_clust, function_val, transforms):
    struct.load_new(d+"mdrun.xtc")
    struct.trajectory.add_transformations(*transforms)
    fval_crd = np.empty([len(struct.trajectory), len(sel), 3])
    crd      = np.empty([len(struct.trajectory), len(sel_clust), 3])
    print("Copying crd")
    for j,ts in enumerate(struct.trajectory):
        fval_crd[j] = sel.positions
        crd[j]      = sel_clust.positions

    print("Calculating fval")
    fval = function_val(fval_crd)
    print("Saving fval_data.npz")
    xtc_mod_t = os.path.getmtime(d+"mdrun.xtc")
    func_hash = utils.hash_func(function_val)
    np.savez_compressed(d+"fval_data.npz",
                        fval=fval,
                        fval_crd=fval_crd,
                        crd=crd,
                        xtc_mod_t=xtc_mod_t,
                        func_hash=func_hash,
                        sel=sel.indices, sel_clust=sel_clust.indices)

    return fval, crd

def load_from_dir(d, struct, sel, sel_clust, function_val, load_fval, transforms):
    if(not load_fval):
        try:
            fval,crd = get_data_from_archive(d, sel, sel_clust, function_val)
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

    return get_data_from_xtc(d, struct, sel, sel_clust, function_val, transforms)

def load_epoch_data(struct, sel, sel_clust, function_val, epoch, load_fval, transforms):
    N = check_num("epoch%02d/rep"%epoch)
    fval = []
    crd  = []
    for i in range(1,N+1):
        d="epoch%02d/rep%02d/"%(epoch, i)
        fv,cr = load_from_dir(d, struct, sel, sel_clust, function_val, load_fval, transforms)
        fval.append(fv)
        crd.append(cr)

    reps = [np.full(f.shape,i+1) for i,f in enumerate(fval)]
    frms = [np.arange(len(f)) for f in fval]

    return np.concatenate(reps), np.concatenate(fval), np.concatenate(frms), np.concatenate(crd)


def load_data(struct, sel, sel_clust, function_val, load_fval=False, transforms=[]):
    epochs = check_num("epoch")
    fval = []
    reps = []
    frms = []
    crds = []
    for i in range(1,epochs+1):
        r,f,fr,crd = load_epoch_data(struct, sel,sel_clust, function_val, i, load_fval, transforms)
        reps.append(r)
        fval.append(f)
        frms.append(fr)
        crds.append(crd)


    epcs = [np.full(f.shape,i+1) for i,f in enumerate(fval)]

    reps = np.concatenate(reps)
    fval = np.concatenate(fval)
    epcs = np.concatenate(epcs)
    frms = np.concatenate(frms)
    crds = np.concatenate(crds)
    return fval, crds, frms, reps, epcs
