# -*- coding: utf-8 -*-
import os
import numpy as np
import mdtraj

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

def load_epoch_data(struct, sel, sel_clust, function_val, epoch, load_fval):
    N = check_num("epoch%02d/rep"%epoch)
    fval = []
    crd  = []
    for i in range(1,N+1):
        d = "epoch%02d/rep%02d/"%(epoch, i)
        if(not load_fval):
            try:
                with np.load(d+"fval_data.npz") as npz:
                    dat = dict(npz)
                if(dat["xtc_mod_t"]!=os.path.getmtime(d+"mdrun.xtc")):
                    raise LoadError("Modification time of %s does not match, reloading"%(d+"mdrun.xtc"))
                if((not np.array_equal(sel,dat["sel"])) or (not np.array_equal(sel_clust,dat["sel_clust"]))):
                    raise LoadError("Selections in %sfval_data.npz do not match, reloading"%d)
                if(dat["func_hash"]!=utils.hash_func(function_val)):
                    print(dat["func_hash"],utils.hash_func(function_val))
                    print("Function hash changed, recalculating fval")
                    fval.append(function_val(dat["fval_crd"]))
                    dat["fval"] = fval[-1]
                    print("Saving modified fval to fval_data.npz")
                    dat["func_hash"]=utils.hash_func(function_val)
                    np.savez_compressed(d+"fval_data.npz",**dat)
                else:
                    fval.append(dat["fval"])

                crd.append(dat["crd"])
                print("Loaded data from %sfval_data.npz"%d)
                continue
            except FileNotFoundError:
                print("Could not open %s, loading from xtc"%(d+"fval_data.npz"))
            except LoadError as e:
                print(e)
            except KeyError as e:
                print(d+"fval_data.npz"+" missing data:")
                print(e)
                print("reloading from xtc")

        print("Loading trajectory %s"%(d+"mdrun.xtc"))
        traj = mdtraj.load(d+"mdrun.xtc", top=struct)
        print("Calculating fval")
        fval_crd = traj.xyz[:,sel,:].copy()
        fval.append(function_val(fval_crd))
        print("Copying crd")
        crd.append(traj.xyz[:,sel_clust,:].copy())
        print("Saving fval_data.npz")
        xtc_mod_t = os.path.getmtime(d+"mdrun.xtc")
        func_hash = utils.hash_func(function_val)
        np.savez_compressed(d+"fval_data.npz",
                            fval=fval[-1],
                            fval_crd=fval_crd,
                            crd=crd[-1],
                            xtc_mod_t=xtc_mod_t,
                            func_hash=func_hash,
                            sel=sel, sel_clust=sel_clust)

    reps = [np.full(f.shape,i+1) for i,f in enumerate(fval)]
    frms = [np.arange(len(f)) for f in fval]

    return np.concatenate(reps), np.concatenate(fval), np.concatenate(frms), np.concatenate(crd)


def load_data(struct, sel, sel_clust, function_val, load_fval=False):
    epochs = check_num("epoch")
    fval = []
    reps = []
    frms = []
    crds = []
    for i in range(1,epochs+1):
        r,f,fr,crd = load_epoch_data(struct, sel,sel_clust, function_val, i, load_fval)
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
