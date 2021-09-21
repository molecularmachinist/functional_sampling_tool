# -*- coding: utf-8 -*-
import os

def check_num(prefix):
    """
    Checks filenames prefix01, prefix02, etc and returns the highest found (0 if none found)
    """
    parts = prefix.split("/")
    dir = "/".join(prefix.split("/")[:-1])
    # if dir is not empty, add trailing slash
    if(dir): dir+="/"
    fprefix = parts[-1]
    files = os.listdir(dir)
    i=0
    while(True):
        i+=1
        if("%s%02d"%(fprefix,i) not in files):
            return i-1

def load_epoch_data(struct, sel, epoch, load_fval=False):
    N = check_num("epoch%02d/rep")
    fval = []
    for i in range(1,N+1):
        d = "epoch%02d/rep%02d/"%(epoch, i)
        if(not load_fval):
            try:
                fval.append(np.load(d+"fval.npy"))
                print("Loaded fval from %sfval.npy"%d)
                continue
            except FileNotFoundError:
                print("Could not open %s"%(d+"fval.npy"))

        print("Loading trajectory %s"%(d+"mdrun.xtc"))
        traj = mdtraj.load(d+"mdrun.xtc", top=struct,atom_indices=sel)
        print("Calculating fval")
        fval.append(function_val(traj.xyz))
        print("Saving fval")
        np.save(d+"fval.npy", fval[-1])

    reps = [np.full(f.shape,i+1) for i,f in enumerate(fval)]
    frms = [np.arange(len(f)) for f in fval]

    return np.concatenate(reps), np.concatenate(fval), np.concatenate(frms)


def load_data():
    epochs = check_num("epoch")
    fval = []
    reps = []
    frms = []
    for i in range(1,epochs+1):
        r,f,fr = load_epoch_data(struct, sel, i, load_fval=False)
        reps.append(r)
        fval.append(f)
        frms.append(fr)


    epcs = [np.full(f.shape,i+1) for i,f in enumerate(fval)]

    reps = np.concatenate(reps)
    fval = np.concatenate(fval)
    epcs = np.concatenate(epcs)
    frms = np.concatenate(frms)
    return reps, fval, epcs, frms
