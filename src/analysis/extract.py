import MDAnalysis as mda
import argparse
import numpy as np
import warnings

from .. import inout
from .. import utils
from .. import transformations
from ..exceptions import WrongSelectionSizeError

# Type hints
from typing import Any, Tuple, Optional, List, Dict
from numpy.typing import ArrayLike
from MDAnalysis.core.groups import AtomGroup
from ..utils import transform_type


def get_cfg_sel(sel: str, cfg: Any, default: Optional[str] = None) -> str:
    if (sel is None):
        return default
    if (sel == "select_str"):
        return cfg.select_str
    elif (sel == "select_str_clust"):
        return cfg.select_str_clust
    return sel


def load_struct(args: argparse.Namespace) -> Tuple[mda.Universe, AtomGroup, List[transform_type]]:
    if (args.index is None):
        args.index = args.cfg.index_file
    indexes = utils.read_ndx(args.index)

    args.selection = get_cfg_sel(args.selection, args.cfg)
    args.sel_unwrap = get_cfg_sel(args.sel_unwrap, args.cfg, args.selection)
    args.sel_superpos = get_cfg_sel(args.sel_superpos, args.cfg,
                                    args.selection)
    args.precenter_atom = get_cfg_sel(args.precenter_atom, args.cfg)

    if (args.unwrap_starters == "unwrap_starters"):
        args.unwrap_starters = args.cfg.unwrap_starters

    print("Loading structure")
    u = mda.Universe(args.cfg.initial_struct[0])
    sel = utils.load_sel(args.selection, u, indexes)
    print("Selected %d atoms for extraction" % len(sel))
    sel_superpos = utils.load_sel(args.sel_superpos, u, indexes)

    if (args.unwrap or args.wrap or args.precenter):
        unwrap_sel = utils.load_sel(args.sel_unwrap, u, indexes)
    if (args.unwrap or args.wrap):
        # Preparing molecule unwrapper
        bonded_struct = mda.Universe(
            "epoch01/rep01/mdrun.tpr",
            args.cfg.initial_struct[0]
        )
        unwrap_sel = bonded_struct.atoms[unwrap_sel.indices]

    traj_transforms = []
    if (args.precenter):
        if (args.precenter_atom is not None):
            centre_atom_group = utils.load_sel(args.precenter_atom,
                                               u, indexes)
            if (len(centre_atom_group) != 1):
                raise WrongSelectionSizeError(f"Selection precenter_atom={repr(args.precenter_atom)} resulted "
                                              f"in {len(centre_atom_group)} atoms. Should be exactly 1!")
            centre_atom = centre_atom_group[0]
        else:
            centre_atom = None
        traj_transforms.append(
            transformations.Precenter(unwrap_sel, centre_atom))
        ca = u.atoms[traj_transforms[-1].centre_atom]
        print(f"Precentering using atom index {ca.index}",
              f"({ca.name} of {ca.resname}:{ca.resid})")

    if (args.unwrap):
        print("Selected %d atoms for unwrapping" % len(unwrap_sel))
        if (args.unwrap_starters is None):
            unwrap_starters = []
        else:
            unwrap_starters = utils.load_sel(
                args.unwrap_starters, unwrap_sel, indexes)
            print("Selected %d atoms as unwrap starters" %
                  len(unwrap_starters))
        traj_transforms.append(
            transformations.Unwrapper(unwrap_sel, unwrap_starters))

    if (args.wrap):
        print("Putting mol COMs back in box for selection of %d atoms" %
              len(unwrap_sel))
        traj_transforms.append(transformations.MolWrapper(unwrap_sel))

    if (args.superposition):
        print("%d atoms will be superpositioned" % (len(sel_superpos)))
    elif (args.superpos_trans):
        print("%d atoms coordinates will be centred on " % (len(sel_superpos)))

    traj_transforms.append(transformations.Superpos(sel_superpos,
                                                    args.superpos_trans,
                                                    args.superposition,
                                                    subselection=sel))

    return u, sel, traj_transforms


def extract_around(u: mda.Universe, sel: AtomGroup, transforms: List[transform_type], args: argparse.Namespace) -> Dict[str, ArrayLike]:
    data, files = inout.load_flat_extract_data(args.cfg, args.doignore)
    mask = (data["frms"] >= args.beginning)*(data["frms"] % args.stride == 0)

    srt_ndx = np.arange(len(data["fval"]))[mask][np.argsort(
        np.abs(data["fval"][mask]-args.around))][:args.number]
    for key in data:
        data[key] = data[key][srt_ndx]

    data["time"] = []
    prev_e = -1
    prev_r = -1
    ndata = len(srt_ndx)
    print("Writing frames to", args.output)
    with mda.Writer(str(args.output), sel.n_atoms) as writer:
        for i, (e, r, f) in enumerate(zip(data["epcs"], data["reps"], data["frms"])):
            print("Writing frame %d/%d" % (i+1, ndata), end="\r")
            if (prev_e != e or prev_r != r):
                u.load_new(files[e][r])
                u.trajectory.add_transformations(*transforms)
                prev_r = r
                prev_e = e

            ts = u.trajectory[f]
            data["time"].append(ts.time)
            writer.write(sel)

    data["time"] = np.array(data["time"])
    return data


def extract_all(u: mda.Universe, sel: AtomGroup, transforms: List[transform_type], args: argparse.Namespace) -> Dict[str, ArrayLike]:
    data = inout.load_extract_data(args.cfg, args.doignore)
    data_out = {"frame": [], "time": [], "epoch": [], "rep": [], "fval": []}
    ntrajs = np.sum([len(data["fnames"][e]) for e in data["fnames"]])

    i = 0
    print("Writing frames to", args.output)
    with mda.Writer(str(args.output), sel.n_atoms) as writer:
        for e in data["fval"]:
            for r in data["fval"][e]:
                i += 1
                print("Writing trajectory %d/%d" % (i, ntrajs), end="\r")
                u.load_new(data["fnames"][e][r])
                u.trajectory.add_transformations(*transforms)
                for ts in u.trajectory[args.beginning::args.stride]:
                    writer.write(sel)
                    data_out["frame"].append(ts.frame)
                    data_out["time"].append(ts.time)
                    data_out["epoch"].append(e)
                    data_out["rep"].append(r)
                    data_out["fval"].append(data["fval"][e][r][ts.frame])

    return data_out


def extract(args: argparse.Namespace) -> None:
    if (args.stride is None):
        args.stride = args.cfg.stride
    if (args.beginning is None):
        args.beginning = args.cfg.ignore_from_start

    u, sel, transforms = load_struct(args)

    u.trajectory.add_transformations(*transforms)
    # The pdb file name from xtc file dir and stem (name without suffix)
    struct_out = args.output.parent / \
        (args.output.stem + args.cfg.initial_struct[0].suffix)
    #
    args.output.parent.mkdir(parents=True, exist_ok=True)
    print("Writing structure to", struct_out)
    # Silencing warning about missing chainIDs
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", "Found missing chainIDs.", UserWarning)
        warnings.filterwarnings(
            "ignore", "Found no information for attr", UserWarning)
        sel.write(struct_out)

    if (args.around is None):
        data = extract_all(u, sel, transforms, args)
    else:
        data = extract_around(u, sel, transforms, args)

    npz_out = args.output.parent / (args.output.stem + ".npz")
    print("\nWriting data to", npz_out)
    np.savez_compressed(npz_out, **data)
