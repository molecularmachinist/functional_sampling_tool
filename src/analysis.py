import os
import MDAnalysis as mda
import argparse
import pathlib
import numpy as np
import warnings

from . import inout
from . import utils
from . import transformations

# Type hints
from typing import Any, Tuple, Optional, List, Dict
from numpy.typing import ArrayLike
from MDAnalysis.core.groups import AtomGroup
from .utils import transform_type


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

    args.selection = get_cfg_sel(args.selection,    args.cfg)
    args.sel_unwrap = get_cfg_sel(args.sel_unwrap,   args.cfg, args.selection)
    args.sel_superpos = get_cfg_sel(
        args.sel_superpos, args.cfg, args.selection)

    if (args.unwrap_starters == "unwrap_starters"):
        args.unwrap_starters = args.cfg.unwrap_starters

    print("Loading structure")
    u = mda.Universe(str(args.cfg.initial_struct))
    sel = utils.load_sel(args.selection, u, indexes)
    print("Selected %d atoms for extraction" % len(sel))
    sel_superpos = utils.load_sel(args.sel_superpos, u, indexes)

    if (args.unwrap or args.wrap):
        # Preparing molecule unwrapper
        bonded_struct = mda.Universe(
            "epoch01/rep01/mdrun.tpr", str(args.cfg.initial_struct))
        unwrap_sel = utils.load_sel(args.sel_unwrap, u, indexes)
        unwrap_sel = bonded_struct.atoms[unwrap_sel.indices]
        print("Selected %d atoms for unwrapping" % len(unwrap_sel))
        if (args.unwrap_starters is None):
            unwrap_starters = []
        else:
            unwrap_starters = utils.load_sel(
                args.unwrap_starters, unwrap_sel, indexes)
            print("Selected %d atoms as unwrap starters" %
                  len(unwrap_starters))

    traj_transforms = []
    if (args.unwrap):
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


def analysis_subparser(parser: argparse.ArgumentParser) -> None:
    """
    Add the analysis options to the provided parser
    """
    subparsers = parser.add_subparsers(
        description="Specify which analysis function to run")
    subparsers.required = True
    subparsers.dest = 'subcommand'

    extract_description = "Only repetitions for which the data " \
        "npz-archive is present can be extracted. Run the choose command with " \
        "the--choose_only option to calculate the functions and make the archives " \
        "without making a new epoch"
    # extractor command
    extr_parser = subparsers.add_parser("extract", help="Extract frames.",
                                        description=extract_description)
    extr_parser.add_argument("--selection", metavar="str",
                             help="The selection to extract. \"select_str\" and \"select_str_clust\" will be read from the config. "
                                  "as with those selection string, this can also be a group in \"index_file\" (default: \"%(default)s\")",
                             default="select_str")
    extr_parser.add_argument("-n", "--index", metavar="<name>.ndx", dest="index",
                             help="Overwrite the value in \"index_file\" to use another one just for this analysis (default: %(default)s)",
                             default=None)
    extr_parser.add_argument("-o", "--output", metavar="<name>.xtc", dest="output",
                             help="Output xtc file for extraction. Another file with the same name, but npz file ending will be made with info on each extracted frame, as well as a structure file. (default: %(default)s)",
                             default=pathlib.Path("analysis/extracted.xtc"),
                             type=pathlib.Path)
    extr_parser.add_argument("--noignore",     action="store_false", dest="doignore",
                             help="Do not ignore the epochs and repetitions defined in the config (default: do ignore)")
    extr_parser.add_argument("--around", metavar="fval",
                             help="A function value around which to extract frames from. Extracted frames will be ordered by function value. "
                                  "By default None, which means all frames are extracted, ordered by epoch, rep and frame.",
                             default=None,
                             type=float)
    extr_parser.add_argument("--number", "-N", metavar="number",
                             help="The number of closest frames to extract around the --around value. Ignored if --around not given. (default: %(default)d)",
                             default=1000,
                             type=int)
    extr_parser.add_argument("--stride", metavar="N",
                             help="Only extract every Nth frame. By default \"stride\" from config.",
                             default=None,  type=int)
    extr_parser.add_argument("--beginning", "-b", metavar="N",
                             help="Ignore this many frames from the beginning of each simulation. By default \"ignore_from_start\" from config.",
                             default=None,  type=int)
    extr_parser.add_argument("--unwrap",   action="store_true",
                             help="Unwrap molecules (default: No)")
    extr_parser.add_argument("--sel_unwrap", metavar="str",
                             help="The selection to unwrap. Same syntax as --selection, and by default uses same selection.",
                             default=None)
    extr_parser.add_argument("--unwrap_starters", metavar="str",
                             help="The selection to use as starters in unwrapping. \"unwrap_starters\" to use from config, by default use none",
                             default=None)
    extr_parser.add_argument("--wrap",     action="store_true",
                             help="Put centre of mass of molecules back in box. "
                                  "Uses the --sel_unwrap. (default: No)")
    extr_parser.add_argument("--superpos_trans", action="store_true",
                             help="Superposition the --sel_superpos selection onto the initial structure, only translating (default: No)")
    extr_parser.add_argument("--superposition",     action="store_true",
                             help="Superposition the --sel_superpos selection onto the initial structure, translating and rotating (default: No)")
    extr_parser.add_argument("--sel_superpos", metavar="str",
                             help="The selection to superposition. Same syntax as --selection, and by default uses same selection.",
                             default=None)
    # Use just config import as config_func, so the selections and such are not loaded
    extr_parser.set_defaults(func=extract, config_func=inout.import_cfg)


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
        (args.output.stem + args.cfg.initial_struct.suffix)
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
