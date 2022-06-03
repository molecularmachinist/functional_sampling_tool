#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import MDAnalysis as mda
import importlib
import argparse
import pathlib
import sys,os
import numpy as np

from . import inout
from . import epoch_starting
from . import clustering
from . import utils
from . import transformations
from . import default_config

def import_cfg(cfgname):
    """
    Only import the config, does not load structs or do anything with it.
    """
    print("Loading config from %s"%cfgname)
    spec = importlib.util.spec_from_file_location("config", cfgname)
    writebytecode = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    cfg = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cfg)
    sys.dont_write_bytecode = writebytecode

    #Load defaults
    for item in dir(default_config):
        if(item.startswith("__")):
            continue
        if(not hasattr(cfg, item)):
            setattr(cfg,item,getattr(default_config,item))

    cfg.ignore_epcs = set(cfg.ignore_epcs)
    cfg.ignore_reps = set(cfg.ignore_reps)
    return cfg

def load_sel(sel_str, struct, ndx):
    if(sel_str in ndx):
        sel = struct.atoms[np.array(ndx[sel_str])-1]
    else:
        sel = struct.select_atoms(sel_str)
    return sel

def load_options(cfgname):
    """
    Imports the config, and loads the structs and selections
    """
    cfg    = import_cfg(cfgname)
    print("Loading structure")
    if(not os.path.isfile("initial/start.pdb")):
        utils.make_pdb(cfg)
    cfg.struct = mda.Universe("initial/start.pdb")
    if(not cfg.index_file is None):
        cfg.indexes = utils.read_ndx(cfg.index_file)
    else:
        cfg.indexes = {}
    cfg.sel = load_sel(cfg.select_str, cfg.struct, cfg.indexes)
    cfg.sel_clust = load_sel(cfg.select_str_clust, cfg.struct, cfg.indexes)

    print("Selected %d atoms"%len(cfg.sel))
    print("Selected %d atoms for clustering"%len(cfg.sel_clust))

    if(cfg.unwrap_mols):
        # Preparing molecule unwrapper
        bonded_struct = mda.Universe("epoch01/rep01/mdrun.tpr", "initial/start.pdb")
        unwrap_sel = load_sel(cfg.unwrap_sel, cfg.struct, cfg.indexes)
        unwrap_sel = bonded_struct.atoms[unwrap_sel.indices]
        if(cfg.unwrap_starters is None):
            unwrap_starters = []
        else:
            unwrap_starters = load_sel(cfg.unwrap_starters, unwrap_sel, cfg.indexes)
        cfg.traj_transforms = [transformations.Unwrapper(unwrap_sel,unwrap_starters)]
        print("Selected %d atoms for unwrapping"%len(unwrap_sel))
        if(cfg.mols_in_box):
            print("Also putting mol COMs back in box")
            cfg.traj_transforms.append(transformations.MolWrapper(unwrap_sel))

    else:
        cfg.traj_transforms = []
    cfg.startval = cfg.function_val(np.array([cfg.sel.positions]))[0]
    print("Initial function value %g"%cfg.startval)

    if(cfg.clust_superpos):
        print("Clustering coordinates will be superpositioned")
    elif(cfg.clust_centre):
        print("Clustering coordinates will be centred")

    cfg.clust_transform = transformations.Superpos(cfg.sel_clust,
                                                   cfg.clust_centre,
                                                   cfg.clust_superpos)


    # Min and maxvals
    if(cfg.minval=="start"):
        cfg.minval=cfg.startval
    elif(cfg.minval is None):
        cfg.minval=-float("inf")

    if(cfg.maxval=="start"):
        cfg.maxval=cfg.startval
    elif(cfg.maxval is None):
        cfg.maxval=float("inf")

    return cfg


def init(args):
    print("Initializing first epoch")
    epoch_starting.start_epoch(1, args.cfg)
    if(args.push):
        utils.rsync_up(args.cfg)


def pushpull(args):
    if(args.pull):
        utils.rsync_down(args.cfg)
    else:
        utils.rsync_up(args.cfg)


def choose(args):
    if(args.pull):
        utils.rsync_down(args.cfg)
    chooser = clustering.ClusterChooser.fromReadData(args.cfg, args.reload_fval)
    chooser.plot_hist()
    val, epc, rep, frm = chooser.make_choices()
    if(not args.choose_only):
        epoch_starting.start_epoch(chooser.nextepoch, args.cfg, val, epc, rep, frm)
    else:
        chooser.print_choices(val, epc, rep, frm)
    if(args.push):
        utils.rsync_up(args.cfg)


def argP():
    """
    Define command line arguments. For each
    """
    description = "Functional Sampling Tool: A tool for functional sampling."
    epilog = "For more detailed explanations, see the README.md"
    parser = argparse.ArgumentParser(prog="fst", description=description, epilog=epilog)

    # The only global option
    parser.add_argument("-c","--config",metavar="<name>.py",default=pathlib.Path("config.py"),
                        type=pathlib.Path, help="Path to config file (default: %(default)s)")
    # Global default for config_func
    parser.set_defaults(config_func=load_options)

    # A subparser object to make subcommands. A subcommand is required.
    subparsers = parser.add_subparsers(description="Specify what the program should do. Only one command should be given per run. 'fst <cmd> -h' to get command specific option.")
    subparsers.required = True
    subparsers.dest = 'command'

    # First epoch initializer command
    init_parser = subparsers.add_parser("init", help="Initialize the first epoch.")
    init_parser.add_argument("--push", action="store_true", help="push to remote after initialization (default: %(default)s)")
    # The initializer has a different config func to not try to load the struct.
    init_parser.set_defaults(func=init, config_func=import_cfg)

    # The chooser command
    choose_parser = subparsers.add_parser("choose", help="Choose frames for next epoch and initialize it.")
    choose_parser.add_argument("--pull", action="store_true", help="pull from remote before choosing (default: %(default)s)")
    choose_parser.add_argument("--push", action="store_true", help="push to remote after choosing (default: %(default)s)")
    choose_parser.add_argument("--choose_only", action="store_true", help="Only make the choices and plots, do not initialize next epoch (default: %(default)s)")
    choose_parser.add_argument("--reload_fval", action="store_true", help="Reload data from mdrun.xtc even if fval_data.npz exists and everything matches (default: %(default)s)")
    choose_parser.set_defaults(func=choose)

    # newepoch command
    epochstarter_parser = subparsers.add_parser("newepoch", help="Shorthand for \"choose --pull --push\".")
    epochstarter_parser.set_defaults(func=choose, push=True, pull=True, choose_only=False)
    epochstarter_parser.add_argument("--reload_fval", action="store_true", help="Reload data from mdrun.xtc even if fval_data.npz exists and everything matches (default: %(default)s)")

    # Push and pull commands
    push_parser = subparsers.add_parser("push", help="rsync from local to remote")
    push_parser.set_defaults(func=pushpull,config_func=import_cfg,pull=False)
    pull_parser = subparsers.add_parser("pull", help="rsync from remote to local")
    pull_parser.set_defaults(func=pushpull,config_func=import_cfg,pull=True)


    arguments = parser.parse_args()
    assert arguments.config.exists(), "Config file does not exist at %s."%arguments.config

    arguments.cfg = arguments.config_func(arguments.config)

    return arguments
