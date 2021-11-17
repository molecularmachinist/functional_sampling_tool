#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import mdtraj
import importlib
import argparse
import pathlib
import sys,os
import numpy as np

from . import inout
from . import epoch_starting
from . import clustering
from . import utils

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
    return cfg


def load_options(cfgname):
    """
    Imports the config, and loads the structs and selections
    """
    cfg    = import_cfg(cfgname)
    print("Loading structure")
    if(not os.path.isfile("initial/start.pdb")):
        util.make_pdb(cfg)
    cfg.struct = mdtraj.load("initial/start.pdb")
    if(not cfg.index_file is None):
        cfg.indexes = utils.read_ndx(cfg.index_file)
    else:
        cfg.indexes = {}
    if(cfg.select_str in cfg.indexes):
        cfg.sel = np.array(cfg.indexes[cfg.select_str])-1
    else:
        cfg.sel = cfg.struct.topology.select(cfg.select_str)-1
    if(cfg.select_str_clust in cfg.indexes):
        cfg.sel_clust = np.array(cfg.indexes[cfg.select_str_clust])
    else:
        cfg.sel_clust = cfg.struct.topology.select(cfg.select_str_clust)
    print("Selected %d atoms"%len(cfg.sel))
    print("Selected %d atoms for clustering"%len(cfg.sel_clust))
    cfg.startval = cfg.function_val(cfg.struct.xyz[:,cfg.sel,:])[0]
    print("Initial function value %g"%cfg.startval)

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
    choose_parser.add_argument("--pull", action="store_true", help="push to remote after initialization (default: %(default)s)")
    choose_parser.add_argument("--push", action="store_true", help="push to remote after initialization (default: %(default)s)")
    choose_parser.add_argument("--choose_only", action="store_true", help="Only make the choices and plots, do not initialize next epoch (default: %(default)s)")
    choose_parser.add_argument("--reload_fval", action="store_true", help="Recalculate fval even if fval.npy exists (default: %(default)s)")
    choose_parser.set_defaults(func=choose)

    # newepoch command
    epochstarter_parser = subparsers.add_parser("newepoch", help="Shorthand for \"choose --pull --push\".")
    epochstarter_parser.set_defaults(func=choose, push=True, pull=True, choose_only=False)
    epochstarter_parser.add_argument("--reload_fval", action="store_true", help="Recalculate fval even if fval.npy exists (default: %(default)s)")

    # Push and pull commands
    push_parser = subparsers.add_parser("push", help="rsync from local to remote")
    push_parser.set_defaults(func=pushpull,config_func=import_cfg,pull=False)
    pull_parser = subparsers.add_parser("pull", help="rsync from remote to local")
    pull_parser.set_defaults(func=pushpull,config_func=import_cfg,pull=True)


    arguments = parser.parse_args()
    assert arguments.config.exists(), "Config file does not exist at %s."%arguments.config

    arguments.cfg = arguments.config_func(arguments.config)

    return arguments
