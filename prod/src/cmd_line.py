#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import mdtraj
import importlib


def import_cfg(cfgname):
    """
    Only import the config, does not load structs or do anything with it.
    """
    print("Loading config from %s"%cfgname)
    spec = importlib.util.spec_from_file_location("config", cfgname)
    cfg = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cfg)
    return cfg


def load_options(cfgname):
    """
    Imports the config, and loads the structs and selections
    """
    cfg    = import_cfg(cfgname)
    print("Loading structure")
    cfg.struct = mdtraj.load("initial/start.pdb")
    cfg.sel    = cfg.struct.topology.select(cfg.select_str)
    print("Selected %d atoms")
    cfg.startval = cfg.function_val(cfg.struct.xyz[:,cfg.sel,:])
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
