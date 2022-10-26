#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import pathlib
import shutil
import sys

from . import inout
from . import epoch_starting
from . import clustering
from . import utils
from . import analysis
from .exceptions import handle_errors

# For python 3.10 (or higher) import from standard lib, older releases use importlib_resources.
if (sys.version_info >= (3, 10)):
    from importlib.resources import files as import_files
else:
    from importlib_resources import files as import_files


@handle_errors
def init(args: argparse.Namespace) -> None:
    print("Initializing first epoch")
    epoch_starting.start_epoch(1, args.cfg)
    if (args.push):
        utils.rsync_up(args.cfg)


@handle_errors
def pushpull(args: argparse.Namespace) -> None:
    if (args.pull):
        utils.rsync_down(args.cfg)
    else:
        utils.rsync_up(args.cfg)


@handle_errors
def choose(args: argparse.Namespace) -> None:
    if (args.pull):
        utils.rsync_down(args.cfg)
    chooser = clustering.ClusterChooser.fromReadData(
        args.cfg, args.reload_fval)
    chooser.plot_hist()
    val, epc, rep, frm = chooser.make_choices()
    if (not args.choose_only):
        epoch_starting.start_epoch(
            chooser.nextepoch, args.cfg, val, epc, rep, frm)
    else:
        chooser.print_choices(val, epc, rep, frm)
    if (args.push):
        utils.rsync_up(args.cfg)


@handle_errors
def copy_templates(args: argparse.Namespace) -> None:
    if (not args.no_config):
        args.config_out.parent.mkdir(parents=True, exist_ok=True)
        fin = import_files("%s.templates" %
                           __package__).joinpath("config_example.py.txt")
        print("Making", args.config_out)
        shutil.copyfile(fin, args.config_out)
    if (not args.no_sbatch):
        args.sbatch_out.parent.mkdir(parents=True, exist_ok=True)
        fin = import_files("%s.templates" %
                           __package__).joinpath("sbatch_template.sh")
        print("Making", args.sbatch_out)
        shutil.copyfile(fin, args.sbatch_out)


@handle_errors
def clean(args: argparse.Namespace) -> None:
    print("Running cleanup")
    inout.clean_latest_epoch(args.force)


def argP() -> argparse.Namespace:
    """
    Define command line arguments. For each
    """
    description = "Functional Sampling Tool: A tool for functional sampling."
    epilog = "For more detailed explanations, see the README.md"
    parser = argparse.ArgumentParser(
        prog="fst", description=description, epilog=epilog)

    # The only global option
    parser.add_argument("-c", "--config", metavar="<name>.py", default=pathlib.Path("config.py"),
                        type=pathlib.Path, help="Path to config file (default: %(default)s)")
    # Global default for config_func
    parser.set_defaults(config_func=inout.load_options)

    # A subparser object to make subcommands. A subcommand is required.
    subparsers = parser.add_subparsers(
        description="Specify what the program should do. Only one command should be given per run. 'fst <cmd> -h' to get command specific option.")
    subparsers.required = True
    subparsers.dest = 'command'

    # First epoch initializer command
    init_parser = subparsers.add_parser(
        "init", help="Initialize the first epoch.")
    init_parser.add_argument("--push", action="store_true",
                             help="push to remote after initialization (default: %(default)s)")
    # The initializer has a different config func to not try to load the struct.
    init_parser.set_defaults(func=init, config_func=inout.import_cfg)

    # The chooser command
    choose_parser = subparsers.add_parser(
        "choose", help="Choose frames for next epoch and initialize it.")
    choose_parser.add_argument("--pull", action="store_true",
                               help="pull from remote before choosing (default: %(default)s)")
    choose_parser.add_argument("--push", action="store_true",
                               help="push to remote after choosing (default: %(default)s)")
    choose_parser.add_argument("--choose_only", action="store_true",
                               help="Only make the choices and plots, do not initialize next epoch (default: %(default)s)")
    choose_parser.add_argument("--reload_fval", action="store_true",
                               help="Reload data from mdrun.xtc even if fval_data.npz exists and everything matches (default: %(default)s)")
    choose_parser.set_defaults(func=choose)

    # newepoch command
    epochstarter_parser = subparsers.add_parser(
        "newepoch", help="Shorthand for \"choose --pull --push\".")
    epochstarter_parser.set_defaults(
        func=choose, push=True, pull=True, choose_only=False)
    epochstarter_parser.add_argument("--reload_fval", action="store_true",
                                     help="Reload data from mdrun.xtc even if fval_data.npz exists and everything matches (default: %(default)s)")

    # Push and pull commands
    push_parser = subparsers.add_parser(
        "push", help="rsync from local to remote")
    push_parser.set_defaults(
        func=pushpull, config_func=inout.import_cfg, pull=False)
    pull_parser = subparsers.add_parser(
        "pull", help="rsync from remote to local")
    pull_parser.set_defaults(
        func=pushpull, config_func=inout.import_cfg, pull=True)

    # clean commands
    push_parser = subparsers.add_parser("clean", help="clean latest epoch")
    push_parser.set_defaults(func=clean, config_func=(lambda cfgpath: None))
    push_parser.add_argument("--force", action="store_true", help="Force the removal of the latest epoch dir, "
                             "even if a mdrun.xtc file can be found in one of the repetition folders. Note, this "
                             "can easily lead to oss of data")

    # Template command
    templ_parser = subparsers.add_parser(
        "make_templates", help="Copy default config.py and sbatch_template.sh files")
    templ_parser.set_defaults(
        func=copy_templates, config_func=(lambda cfgpath: None))
    templ_parser.add_argument("--config_out", metavar="<name>.py", default=pathlib.Path("config.py"),
                              help="Filename of produced config file (default: %(default)s)", type=pathlib.Path)
    templ_parser.add_argument("--sbatch_out", metavar="<name>.sh", default=pathlib.Path(
        "sbatch_launch.sh"), help="Filename of produced sbatch file (default: %(default)s)", type=pathlib.Path)
    templ_parser.add_argument("--no_config",  action="store_true",
                              help="Do not produce config file (default: %(default)s)")
    templ_parser.add_argument("--no_sbatch",  action="store_true",
                              help="Do not produce sbatch file (default: %(default)s)")

    # Analysis command
    analysis_parser = subparsers.add_parser(
        "analysis", help="Utilities to help with analysis")
    analysis.analysis_subparser(analysis_parser)

    arguments = parser.parse_args()
    arguments.cfg = handle_errors(arguments.config_func)(arguments.config)

    return arguments


def _run_tool() -> None:
    args = argP()
    args.func(args)
