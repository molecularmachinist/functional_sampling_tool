#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import pathlib
import shutil
import sys

from . import __version__

from .exceptions import handle_errors

# For python 3.10 (or higher) import from standard lib, older releases use importlib_resources.
if (sys.version_info >= (3, 10)):
    from importlib.resources import files as import_files
else:
    from importlib_resources import files as import_files


@handle_errors
def init(args: argparse.Namespace) -> None:
    from .epoch_starting import start_epoch
    from .utils import rsync_up
    from .inout import config_validation
    config_validation(args.cfg)
    print("Initializing first epoch")
    start_epoch(1, args.cfg, numproc=args.num_proc)
    if (args.push):
        rsync_up(args.cfg)


@handle_errors
def pushpull(args: argparse.Namespace) -> None:
    from .utils import rsync_up, rsync_down
    if (args.pull):
        rsync_down(args.cfg, args.sync_all)
    else:
        rsync_up(args.cfg, args.sync_all)


@handle_errors
def choose(args: argparse.Namespace) -> None:
    from .clustering import ClusterChooser
    from .epoch_starting import start_epoch
    from .utils import rsync_up, rsync_down
    if (args.pull):
        rsync_down(args.cfg)
    chooser = ClusterChooser.fromReadData(
        args.cfg, args.reload_fval)
    chooser.plot_hist()
    val, epc, rep, frm = chooser.make_choices()
    if (not args.choose_only):
        start_epoch(
            chooser.nextepoch,
            args.cfg,
            val,
            epc,
            rep,
            frm,
            numproc=args.num_proc
        )
    else:
        chooser.print_choices(val, epc, rep, frm)
    if (args.push):
        rsync_up(args.cfg)


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
    from .inout import clean_latest_epoch
    print("Running cleanup")
    clean_latest_epoch(args.force)


@handle_errors
def config_importer(config: pathlib.Path):
    from .inout import import_cfg
    return import_cfg(config)


@handle_errors
def config_loader(config: pathlib.Path):
    from .inout import load_options
    return load_options(config)


def argP() -> argparse.Namespace:
    """
    Define command line arguments. For each
    """
    # This should be importable without importing the whole lot,
    # but importing it at the very beginning makes for a circular
    # import.
    from .analysis import parsers as analysis_parsers

    description = f"Functional Sampling Tool: A tool for functional sampling."
    epilog = "For more detailed explanations, see the README.md"
    parser = argparse.ArgumentParser(
        prog="fst", description=description, epilog=epilog)

    # Option to print program version
    parser.add_argument("-V", "--version", action="version",
                        version=f"{__package__} {__version__}")
    # The only global option
    parser.add_argument("-c", "--config", metavar="<name>.py", default=pathlib.Path("config.py"),
                        type=pathlib.Path, help="Path to config file (default: %(default)s)")
    # Global default for config_func
    parser.set_defaults(config_func=config_loader)

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
    init_parser.add_argument("--num-proc", type=int, default=None,
                             help="Maximum number of parallel processes to use for grompping. By default uses at most half the number of logical "
                             "CPUs. The actual number of threads is limited by the number of repetitions to start.")
    # The initializer has a different config func to not try to load the struct.
    init_parser.set_defaults(func=init, config_func=config_importer)

    # The chooser command
    choose_parser = subparsers.add_parser(
        "choose", help="Choose frames for next epoch and initialize it.")
    choose_parser.add_argument("--pull", action="store_true",
                               help="pull from remote before choosing (default: %(default)s)")
    choose_parser.add_argument("--push", action="store_true",
                               help="push to remote after choosing (default: %(default)s)")
    choose_parser.add_argument("--choose-only", action="store_true",
                               help="Only make the choices and plots, do not initialize next epoch (default: %(default)s)")
    choose_parser.add_argument("--reload-fval", action="store_true",
                               help="Reload data from mdrun.xtc even if fval_data.npz exists and everything matches (default: %(default)s)")
    choose_parser.add_argument("--num-proc", type=int, default=None,
                               help="Maximum number of parallel threads to use for grompping. By default uses at most half the number of logical "
                               "CPUs. The actual number of threads is limited by the number of repetitions to start.")
    choose_parser.set_defaults(func=choose)

    # newepoch command
    epochstarter_parser = subparsers.add_parser(
        "newepoch", help="Shorthand for \"choose --pull --push\".")
    epochstarter_parser.set_defaults(
        func=choose, push=True, pull=True, choose_only=False)
    epochstarter_parser.add_argument("--reload-fval", action="store_true",
                                     help="Reload data from mdrun.xtc even if fval_data.npz exists and everything matches (default: %(default)s)")
    epochstarter_parser.add_argument("--num-proc", type=int, default=None,
                                     help="Maximum number of parallel threads to use for grompping. By default uses at most half the number of logical "
                                     "CPUs. The actual number of threads is limited by the number of repetitions to start.")

    # Push and pull commands
    push_parser = subparsers.add_parser(
        "push", help="rsync from local to remote")
    push_parser.add_argument("--sync-all", action="store_true",
                             help="Push all epochs, not just latest (default: %(default)s)")
    push_parser.set_defaults(
        func=pushpull, config_func=config_importer, pull=False)
    pull_parser = subparsers.add_parser(
        "pull", help="rsync from remote to local")
    pull_parser.add_argument("--sync-all", action="store_true",
                             help="Pull all epochs, not just latest (default: %(default)s)")
    pull_parser.set_defaults(
        func=pushpull, config_func=config_importer, pull=True)

    # clean commands
    clean_parser = subparsers.add_parser("clean", help="clean latest epoch")
    clean_parser.set_defaults(func=clean, config_func=(lambda cfgpath: None))
    clean_parser.add_argument("--force", action="store_true", help="Force the removal of the latest epoch dir, "
                              "even if a mdrun.xtc file can be found in one of the repetition folders. Note, this "
                              "can easily lead to loss of data")

    # Template command
    templ_parser = subparsers.add_parser(
        "make-templates", help="Copy default config.py and sbatch_template.sh files")
    templ_parser.set_defaults(
        func=copy_templates, config_func=(lambda cfgpath: None))
    templ_parser.add_argument("--config-out", metavar="<name>.py", default=pathlib.Path("config.py"),
                              help="Filename of produced config file (default: %(default)s)", type=pathlib.Path)
    templ_parser.add_argument("--sbatch-out", metavar="<name>.sh", default=pathlib.Path(
        "sbatch_launch.sh"), help="Filename of produced sbatch file (default: %(default)s)", type=pathlib.Path)
    templ_parser.add_argument("--no-config",  action="store_true",
                              help="Do not produce config file (default: %(default)s)")
    templ_parser.add_argument("--no-sbatch",  action="store_true",
                              help="Do not produce sbatch file (default: %(default)s)")

    # Analysis command
    analysis_parser = subparsers.add_parser(
        "analysis", help="Utilities to help with analysis")
    analysis_parsers.analysis_subparser(analysis_parser)

    arguments = parser.parse_args()
    arguments.cfg = arguments.config_func(arguments.config)

    return arguments


def _run_tool() -> None:
    args = argP()
    args.func(args)
