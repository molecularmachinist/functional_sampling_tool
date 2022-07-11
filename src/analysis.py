import argparse

from . import inout

def analysis_subparser(parser):
    """
    Add the analysis options to the provided parser
    """
    subparsers = parser.add_subparsers(description = "Specify which analysis function to run")
    subparsers.required = True
    subparsers.dest = 'subcommand'


    # First epoch initializer command
    extr_parser = subparsers.add_parser("extract", help="Extract frames.")
    extr_parser.add_argument("--push", action="store_true", help="push to remote after initialization (default: %(default)s)")
    # The initializer has a different config func to not try to load the struct.
    extr_parser.set_defaults(func=extract, config_func=inout.import_cfg)


def extract(args):
    print("I'm extracting")