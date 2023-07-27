import argparse
import pathlib

from ..exceptions import handle_errors
from ..cmd_line import config_importer


@handle_errors
def ancestry_graphing(args: argparse.Namespace):
    from .ancestry import ancestry
    ancestry(args)


@handle_errors
def extraction(args: argparse.Namespace):
    from .extract import extract
    extract(args)


def extract_subparser(subparsers: "argparse._SubParsersAction[argparse.ArgumentParser]") -> None:
    """
    Add the extraction options to the provided subparser
    """

    extract_description = "Only repetitions for which the data " \
        "npz-archive is present can be extracted. Run the choose command with " \
        "the --choose-only option to calculate the functions and make the archives " \
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
    extr_parser.add_argument("--precenter",   action="store_true",
                             help="Precenter the --sel-unwrap by shifting it such that --precenter-atom is in the centre of the box  (default: No)")
    extr_parser.add_argument("--precenter-atom", metavar="str",
                             help="The selection to unwrap. Same syntax as --selection, but should result in a single atom being selected. "
                             "By default None, which results in using whichever atom in --sel-unwrap is closest to box center initially.",
                             default=None)
    extr_parser.add_argument("--unwrap",   action="store_true",
                             help="Unwrap molecules (default: No)")
    extr_parser.add_argument("--sel-unwrap", metavar="str",
                             help="The selection to unwrap. Same syntax as --selection, and by default uses same selection.",
                             default=None)
    extr_parser.add_argument("--unwrap-starters", metavar="str",
                             help="The selection to use as starters in unwrapping. \"unwrap_starters\" to use from config, by default use none",
                             default=None)
    extr_parser.add_argument("--wrap",     action="store_true",
                             help="Put centre of mass of molecules back in box. "
                                  "Uses the --sel-unwrap. (default: No)")
    extr_parser.add_argument("--superpos-trans", action="store_true",
                             help="Superposition the --sel-superpos selection onto the initial structure, only translating (default: No)")
    extr_parser.add_argument("--superposition",     action="store_true",
                             help="Superposition the --sel-superpos selection onto the initial structure, translating and rotating (default: No)")
    extr_parser.add_argument("--sel-superpos", metavar="str",
                             help="The selection to superposition. Same syntax as --selection, and by default uses same selection.",
                             default=None)
    # Use just config import as config_func, so the selections and such are not loaded
    extr_parser.set_defaults(func=extraction,
                             config_func=config_importer)


def ancestry_subparser(subparsers: "argparse._SubParsersAction[argparse.ArgumentParser]") -> None:
    """
    Add the ancestry graph options to the provided subparser
    """

    ancestry_description = "Make a graph of which simulation each new simulation is spawned from."

    # extractor command
    anc_parser = subparsers.add_parser("ancestry", help="Make an ancestry graph of the simulations.",
                                       description=ancestry_description)
    anc_parser.add_argument("-o", "--output", metavar="<name>.pdf", dest="output",
                            help="Output file for the graph. The parent directory will be created if it does not yet exist (default: %(default)s)",
                            default=pathlib.Path("analysis/figs/ancestry.pdf"),
                            type=pathlib.Path)
    anc_parser.add_argument("--include-start",     action="store_true",
                            help="Add the starting structure(s) also to the graph. (default: do not include)")
    anc_parser.add_argument("--doignore",     action="store_true", dest="doignore",
                            help="Ignore the epochs and repetitions defined in the config. "
                            "Note that this can result in problems if any of the ignored simulations "
                            "is a parent to another one. (default: do not ignore)")

    # Use just config import as config_func, so the selections and such are not loaded
    anc_parser.set_defaults(func=ancestry_graphing,
                            config_func=config_importer)


def analysis_subparser(parser: argparse.ArgumentParser):
    """
    Add the analysis options to the provided parser
    """

    subparsers = parser.add_subparsers(
        description="Specify which analysis function to run")
    subparsers.required = True
    subparsers.dest = 'subcommand'

    extract_subparser(subparsers)
    ancestry_subparser(subparsers)
