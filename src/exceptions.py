# -*- coding: utf-8 -*-

import argparse
import sys
from typing import Callable


class FSTException(Exception):
    """
    A base exception differnetiate between fst-internal exception,
    to be caught nicely, and some other code's exceptions, which should not be caught.
    """
    code = 1
    pass


def handle_errors(func: Callable[[argparse.Namespace], None]) -> Callable[[argparse.Namespace], None]:
    def wrapped_func(args: argparse.Namespace) -> None:
        try:
            return func(args)
        except FSTException as e:
            print("--------------------------\n", file=sys.stderr)
            print(e.__class__.__name__, file=sys.stderr)
            print("", file=sys.stderr)
            print(e,  file=sys.stderr)
            print("\n--------------------------", file=sys.stderr)
            sys.exit(e.code)
    return wrapped_func


class DeprecatedUsageWarning(UserWarning):
    """
    A warning to let the user know, that what they are doing will be deprecated.
    """
    pass


class NoEpochsFoundError(FSTException, FileNotFoundError):
    """
    An exception raised, when no Epoch data is found, even though the command requires it.
    """
    pass


class RepExistsError(FSTException, FileExistsError):
    """
    An exception raised, when a rep folder is already present.
    """
    pass


class NonzeroReturnError(FSTException, RuntimeError):
    """
    An exception raised, when a subprocess returns a nonzero return code.
    """
    pass


class NoConfigError(FSTException, FileNotFoundError):
    """
    An exception raised, when the config file does not exists.
    """
    pass


class NoSbatchLaunchError(FSTException, FileNotFoundError):
    """
    An exception raised, when the config file does not exists.
    """
    pass
