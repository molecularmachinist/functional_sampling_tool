# -*- coding: utf-8 -*-

import argparse
import sys
from typing import Callable, Any


class FSTException(Exception):
    """
    A base exception to differentiate between fst-internal exception (to be caught
    nicely) and some other code's exceptions (which should not be caught).
    """

    def __init__(self, *args: Any, code=1) -> None:
        super().__init__(*args)
        self.code = code


def handle_errors(func: Callable[[argparse.Namespace], None]) -> Callable[[argparse.Namespace], None]:
    """
    A function decorator that will handle FSTException and KeyboardInterrupt nicely
    while not catching any other Exceptions.

    parameters:
        func: Callable(Namespace) -> None
            The function to wrap with the error handling.

    return:
        wrapped_func: Callable(Namespace) -> None
            The wrapped function.
    """
    def wrapped_func(args: argparse.Namespace) -> None:
        try:
            return func(args)
        except FSTException as e:
            print("\n------------------------------------------------------\n",
                  file=sys.stderr)
            print(e.__class__.__name__, file=sys.stderr)
            print("", file=sys.stderr)
            print(e,  file=sys.stderr)
            print("\n----------------------------------------------------",
                  file=sys.stderr)
            sys.exit(e.code)
        except KeyboardInterrupt as e:
            print("\n------------------------------------------------------\n",
                  file=sys.stderr)
            print("KeyboardInterrupt")
            print("", file=sys.stderr)
            print("Stopping due to user interrupt",  file=sys.stderr)
            print("\n----------------------------------------------------",
                  file=sys.stderr)
            sys.exit(130)
    return wrapped_func


class DeprecatedUsageWarning(UserWarning):
    """
    A warning to let the user know, that what they are doing will be deprecated.
    """
    pass


class FunctionDimensionError(FSTException, ValueError):
    """
    An exception raised, when the function value is wrong length or dimensionality.
    """
    pass


class NotEnoughDataError(FSTException, ValueError):
    """
    An exception raised, there is not enough data to make histograms or choose new epoch.
    """
    pass


class WrongSelectionSizeError(FSTException, ValueError):
    """
    An exception raised, when a selection should be a certain size, but is not.
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


class ExternalProgramMissingError(FSTException, FileNotFoundError):
    """
    An exception raised, when a subprocess does not find the command.
    """

    def __init__(self, *args: Any) -> None:
        super().__init__(*args, code=127)


class NonzeroReturnError(FSTException, RuntimeError):
    """
    An exception raised, when a subprocess returns a nonzero return code.
    """
    pass


class RequiredFileMissingError(FSTException, FileNotFoundError):
    """
    An exception raised, when a required file does not exists.
    """
    pass


class NoConfigError(RequiredFileMissingError):
    """
    An exception raised, when the config file does not exists.
    """
    pass


class NoSbatchLaunchError(RequiredFileMissingError):
    """
    An exception raised, when the sbatch file does not exists.
    """
    pass
