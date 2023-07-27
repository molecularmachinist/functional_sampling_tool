#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from MDAnalysis.analysis import align
from MDAnalysis.lib.mdamath import triclinic_vectors
import warnings

from numpy.typing import NDArray

from . import _ctransformations

# Type hints
from typing import Union, Optional, List
from MDAnalysis.core.groups import AtomGroup, Atom
from MDAnalysis.coordinates.base import Timestep


"""
This module includes functions to make on the fly transformations for MDAnalysis trajectories
"""


def make_whole(ts: Timestep, bonds: NDArray[np.int_], sel: NDArray[np.int_]):
    box = ts.triclinic_dimensions
    ts.positions[sel] = _ctransformations.make_whole(
        ts.positions[sel], bonds, box)

    return ts


class Unwrapper:
    """ Make molecules in ag whole over the pbc. only considers
        continuosly bonded molecules, so if molecules are in many parts,
        the nonbonded parts will be ignored and can be broken.
            Starters is a list (or atom group) of atoms to use as starting points.
        They will be put into the box if they are outside and each consecutive
        bonded atom will be moved by a box vector if it is more than half the
        length of the box.

        parameters:
            ag:        Atom group of molecules to make whole
            starters:  Atom goup (or list) of atoms that are guaranteed to stay
                       in box. If multiple atoms are part of same molecule, only
                       first is guaranteed.
            initsetup: If True, setup is run when initializing object.

        returns:
            transformation function
    """

    def __init__(self,
                 ag: AtomGroup,
                 starters: Union[List[Atom], AtomGroup] = [],
                 initsetup: bool = False):
        self.starters = starters
        self.ag = ag
        self.sel = self.ag.indices
        self.__setup_run = False
        if (initsetup):
            self._setup()

    def _setup(self) -> None:
        if (self.__setup_run):
            warnings.warn("Unwrapper setup being run after __setup_run is already. "
                          "Continuing without new setup.", RuntimeWarning)
        self.bonds = _ctransformations.traverse_mol(self.ag.indices,
                                                    self.ag.universe.bonds.to_indices(),
                                                    np.array([s.index for s in self.starters], dtype=int))
        self.__setup_run = True
        # No need to remember atom group or starters after setup
        del self.ag
        del self.starters

    def __call__(self, ts: Timestep) -> Timestep:
        if (not self.__setup_run):
            self._setup()
        return make_whole(ts, self.bonds, self.sel)


class MolWrapper:
    """
    Put centre of mass of molecules in selection to box

    parameters:
        ag:        Atom group of molecules to put in box
        initsetup: If True, setup is run when initializing object.

    returns:
        transformation function
    """

    def __init__(self, ag: AtomGroup, initsetup: bool = False):
        self.ag = ag
        self.selection = self.ag.indices
        self.__setup_run = False
        if (initsetup):
            self._setup()

    def _setup(self) -> None:
        if (self.__setup_run):
            warnings.warn("Unwrapper setup being run after __setup_run is already True."
                          "Continuing without new setup.", RuntimeWarning)

        self.mols, self.nmols = _ctransformations.find_frags(
            self.selection,
            self.ag.universe.bonds.to_indices()
        )
        self.weights = self.ag.masses.astype(self.ag.positions.dtype)
        self.__setup_run = True
        # No need to remember atom group after setup
        del self.ag

    def __call__(self, ts: Timestep) -> Timestep:
        if (not self.__setup_run):
            self._setup()

        box = ts.triclinic_dimensions
        ts.positions[self.selection] = _ctransformations.wrap_mols(
            ts.positions[self.selection],
            self.weights,
            self.mols,
            self.nmols,
            box
        )
        return ts


class Superpos:
    """
    Superposition for optimal mass weighted rmsd

    parameters:
        ag:            Atom group of atoms to fit, from the reference universe
        centre:        Boolean of whether to centre the selection
        superposition: Boolean of whether to centre the selection and fit rotationally
        subselection:  The atom group to move and/or rotate, None to use ag. [default: None]

    returns:
        transformation function
    """

    def __init__(self,
                 ag: AtomGroup,
                 centre: bool,
                 superposition: bool,
                 subselection: Optional[AtomGroup] = None):
        self.seli = ag.indices.copy()
        self.ref = ag.positions.copy()
        self.ref_com = ag.center_of_mass()
        self.w = ag.masses.copy()
        self.totw = self.w.sum()
        if (subselection is None):
            self.subsel = self.seli
        else:
            self.subsel = subselection.indices
        if (superposition):
            self.func = self.superpos
        elif (centre):
            self.func = self.centre
        else:
            self.func = self.nothing

    def __call__(self, ts: Timestep) -> Timestep:
        return self.func(ts)

    def superpos(self, ts: Timestep) -> Timestep:
        sel = ts.positions[self.seli].copy()
        sel_com = np.sum(self.w*sel.T, axis=-1)/self.totw
        sel -= sel_com
        R, rmsd = align.rotation_matrix(
            sel, self.ref-self.ref_com, weights=self.w)
        sel = sel @ R.T
        ts.positions[self.subsel] = (
            (ts.positions[self.subsel]-sel_com) @ R.T)+self.ref_com
        return ts

    def centre(self, ts: Timestep) -> Timestep:
        sel = ts.positions[self.seli]
        sel_com = np.sum(self.w*sel.T, axis=-1)/self.totw
        ts.positions[self.subsel] += self.ref_com-sel_com
        return ts

    def nothing(self, ts: Timestep) -> Timestep:
        return ts


class Precenter:
    """
    Center a single atom before putting molecules back in box 

    parameters:
        ag:             Atom group of atoms to move, from the reference universe
        centre_atom:    The Atom to centre. If None the closest atom to box centre
                        is used. [default: None]
        subselection:   The atom group to move, None to use ag. [default: None]

    returns:
        transformation function
    """

    def __init__(self,
                 ag: AtomGroup,
                 centre_atom: Optional[Atom] = None,
                 subselection: Optional[AtomGroup] = None):
        self.seli = ag.indices.copy()
        if (centre_atom is None):
            box = triclinic_vectors(ag.universe.dimensions)
            box_center = (box/2).sum(axis=0)
            dists = np.linalg.norm(ag.positions-box_center, axis=-1)
            self.centre_atom = ag[dists.argmin()].index
        else:
            self.centre_atom = centre_atom.index

        if (subselection is None):
            self.sel = ag.indices
        else:
            self.sel = subselection.indices

    def __call__(self, ts: Timestep) -> Timestep:
        box = ts.triclinic_dimensions
        box_center = (box/2).sum(axis=0)
        diff = box_center-ts.positions[self.centre_atom]
        ts.positions[self.sel] = ts.positions[self.sel]+diff
        return ts
