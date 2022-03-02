#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from MDAnalysis.analysis import align
import warnings
from numba import njit
from numba.typed import List,Dict

"""
This module includes functions to make on the fly transformations for MDAnalysis trajectories
"""


def find_frags(ag):
    """ Traverse all molecules that atoms in ag are part of.
        Parameters:
            ag:       Atom group of atoms to check
            starters: List of atoms or atom group to use as starting points. (optional)
        Returns:
            List of indices for each molecule, only including indices in ag.
    """
    # A dict to map the indices between system and selection
    agi = set(ag.indices)
    mols = []

    # We iterate until there are no more atoms to search
    while(agi):
        # A set to hold info on whether we have "visited" atom
        done = set()
        start_i = agi.pop()
        start = ag.universe.atoms[start_i]

        # Make a stack to hold atoms to visit and add start atom
        stack = [start]
        mol = [start.index]
        done.add(start.index)

        # Continue until stack is empty
        while(stack):
            # Pop from top of stack
            curr = stack.pop()
            # Iterate atoms bonded to popped atom
            for a in curr.bonded_atoms:
                # Ignore if atom visited or if it is not in selection
                if(a.index in done):
                    continue

                # Add atom to list
                if(a.index in agi):
                    mol.append(a.index)
                    agi.remove(a.index)
                # Add atom on top of stack and mark it as visited
                stack.append(a)
                done.add(a.index)


        mols.append(np.unique(mol))

    return mols


def traverse_mol(ag, starters=[]):
    """ Traverse all graphs of bonded atoms in ag, using atoms in starters (or the
        smallest indexes) as starting points for a DFS.
        Parameters:
            ag:       Atom group of atoms to check
            starters: List of atoms or atom group to use as starting points. (optional)
        Returns:
            Indexes for each bonds of the bonded atoms within the selection
            (i.e. indexes within selection, not within whole universe). int[n,2].
    """
    # make sure this is a list, not ag and flip it to pop from the "start"
    starters=list(starters)[::-1]
    # A dict to map the indices between system and selection
    agi_indexes = {}
    for i,a in enumerate(ag):
        agi_indexes[a.index] = i

    # A list to hold the bonds
    bonds = []

    # We iterate until there are no more atoms to search
    while(agi_indexes):
        # A set to hold info on whether we have "visited" atom
        done = set()
        # Decide next starting atom
        start = None
        # Go thorugh starter atoms until it is empty or we choose one
        while(start is None and starters):
            a = starters.pop()
            if a.index in agi_indexes:
                start = a

        # If we didn't choose one, taek whichever has smallest index
        if(start is None):
            start_i = min(agi_indexes.keys())
            start = ag.universe.atoms[start_i]

        # Make a stack to hold atoms to visit and add start atom
        stack = [start]
        done.add(start.index)
        # Continue until stack is empty
        while(stack):
            # Pop from top of stack
            curr = stack.pop()
            # Iterate atoms bonded to popped atom
            for a in curr.bonded_atoms:
                # Ignore if atom visited or if it is not in selection
                if(a.index in done or not a.index in agi_indexes):
                    continue

                # Add bond to list
                bonds.append((agi_indexes[curr.index], agi_indexes[a.index]))
                # Add atom on top of stack and mark it as visited
                stack.append(a)
                done.add(a.index)

        # Remove visited atoms from selction to keep track of whether we are done
        for a in done:
            agi_indexes.pop(a)

    if(not bonds):
        warnings.warn("No bonds found in unwrap selection")
        return np.zeros((0,2),dtype=int)
    return np.array(bonds)


@njit(["f8[:,:](f8[:,:],i8[:,:])",
      "f4[:,:](f4[:,:],i8[:,:])"])
def _iterate_bonds(pos, bonds):
    """
    Fixes molecules whole over pbc.
    Parameters:
        pos:   Positions in reciprocal space float[n,d].
        bonds: array of atom indices corresponding to bonds.
               The second one will be fixed. int[m,2].
    Returns:
        Fixed positions in reciprocal space. float[n,d].
    """
    # Iterate over bonds
    for atoms in bonds:
        a1 = atoms[0]
        a2 = atoms[1]
        # Difference vector
        diff_v = pos[a2]-pos[a1]
        # If distance is more than half box vector, translate by box vector (1.0)
        pos[a2] -= (np.abs(diff_v)>0.5)*np.sign(diff_v)
    # Return fixed coordinates
    return pos


def make_whole(ts, bonds, sel):
    """ Makes molecules in sel whole using bond info from bonds.
        Coordinates are first switched to reciprocal space, then passed to
        _iterate_bonds to actually make them whole. The fixed positions are
        switched back to normal space and set as the current positions for sel.
    """
    box = ts.triclinic_dimensions
    # Inverse box
    invbox = np.linalg.inv(box)
    # Transfer coordinates to unit cell and put in box
    unitpos = (ts.positions[sel] @ invbox) % 1
    # Fix box in reciprocal space
    unitpos = _iterate_bonds(unitpos, bonds)

    ts.positions[sel] = unitpos @ box

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
            ag:       Atom group of molecules to make whole
            starters: Atom goup (or list) of atoms that are guaranteed to stay
                      in box. If multiple atoms are part of same molecule, only
                      first is guaranteed.
        returns:
            transformation function
    """
    def __init__(self, ag, starters=[]):
        self.bonds = traverse_mol(ag, starters)
        self.sel = ag.indices

    def __call__(self, ts):
        return make_whole(ts, self.bonds, self.sel)




class MolWrapper:
    """
    Put centre of mass of molecules in selection to box
    parameters:
        ag:       Atom group of molecules to put in box
    returns:
        transformation function
    """
    def __init__(self,ag):
        self.selection = ag.indices
        self.mols = []
        self.weights = []
        for frag in find_frags(ag):
            self.mols.append(frag)
            self.weights.append(ag.universe.atoms[frag].masses)

        self.totw = np.array([np.sum(w) for w in self.weights])

    def __call__(self,ts):
        box = ts.triclinic_dimensions
        # Inverse box
        invbox = np.linalg.inv(box)
        for i,m in enumerate(self.mols):
            pos = np.sum(self.weights[i]*ts.positions[m].T,axis=-1)/self.totw[i]
            # Transfer coordinates to unit cell and put in box
            unitpos = (pos @ invbox) % 1
            newpos = unitpos @ box
            trans = newpos-pos
            ts.positions[m] += trans
        return ts



class Superpos:
    """
    Superposition for optimal mass weighted rmsd
    parameters:
        ag:            Atom group of atoms to fit, from the reference universe
        centre:        Boolean of whether to centre the selection
        superposition: Boolean of whether to centre the selection and fit rotationally
    returns:
        transformation function
    """
    def __init__(self,ag, centre, superposition):
        self.seli = ag.indices.copy()
        self.ref = ag.positions.copy()
        self.ref_com = ag.center_of_mass()
        self.w   = ag.masses.copy()
        self.totw = self.w.sum()
        if(superposition):
            self.func = self.superpos
        elif(centre):
            self.func = self.centre
        else:
            self.func = self.nothing

    def __call__(self, *args, **kwargs):
        return self.func(*args, **kwargs)

    def superpos(self,ts):
        sel = ts.positions[self.seli]
        sel_com = np.sum(self.w*sel.T, axis=-1)/self.totw
        sel -= sel_com
        R, rmsd = align.rotation_matrix(sel, self.ref-self.ref_com, weights=self.w)
        sel = sel @ R.T
        ts.positions[self.seli] = sel+self.ref_com
        return ts

    def centre(self,ts):
        sel = ts.positions[self.seli]
        sel_com = np.sum(self.w*sel.T, axis=-1)/self.totw
        ts.positions[self.seli] += self.ref_com-sel_com
        return ts

    def nothing(self,ts):
        return ts
