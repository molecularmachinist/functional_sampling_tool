#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import mdtraj


def function_val(positions):
    """
      Write here your analysis function. Positions will be
        numpy array of shape (n,m,3) for n frames of m atoms.
        Note that m is the number of atoms in the selection,
        not the whole trajectory.
      The function should return a numpy array of shape (n).
    """
    m = positions.shape[-2]
    # 4 residues match, so we separate them
    m_res = m//4
    # selections for each res
    ind1 = np.array(list(range(m_res)))
    ind2 = m_res+ind1
    ind3 = m_res+ind2
    ind4 = m_res+ind3

    # distances between average position of res1 to average position of res3
    dist1 = np.linalg.norm(np.mean(positions[:,ind1,:], axis=-2)-np.mean(positions[:,ind3,:], axis=-2), axis=-1)
    # same for res 2 and 4
    dist2 = np.linalg.norm(np.mean(positions[:,ind2,:], axis=-2)-np.mean(positions[:,ind4,:], axis=-2), axis=-1)

    return np.minimum(dist1, dist2)

def function_val2(positions):
    """
      Write here your analysis function. Positions will be
        numpy array of shape (n,m,3) for n frames of m atoms.
        Note that m is the number of atoms in the selection,
        not the whole trajectory.
      The function should return a numpy array of shape (n).
    """
    m = positions.shape[-2]
    # 4 residues match, so we separate them
    m_res = m//4
    # selections for each res
    ind1 = np.array(list(range(m_res)))
    ind2 = m_res+ind1
    ind3 = m_res+ind2
    ind4 = m_res+ind3

    # Vectors from res1 to res3 and res2 to res4
    v13 = np.mean(positions[:,ind1,:], axis=-2)-np.mean(positions[:,ind3,:], axis=-2)
    v24 = np.mean(positions[:,ind2,:], axis=-2)-np.mean(positions[:,ind4,:], axis=-2)

    # Area is the length of the cross product
    A = np.linalg.norm(np.cross(v13,v24), axis=-1)

    return A


def calcval(fname,function_val=function_val):

    select_str = "protein and residue 617 and resname THR"
    print("Loading structure and making selection")
    struct = mdtraj.load(fname)
    sel = struct.topology.select(select_str)

    return function_val(struct.xyz[:,sel,:])

if __name__ == '__main__':
    o_struct = "5weo.pdb"
    c_struct = "closed/5wek.pdb"
    print(f"fval  for open:   {calcval(o_struct,function_val)}")
    print(f"fval2 for open:   {calcval(o_struct,function_val2)}")
    print(f"fval  for closed: {calcval(c_struct,function_val)}")
    print(f"fval2 for closed: {calcval(c_struct,function_val2)}")
