# cython: profile=False
import cython
import numpy as np
cimport numpy as np
import warnings

from libcpp.unordered_set cimport unordered_set as cunset
from libcpp.unordered_map cimport unordered_map as cunmap
from libcpp.stack cimport stack as cstack
from libcpp.vector cimport vector as cvector
from libc.math cimport rint

from cython          cimport floating 
from cython.operator cimport dereference

"""
This module includes functions to make on the fly transformations for MDAnalysis trajectories
"""

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef void cmatmul3D_mod(floating[:,:] A, floating[:,:] B, floating[:,:] out):
    """
    Perform matrix multiplication of (n,3) matrix A and (3,3) matrix B,
    taking modulo of one of the result. Same as (A@B)%1 for numpy arrays.
    Out must be (n,3) matrix and will be overwritten.
    """
    cdef int n,k,m,oint
    for n in range(A.shape[0]): #,schedule="static", nogil=True):
        for m in range(3):
            out[n,m] = 0
            for k in range(3):
                out[n,m] += A[n,k]*B[k,m]
            oint = int(out[n,m])
            if(oint!=0):
                out[n,m] -= oint
                if(oint<0):
                    out[n,m] += 1.0


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef void cmatmul3D(floating[:,:] A, floating[:,:] B, floating[:,:] out):
    """
    Perform matrix multiplication of (n,3) matrix A and (3,3) matrix B,
    same as A@B for numpy arrays. Out must be (n,3) matrix and will be overwritten.
    """
    cdef int n,k,m
    for n in range(A.shape[0]): #, nogil=True,schedule="static"):
        for m in range(3):
            out[n,m] = 0
            for k in range(3):
                out[n,m] += A[n,k]*B[k,m]



def _matmul3D(floating[:,:] A, floating[:,:] B, bint modulo=False):
    """
    Perform matrix multiplication of (n,3) matrix A and (3,3) matrix B,
    same as A@B for numpy arrays. Out must be (n,3) matrix and will be overwritten.
    """
    cdef cython.floating[:,:] out
    if(floating is double):
        out = np.empty(A.shape, dtype=np.float64)
    else: 
        out = np.empty(A.shape, dtype=np.float32)
    if(modulo):
        cmatmul3D_mod(A,B,out)
    else:
        cmatmul3D(A,B,out)

    return out


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def find_frags(np.intp_t[:] ag, int[:,:] bond_ind):
    """ Traverse all molecules that atoms in ag are part of.
        Parameters:
            ag:       Indices of atoms to check
            bonds:    int(n,2) array of atom indices for n bonds.

        Returns:
            mols:  Array of mol indices for each atom in ag.
            nmols: Number (int) of different molecules. Will satisfy np.arange(nmols)=np.unique(mols).
    """

    cdef int[:]                  mols = np.full(ag.shape[0],-1,dtype=np.int32)
    cdef cunset[int]             done
    cdef cunmap[int,int]         agi
    cdef cunmap[int,cunset[int]] bonds
    cdef int                     i,a,start,mol_ind=0
    cdef cstack[int]             stack

    # A dict to hold the bond info
    for i in range(bond_ind.shape[0]):
        bonds[bond_ind[i,0]].insert(bond_ind[i,1])
        bonds[bond_ind[i,1]].insert(bond_ind[i,0])

    # A dict to map the indices between system and selection
    for i in range(ag.shape[0]):
        agi[ag[i]] = i

    # We iterate until there are no more atoms to search
    while(not agi.empty()):
        # A set to hold info on whether we have "visited" atom
        done.clear()
        # Take any atom as the starting point for next search
        start = dereference(agi.begin()).first
        mols[agi[start]]=mol_ind
        #agi.erase(start)

        # Put start in stack
        stack.push(start)
        done.insert(start)

        # Continue until stack is empty
        while(not stack.empty()):
            # Pop from top of stack
            curr = stack.top()
            stack.pop()
            # Iterate atoms bonded to popped atom
            for a in bonds[curr]:
                # Ignore if atom visited
                if(done.count(a)):
                    continue

                # Add atom to current mol if it is in selection
                if(agi.count(a)):
                    mols[agi[a]]=mol_ind
            
                # Add atom on top of stack and mark it as visited
                stack.push(a)
                done.insert(a)

        # Remove visited atoms from selction to keep track of whether we are done
        for a in done:
            agi.erase(a)

        mol_ind+=1

    return mols, mol_ind


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def traverse_mol(np.intp_t[:] ag, int[:,:] bond_ind, np.intp_t[:] starters):
    """ Traverse all graphs of bonded atoms in ag, using atoms in starters (or the
        smallest indexes) as starting points for a DFS.
        Parameters:
            ag:       Indices of atoms to check
            bonds:    int(n,2) array of atom indices for n bonds.
            starters: Indices of atoms to use as starting points.
        Returns:
            Indexes for each bonds of the bonded atoms within the selection
            (i.e. indexes within selection, not within whole universe). int[n,2].
    """
    cdef int[:,:]                 bonds = np.empty((ag.shape[0],2),dtype=np.int32)
    cdef cunmap[int,int]          agi
    cdef cunset[int]              done
    cdef cunmap[int,cvector[int]] bond_map
    cdef int                      i, a, next_starter=0,next_first=0,next_bond=0
    cdef cstack[int]              stack

    # A dict to hold the bond info
    for i in range(bond_ind.shape[0]):
        bond_map[bond_ind[i,0]].push_back(bond_ind[i,1])
        bond_map[bond_ind[i,1]].push_back(bond_ind[i,0])
        
    # A dict to map the indices between system and selection
    for i in range(ag.shape[0]):
        agi[ag[i]] = i


    # We iterate until there are no more atoms to search
    while(not agi.empty()):
        # A set to hold info on whether we have "visited" atom
        done.clear()
        # Decide next starting atom
        start = -1
        # Go thorugh starter atoms until done or we choose one
        while(start < 0 and next_starter < starters.shape[0]):
            a = starters[next_starter]
            next_starter += 1
            if(agi.count(a)):
                start = a

        # If we didn't choose one take the next from ag still in agi
        while(start < 0):
            a = ag[next_first]
            next_first += 1
            if(agi.count(a)):
                start = a

        # Make a stack to hold atoms to visit and add start atom
        stack.push(start)
        done.insert(start)
        # Continue until stack is empty
        while(not stack.empty()):
            # Pop from top of stack
            curr = stack.top()
            stack.pop()
            # Iterate atoms bonded to popped atom
            for a in bond_map[curr]:
                # Ignore if atom visited or if it is not in selection
                if(done.count(a) or not agi.count(a)):
                    continue

                # Add bond to list
                bonds[next_bond,0] = agi[curr]
                bonds[next_bond,1] = agi[a]
                next_bond += 1
                # Add atom on top of stack and mark it as visited
                stack.push(a)
                done.insert(a)

        # Remove visited atoms from selction to keep track of whether we are done
        for a in done:
            agi.erase(a)

    if(next_bond==0):
        warnings.warn("No bonds found in unwrap selection")
    return bonds[:next_bond]


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef void iterate_bonds(floating[:,:] pos, int[:,:] bonds):
    """
    Fixes molecules whole over pbc.
    Parameters:
        pos:   Positions in reciprocal space float[n,d].
        bonds: array of atom indices corresponding to bonds.
               The second one will be fixed. int[m,2].
    Returns:
        Fixed positions in reciprocal space. float[n,d].
    """
    cdef int    i,j, a1, a2
    cdef floating diff
    # Iterate over bonds
    for i in range(bonds.shape[0]):
        a1 = bonds[i,0]
        a2 = bonds[i,1]
        for j in range(pos.shape[1]):
            # Difference vector
            diff = pos[a2,j]-pos[a1,j]
            # a2 is translated by the nearest integer amount of box vectors.
            # If diff is between [-0.5,0.5], nothing happens, [-1.5,-0.5] or [0.5,1.5]
            # will be shifted by one, etc. It should always end up within [-0.5,0.5] of a1  
            pos[a2,j] -= rint(diff)
                


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def make_whole(floating[:,:] pos, int[:,:] bonds, floating[:,:] box):
    """ Makes molecules in sel whole using bond info from bonds.
        Coordinates are first switched to reciprocal space, then passed to
        _iterate_bonds to actually make them whole. The fixed positions are
        switched back to normal space and set as the current positions for sel.
    """
    cdef floating[:,:] invbox
    cdef floating[:,:] unitpos
    # Inverse box
    invbox = np.linalg.inv(box)
    unitpos = np.empty_like(pos)
    # Transfer coordinates to reciprocal space (relative to box vec)
    # Atoms are also put back into box. This is useful to force starter
    # atoms to always be in box 
    cmatmul3D_mod(pos, invbox,unitpos)
    # Fix box in reciprocal space
    iterate_bonds(unitpos, bonds)
    # Bring fixed positions back to real space
    cmatmul3D(unitpos, box, pos)

    return pos


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def wrap_mols(floating[:,:] pos, floating[:] weights, int[:] mols, int nmols, floating[:,:] box):
    """ Puts molecule COMs back in box.
        Parameters:
            pos:     Positions of atoms. float[n,3].
            weights: weights of atoms (masses). float[n].
            mols:    mol indices of each atom. int[n].
            nmols:   number of molecules, np.unique(mols)=np.arange(nmols). int.
            box:     box vectors. float[3,3].
        Returns:
            pos:     positions in pos are moved, and a reference to the original array is returned. float[n,3].
    """
    cdef floating[:,:] invbox
    cdef floating[:,:] molcoms
    cdef floating[:,:] moltrans
    cdef floating[:]   molwsums
    cdef floating[:,:] unitpos  = np.empty_like(pos)
    cdef floating[:]   tmp
    cdef int i,j,k,m,oint
    cdef bint mv

    if(floating is double):
        molcoms  = np.zeros((nmols,3), dtype=np.float64) 
        moltrans = np.zeros((nmols,3), dtype=np.float64) 
        molwsums = np.zeros((nmols),   dtype=np.float64)
        tmp      = np.empty( 3,        dtype=np.float64)
    else: 
        molcoms  = np.zeros((nmols,3), dtype=np.float32) 
        moltrans = np.zeros((nmols,3), dtype=np.float32) 
        molwsums = np.zeros((nmols),   dtype=np.float32)
        tmp      = np.empty( 3,        dtype=np.float32)

    # Inverse box
    invbox = np.linalg.inv(box)

    for i in range(pos.shape[0]):
        for j in range(3):
            molcoms[mols[i],j] += weights[i]*pos[i,j]
        
        molwsums[mols[i]] += weights[i]

    for i in range(nmols):
        for j in range(3):
            molcoms[i,j] /= molwsums[i]

        mv=False
        for m in range(3):
            tmp[m] = 0
            for k in range(3):
                tmp[m] += molcoms[i,k]*invbox[k,m]
            oint = int(tmp[m])
            if(oint!=0 or tmp[m]<0):
                mv=True
                tmp[m] -= oint
                if(tmp[m]<0):
                    tmp[m] += 1.0

        if(mv):
            for m in range(3):
                for k in range(3):
                    moltrans[i,m] += tmp[k]*box[k,m]
                moltrans[i,m] -= molcoms[i,m]

    for i in range(pos.shape[0]):
        for j in range(3):
            pos[i,j] += moltrans[mols[i],j]


    return pos
