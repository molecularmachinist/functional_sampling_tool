# -*- coding: utf-8 -*-

"""
This module is meant to be modified by user to set the config variables.
"""

# Add any import you need for function_val
import numpy as np


########################### Remote options #####################################
# remote dir corresponding to dir of this project
# For the example case we just make a local subfolder
# DO NOT use in production
import os
remote_dir=os.getcwd()+"/../remote_dir"
# remote name, either "<host>" or "<user>@<host>".
remote_name="localhost"
# Set dirs/files/patterns to exclude from rsync command (uncomment to override defaults)
#rsync_excludes = ["config.py", "initial", "templates", "fval_data.npz", "dump", "figs", "fst", "src", "remote_dir"]


########################### Sampling variables
# number of simulations per epoch, only affects new epochs
N = 16

########################### Running simulations
maxwarn=1

############################## Function calcs ##################################
# To make next run faster we save this selection to disk
# The coordinates of the selction are used in function_val
# Should be a valid MDAnalysis selection string OR a group in index_file
# For MDAnalysis selection strings, see
# https://docs.mdanalysis.org/stable/documentation_pages/selections.html
select_str = "protein and resid 1 14 and name CA"
# Same as above, but selection for clustering
select_str_clust = "protein and name CA"
# index file for selections
index_file = None

# Minimum and maximum values to sample from between. None to ignore boundary, "start" string to use the
# initial value of the starting structure
minval=0
maxval="start"


def function_val(positions):
    """
      Write here your analysis function. Positions will be
        numpy array of shape (n,m,3) for n frames of m atoms.
        Note that m is the number of atoms in the selection,
        not the whole trajectory.
      The function should return a numpy array of shape (n).
    """
    # In our case the only two atoms chosen, so m=2
    p1 = positions[:,0,:]
    p2 = positions[:,1,:]
    # Calculate distance between atoms
    dist = np.linalg.norm(p1-p2, axis=-1)
    # p1 and p2 were shape (n,3), and the bove is calculated over the last (length 3) axis.
    # Therefore dist is now shape (n)
    return dist


############################## Unwrapping options ##############################

# whether to make molecules whole from being broken over pbc
unwrap_mols = True
# The selection to make whole
unwrap_sel = "protein"
# A selection of atoms to use as "starters" for unwrapping. Unless mols_in_box=True,
# these atoms are guaranteed to end up inside the box. Makes no difference to outcome
# if mols_in_box=True. None means to use the smallest indexed atom of each mol.
unwrap_starters = None
# Put the molecule COM back into the box (ignored if not unwrapping).
# Only considers atoms in unwrap_sel
mols_in_box = True


############################## Coordinate fitting ##############################

# whether to fit clustering coordinates to initial structure, rotationally AND translationally
clust_superpos = True



############################## Advanced options ################################

# Histogram building
# At least this much data in each bin of the histogram
data_per_bin = 50
