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


########################### sbatch template variables ##########################
email= "erkki.esimerkki@helsinki.fi"
account="project_12345"


########################### Sampling variables
# number of simulations per epoch, only affects new epochs
N = 16


########################### Running simulations
#maxwarn=1

############################## Function calcs ##################################
# To make next run faster we save this selection to disk
# The coordinates of the selction are used in function_val
# Should be a valid MDAnalysis selection string OR a group in index_file
# For MDAnalysis selection strings, see
# https://docs.mdanalysis.org/stable/documentation_pages/selections.html
select_str = "protein and resid 638 and not ( type H )"
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
#maxwarn_add=True
