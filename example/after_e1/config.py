# -*- coding: utf-8 -*-

"""
This module is meant to be modified by user to set the config variables.
"""

# Add any import you need for function_val
import numpy as np


########################### Remote options #####################################
# remote dir corresponding to dir of this notebook
# Just make a local subfolder
import os
remote_dir=os.getcwd()+"/remote_dir"
# remote name, either "host" or "user@host". Must be setup for passwordless connect.
remote_name="localhost"
# Set dirs/files/patterns to exclude from rsync command
rsync_excludes = ["config.py", "initial", "templates", "fval_data.npz", "dump", "figs", "fst", "src", "remote_dir"]


########################### sbatch template variables ##########################
email= "erkki.esimerkki@helsinki.fi"
account="project_12345"


########################### Sampling variables
# number of simulations per epoch, only affects new epochs
N = 16


########################### Running simulations
gmx = "gmx"
# input files for grompping
mdp   = "mdrun.mdp"
topol = "topol.top"
ndx   = "index_grompp.ndx"
maxwarn=1

############################## Function calcs ##################################
# To make next run faster we save this selection to disk
# The coordinates of the selction are used in function_val
# Should be a valid mdtraj selection string OR a group in index_file
# mdtraj selection string https://www.mdtraj.org/1.9.5/atom_selection.html
select_str = "protein and resid 638 and not ( type H )"
# Same as above, but selection for clustering
select_str_clust = "protein and name CA"
# index file for selections
index_file = None

# Minimum and maximum values to sample from between. None to ignore boundary, "start" string to use the
# initial value of the starting structure
minval=7.9
maxval=14


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
mols_in_box = False




############################## Advanced options ################################
maxwarn_add=True


# Histogram building
# At least this much data in each bin of the histogram
data_per_bin = 100
# Maximum amount of bins between boundaries
maxbins = 100

# Cluster histogram
clust_data_per_bin=1000
clust_maxbins=10

#Max number of clusters per bin
maxclust=15
# choose at most this fraction of choices from clustering
clust_choice_frac=0.5
# Tolerace for number of clusters in clustering
clust_tol=0.1
# Number of epoch before clustering
epochs_pre_clust=3
