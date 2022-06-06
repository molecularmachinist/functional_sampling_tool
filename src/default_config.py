# -*- coding: utf-8 -*-

"""
This module is NOT meant to be modified by user.
It provides default values for config, and will be rewritten in updates.
"""


########################### Remote options #####################################
# Set dirs/files/patterns to exclude from rsync command
rsync_excludes = ["config.py", "initial", "templates", "fval_data.npz", ".*.xtc_offsets.npz", "figs", "fst", "src"]



########################### Running simulations
gmx = "gmx"
# input files for grompping
mdp   = "mdrun.mdp"
topol = "topol.top"
ndx   = "index_grompp.ndx"
maxwarn = 0

############################## Function calcs ##################################
# index file for selections
index_file = None

############################## Trajectrory transformations #####################

# whether to make molecules whole from being broken over pbc
unwrap_mols = False
# The selection to make whole
unwrap_sel = "protein"
# A selection of atoms to use as "starters" for unwrapping. Unless mols_in_box=True,
# these atoms are guaranteed to end up inside the box. Makes no difference to outcome
# if mols_in_box=True. None means to use the smallest indexed atom of each mol.
unwrap_starters = None
# Put the molecule COM back into the box (ignored if not unwrapping).
# Only considers atoms in unwrap_sel
mols_in_box = False

# whether to translate clustering coordinates to initial structure
clust_centre = True
# whether to fit clustering coordinates to initial structure, rotationally AND translationally
clust_superpos = False


############################## Advanced options ################################

# List of reps to ignore. Should be a list of tuples (epc,rep),
# e.g. to ignore rep 1 of ep 11 and rep 3 of ep 5: ignore_reps = [(11,1),(5,3)]
ignore_reps = []
# Like above, but just a list of epochs to ignore,
# e.g. to ignore epochs 6 and 7: ignore_epcs = [6,7]
ignore_epcs = []


# whether to add one to maxwarn on the subsequent epochs (due to atom names changing)
maxwarn_add = False

# Histogram building
# At least this much data in each bin of the histogram
data_per_bin = 100
# Maximum amount of bins between boundaries
maxbins = 100
# Minimum amount to choose a frame from. If less frames are in the bin, this many closest frames in value will be used.
# Must be higher or equal to N if allow_choice_duplicates=True
minchoice = data_per_bin
# Allow choosing the same frame more than once in the same epoch.
# In some edge cases of multiple choices within a region of low sampling,
# the minchoice values might overlap and still produce duplicates.
allow_choice_duplicates=False

# Cluster histogram
clust_data_per_bin=1000
clust_maxbins=10

# Number of epoch before clustering
epochs_pre_clust=3
#Max number of clusters per bin
maxclust=15
# choose at most this fraction of choices from clustering
clust_choice_frac=0.5
# Tolerace for number of clusters in clustering
clust_tol=0.1

# These options allow running multiple runs with separate configs suing the same data
# Name of archivefile
npz_file_name = "fval_data.npz"
# number of first rep
first_rep_num = 1
