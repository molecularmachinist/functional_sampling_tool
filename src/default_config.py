# -*- coding: utf-8 -*-

"""
This module is NOT meant to be modified by user.
It provides default values for config, and will be rewritten in updates.
"""


########################### Remote options #####################################
# Set dirs/files/patterns to exclude from rsync command
rsync_excludes = ["fval_data.npz",
                  ".mdrun.xtc_offsets.lock", ".mdrun.xtc_offsets.npz"]


# Running simulations
gmx = "gmx"
# input files for grompping
mdp = "mdrun.mdp"
topol = "topol.top"
ndx = "index_grompp.ndx"
sbatch = "sbatch_launch.sh"
# file for -r option for grompping. False to not use the option, "initial" to use
# the inital structure, "start" to use the starting structrue of the repetition and any
# other string to give the file manually (relative to the repetition folder)
restraint_file = False
maxwarn = 0

############################## Function calcs ##################################
# Initial structure file(s)
initial_struct = "start.gro"
# index file for selections
index_file = None

############################## Trajectrory transformations #####################

# whether to centre by a single atom before making anything whole
precenter = False
# The atom to use for precentering. By default use the one closest to box centre.
precenter_atom = None
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

# the random number seed. None means a totally random run.
rng_seed = None

# List of reps to ignore. Should be a list of tuples (epc,rep),
# e.g. to ignore rep 1 of ep 11 and rep 3 of ep 5: ignore_reps = [(11,1),(5,3)]
ignore_reps = []
# Like above, but just a list of epochs to ignore,
# e.g. to ignore epochs 6 and 7: ignore_epcs = [6,7]
ignore_epcs = []

# Options to handle using only some of the data
# The function values will be recalculated for every epoch and repetition if these options are changed
# Amount of data (frames) to ignore from start of each repetition
ignore_from_start = 0
# The stride for only using every n'th frame
stride = 1

# Histogram building
# At least this much data in each bin of the histogram
data_per_bin = 100
# Maximum amount of bins between boundaries
maxbins = 100
# The window to use for histogram smoothing
smooth_window = 10
# Minimum amount to choose a frame from. If less frames are in the bin, this many closest frames in value will be used.
# Must be higher or equal to N if allow_choice_duplicates=True
minchoice = data_per_bin
# maximum height to choose from, as fraction from minimum found height to maximum
# (maxh-minh)*choice_crit+minh will be the highest allowed histogram bin.
choice_crit = 0.5
# Allow choosing the same frame more than once in the same epoch.
# In some edge cases of multiple choices within a region of low sampling,
# the minchoice values might overlap and still produce duplicates.
allow_choice_duplicates = False
# A dictionary of keyword arguments for scipy.signal.find_peaks https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html
peak_options = {"width": 5, "distance": 10}
# Make plot histogram as cumulative, or each epoch separately
cumulative_histogram = True
# Maximum number of (latest) epochs to plot in histogram, for cumulative only changes
# how many are added to legend.
histogram_max_epochs = 15


# Cluster histogram
clust_data_per_bin = 1000
clust_maxbins = 10

# Number of epoch before clustering
epochs_pre_clust = 3
# Max number of clusters per bin
maxclust = 15
# choose at most this fraction of choices from clustering
clust_choice_frac = 0.5
# Tolerace for number of clusters in clustering
clust_tol = 0.1

# These options allow running multiple runs with separate configs suing the same data
# Name of archivefile
npz_file_name = "fval_data.npz"
# number of first rep
first_rep_num = 1
# figure output dir
fig_output_dir = "figs"
