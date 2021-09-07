# -*- coding: utf-8 -*-

"""
This module is meant to be modified by user to set the config variables.
"""

# Add any import you need for function_val
import numpy as np


########################### Remote options #####################################
# remote dir corresponding to dir of this notebook
remote_dir="/scratch/project_2004581/functional_sampling_tool/prod"
# remote name, either "host" or "user@host". Must be setup for passwordless connect.
remote_name="mahti"
# Set dirs/files/patterns to exclude from rsync command
rsync_excludes = ["naive", "initial", "templates", "fval.npy", "*.ipynb", "dump"]


########################### sbatch template variables ##########################
email= "santeri.e.paajanen@helsinki.fi"
account="project_2004581"


########################### Sampling variables
# number of simulations per epoch
N = 16


########################### Running simulations
gmx = "gmx"
# input files for grompping
mdp   = "mdrun.mdp"
topol = "topol.top"
ndx   = "index_grompp.ndx"


############################## Function calcs ##################################
# To save memory we only load this selection
# Should be a valid mdtraj selection string
# mdtraj selection string https://www.mdtraj.org/1.9.5/atom_selection.html
select_str = "protein and residue 638 and not (name =~ 'H.*')"

# Minimum and maximum values to sample from between. None to ignore boundary, "start" string to use the
# initial value of the starting structure
minval=0.79
maxval=1.4


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
