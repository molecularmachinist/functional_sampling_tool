#!/bin/bash
## Any line starting with "##" will not be copied.
##
## From the below, change AT LEAST the parts with left pointing arrows (and remove the arrows)
## This should be fine for a run on Mahti. Change other variables as needed.
## Writing {epoch_number} anywhere will be substituted with the epoch number.
#SBATCH --nodes=20
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=36:00:00
#SBATCH --partition=medium
#SBATCH --account=project_<1234567>                  <----------------
#SBATCH -o output.txt
#SBATCH -J fst_epoch{epoch_number}
#SBATCH --mail-type=END
#SBATCH --mail-user=<erkki.esimerkki@domain.com>     <----------------



# Load the necessary modules
## Currently these work to load the 2020.5 version of gromacs on Mahti
module load gcc/9.4.0
module load openmpi/4.1.2
module load gromacs/2020.5


# OMP environment variables
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores

# This runs all epochs using the multidir option of gromacs.
## If you have enough resources per single job to do so, this should be the best option.
## Otherwise use e.g. array jobs
srun gmx_mpi mdrun  -deffnm mdrun -multidir rep* -maxh 36

