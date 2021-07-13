#!/bin/bash
#SBATCH --nodes=20
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=36:00:00
#SBATCH --partition=medium
#SBATCH --account={account}
#SBATCH -o output.txt
#SBATCH -J ampa_epoch{i}
#SBATCH --mail-type=END
#SBATCH --mail-user={email}


hostname

module load gcc/9.3.0
module load openmpi/4.0.3
module load gromacs/2020.5




export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores


srun gmx_mpi mdrun  -deffnm mdrun -multidir rep* -maxh 36

