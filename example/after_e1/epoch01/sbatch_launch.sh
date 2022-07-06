#!/bin/bash
#SBATCH --nodes=20
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=36:00:00
#SBATCH --partition=medium
#SBATCH --account=project_12345
#SBATCH -o output.txt
#SBATCH -J ampa_epoch1
#SBATCH --mail-type=END
#SBATCH --mail-user=erkki.esimerkki@helsinki.fi

set -e

for ((i=1; i<17; i++)); do
  cd rep$(printf "%02d" $i)
  gmx mdrun  -deffnm mdrun
  cd ..
done
