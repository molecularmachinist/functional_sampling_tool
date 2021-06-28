#!/bin/bash
#SBATCH --nodes=1
#SBATCH --hint=multithread
#SBATCH --ntasks-per-node=256
#SBATCH --cpus-per-task=1
#SBATCH --time=36:00:00
#SBATCH --partition=medium
#SBATCH --account=project_2004581
#SBATCH -o rep%02a/output.txt
#SBATCH -J ampa_equil
#SBATCH --mail-type=END
#SBATCH --mail-user=santeri.e.paajanen@helsinki.fi
#SBATCH --array=1-16


# stop on any error (easier debug)
set -e

hostname

module load gcc/9.3.0
module load openmpi/4.0.3
module load gromacs/2020.5




export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores

cd rep$(printf "%02d" ${SLURM_ARRAY_TASK_ID})
pwd

# Generated by CHARMM-GUI (http://www.charmm-gui.org) v3.5, modified by Santeri Paajanen



init=../solv/system_solv_ions
rest_prefix=../solv/system_solv_ions
mini_prefix=step6.0_minimization
equi_prefix=step6.%d_equilibration
prod_prefix=step7_production
prod_step=step7
index=../solv/index_grompp
topol=../solv/topol

mdrun_cmd="srun gmx_mpi mdrun"
grompp_cmd="gmx_mpi grompp"

# Minimization
# In the case that there is a problem during minimization using a single precision of GROMACS, please try to use 
# a double precision of GROMACS only for the minimization step.
${grompp_cmd} -f ../charmm_model/${mini_prefix}.mdp -o ${mini_prefix}.tpr -c ${init}.gro -r ${rest_prefix}.gro -p ${topol}.top -n ${index}.ndx
${mdrun_cmd}  -deffnm ${mini_prefix}

# Equilibration
cntmax=6

for ((cnt=1;cnt<=cntmax;cnt++)); do
    pcnt=$(( cnt - 1 ))
    istep=`printf ${equi_prefix} ${cnt}`
    pstep=`printf ${equi_prefix} ${pcnt}`
    if [[ ${cnt} == 1 ]]; then pstep=${mini_prefix}; fi

    ${grompp_cmd} -f ../charmm_model/${istep}.mdp -o ${istep}.tpr -c ${pstep}.gro -r ${rest_prefix}.gro -p ${topol}.top -n ${index}.ndx
    ${mdrun_cmd} -deffnm ${istep}
done
