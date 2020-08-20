#!/bin/bash
#SBATCH -p short
#SBATCH -t 10
#SBATCH -n 4

echo "efit.sbatch script --begin" >> efit.log

##cd $SLURM_SUBMIT_DIR

module list >> efit.log 2>&1

srun --mpi=pmi2 -n 4 efit_exe 65 >> efit.log 2>&1
