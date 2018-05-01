#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0:05:00
#SBATCH --job-name test_ic
#SBATCH --output=output_ic.txt

cd $SLURM_SUBMIT_DIR
source module_load_niagara_intel.sh
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
./ic.x
