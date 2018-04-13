#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --time=0:05:00
#SBATCH --job-name test_main
#SBATCH --output=output_main.txt

cd $SLURM_SUBMIT_DIR
source ../utilities/module_load_niagara_intel.sh
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#srun -N 1  ./main.x
mpirun -bind-to node -np 1 ./main.x 
