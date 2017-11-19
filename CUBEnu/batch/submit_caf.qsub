#!/bin/bash
# MOAB/Torque submission script for SciNet GPC
#
#PBS -l nodes=8:ppn=8,walltime=1:00:00
#PBS -N CUBE

# DIRECTORY TO RUN - $PBS_O_WORKDIR is directory job was submitted from
cd $PBS_O_WORKDIR

# EXECUTION COMMAND; -np = nodes*ppn
export FOR_COARRAY_NUM_IMAGES=8

cd ../utilities/
export I_MPI_PROCESS_MANAGER=mpd
source module_load_intel.sh
module load caf/intel/any

#mpirun ../batch/many/ic_universe1.x > ../batch/many/log_ic_universe1

cd ../main/
cafrun -np 8 -N 1 ../batch/many/cube_universe1.x > ../batch/many/log_cube_universe1

#cd ../utilities/
#../batch/many/dsp_universe1.x > ../batch/many/log_dsp_universe1
#../batch/many/convert_universe1.x > ../batch/many/log_convert_universe1

cd ../batch/
