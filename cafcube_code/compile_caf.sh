module purge
module load gcc/5.4.0
module load openmpi/2.0.0-gcc-5.4.0
module load fftw/3.3.5-gcc-5.4.0-openmpi-2.0.0

make clean
make
