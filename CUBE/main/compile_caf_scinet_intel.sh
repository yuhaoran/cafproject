#scinet intel16
module purge
module load intel/16.0.3 intelmpi/5.0.3.048 fftw/3.3.4-intel-impi 
module list

make clean
make -f Makefile_intel

export FOR_COARRAY_NUM_IMAGES=1
