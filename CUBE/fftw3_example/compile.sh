module purge
module load intel/16.0.3 intelmpi/5.0.3.048 fftw/3.3.4-intel-impi
module list

rm -f *.o *.out

ifort -O3 -xHost -fpp -mcmodel=medium test.f90 -lfftw3f -I -L -lm -ldl
