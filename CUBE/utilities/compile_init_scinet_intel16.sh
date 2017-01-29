module purge
module load intel/16.0.3 intelmpi/5.0.3.048 fftw/3.3.4-intel-impi
module list
#/opt/fftw/3.3.5-gcc-5.4.0-openmpi-2.0.0/lib/
#/opt/fftw/3.3.5-gcc-5.4.0-openmpi-2.0.0/include/

rm -f *.o *.out
ifort -O3 -xHost -fpp -coarray=shared -mcmodel=medium -c ../parameters.f90 -o parameters.o

ifort -O3 -xHost -fpp -coarray=shared -mcmodel=medium -Dpenfft_4x -c penfft_config.f90

ifort -O3 -xHost -fpp -coarray=shared -mcmodel=medium -c penfft_fine.f90 -o penfft_fine.o -lfftw3f -I -L -lm -ldl

ifort -O3 -xHost -fpp -coarray=shared -mcmodel=medium -c powerspectrum.f90 -lfftw3f -I -L -lm -ldl

ifort -O3 -xHost -fpp -coarray=shared -mcmodel=medium -c initial_conditions.f90 -o initial_conditions.o  -lfftw3f -I -L -lm -ldl

ifort -O3 -xHost -fpp -coarray=shared -mcmodel=medium initial_conditions.o parameters.o penfft_config.o powerspectrum.o penfft_fine.o -o a.out  -lfftw3f -I -L -lm -ldl

export FOR_COARRAY_NUM_IMAGES=1
