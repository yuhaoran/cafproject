#/opt/fftw/3.3.5-gcc-5.4.0-openmpi-2.0.0/lib/
#/opt/fftw/3.3.5-gcc-5.4.0-openmpi-2.0.0/include/
module purge
module load intel/intel-16 intelmpi/5.1.3

rm *.o *.out
ifort -O3 -fpp -coarray=shared -mcmodel=medium -c ../parameters.f90 -o parameters.o

ifort -O3 -fpp -coarray=shared -mcmodel=medium -Dpenfft_4x -c penfft_config.f90

ifort -O3 -fpp -coarray=shared -mcmodel=medium -c penfft_fine.f90 -o penfft_fine.o -lfftw3f -I/cita/h/home-1/haoran/fftw3_lib/include/ -L/cita/h/home-1/haoran/fftw3_lib/lib/ -lm -ldl

ifort -O3 -fpp -coarray=shared -mcmodel=medium -c powerspectrum.f90 -lfftw3f -I/cita/h/home-1/haoran/fftw3_lib/include/ -L/cita/h/home-1/haoran/fftw3_lib/lib/ -lm -ldl

ifort -O3 -fpp -coarray=shared -mcmodel=medium -c initial_conditions.f90 -o initial_conditions.o  -lfftw3f -I/cita/h/home-1/haoran/fftw3_lib/include/ -L/cita/h/home-1/haoran/fftw3_lib/lib/ -lm -ldl

ifort -O3 -fpp -coarray=shared -mcmodel=medium initial_conditions.o parameters.o penfft_config.o powerspectrum.o penfft_fine.o -o a.out  -lfftw3f -I/cita/h/home-1/haoran/fftw3_lib/include/ -L/cita/h/home-1/haoran/fftw3_lib/lib/ -lm -ldl

export FOR_COARRAY_NUM_IMAGES=1
