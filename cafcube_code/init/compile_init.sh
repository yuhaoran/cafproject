#/opt/fftw/3.3.5-gcc-5.4.0-openmpi-2.0.0/lib/
#/opt/fftw/3.3.5-gcc-5.4.0-openmpi-2.0.0/include/
module purge
module load gcc/5.4.0
module load openmpi/2.0.0-gcc-5.4.0
module load fftw/3.3.5-gcc-5.4.0-openmpi-2.0.0

module list

rm -f *.o *.out *.mod *~

mpif90 -O3 -cpp -fcoarray=single -mcmodel=medium -c ../parameters.f90 -o parameters.o

mpif90 -O3 -cpp -fcoarray=single -mcmodel=medium -Dpenfft_4x -c penfft_config.f90

mpif90 -O3 -cpp -fcoarray=single -mcmodel=medium -c penfft_fine.f90 -o penfft_fine.o -lfftw3f -I/opt/fftw/3.3.5-gcc-5.4.0-openmpi-2.0.0/include/ -L/opt/fftw/3.3.5-gcc-5.4.0-openmpi-2.0.0/lib/ -lm -ldl

mpif90 -O3 -cpp -fcoarray=single -mcmodel=medium -c powerspectrum.f90 -lfftw3f -I/opt/fftw/3.3.5-gcc-5.4.0-openmpi-2.0.0/include/ -L/opt/fftw/3.3.5-gcc-5.4.0-openmpi-2.0.0/lib/ -lm -ldl

mpif90 -O3 -cpp -fcoarray=single -c initial_conditions.f90 -o initial_conditions.o  -lfftw3f -I/opt/fftw/3.3.5-gcc-5.4.0-openmpi-2.0.0/include/ -L/opt/fftw/3.3.5-gcc-5.4.0-openmpi-2.0.0/lib/ -lm -ldl

mpif90 -O3 -cpp -fcoarray=single initial_conditions.o parameters.o penfft_config.o powerspectrum.o penfft_fine.o -o a.out  -lfftw3f -I/opt/fftw/3.3.5-gcc-5.4.0-openmpi-2.0.0/include/ -L/opt/fftw/3.3.5-gcc-5.4.0-openmpi-2.0.0/lib/ -lm -ldl
