module purge
module load gcc/5.2.0 openmpi/gcc/1.8.3 use.experimental caf/gcc/5.2.0-openmpi
module load fftw/3.3.0-gcc-openmpi

rm -f *.o *.out

caf -O3 -cpp -march=native -mcmodel=medium -c ../parameters.f90 -o parameters.o

caf -O3 -cpp -march=native -mcmodel=medium -c rhof.f90 -o rhof.o

caf -O3 -cpp -march=native rhof.o parameters.o -mcmodel=medium -o a.out
