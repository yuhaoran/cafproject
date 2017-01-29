module purge
module load gcc/5.2.0 openmpi/gcc/1.8.3 use.experimental caf/gcc/5.2.0-openmpi
module load fftw/3.3.0-gcc-openmpi

module list
rm -f *.o *.out *.mod

FC=caf
XFLAG='-O3 -cpp -march=native -mcmodel=medium'
OFLAG=${XFLAG}' -c'

$FC $OFLAG ../parameters.f90

$FC $OFLAG -Dpenfft_4x penfft_config.f90

$FC $OFLAG voronoi.f90

$FC $XFLAG $FFTFLAG voronoi.o parameters.o penfft_config.o

