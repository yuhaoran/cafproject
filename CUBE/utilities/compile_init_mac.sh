rm *.o *.out

FC=mpif90
XFLAG='-O3 -cpp -fcoarray=single -mcmodel=medium'
OFLAG=${XFLAG}' -c'
FFTFLAG='-lfftw3f -lm -ldl'

$FC $OFLAG ../main/parameters.f90

$FC $OFLAG $FFTFLAG -Dpenfft_4x penfft_config.f90

$FC $OFLAG $FFTFLAG penfft_fine.f90

$FC $OFLAG $FFTFLAG cross_power.f90

$FC $OFLAG $FFTFLAG initial_conditions.f90

$FC $XFLAG $FFTFLAG initial_conditions.o parameters.o penfft_config.o cross_power.o penfft_fine.o
