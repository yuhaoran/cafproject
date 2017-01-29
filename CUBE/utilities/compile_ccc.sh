rm -f *.o *.out *.mod
FC=mpif90
XFLAG='-O3 -cpp -fcoarray=single -mcmodel=medium'
OFLAG=${XFLAG}' -c'
FFTFLAG='-I -L -lfftw3f -lm -ldl'

$FC $OFLAG ../parameters.f90
$FC $OFLAG -Dpenfft_4x penfft_config.f90
$FC $OFLAG $FFTFLAG penfft_fine.f90
$FC $OFLAG $FFTFLAG powerspectrum.f90

$FC $OFLAG $FFTFLAG ccc.f90
$FC $XFLAG $FFTFLAG ccc.o parameters.o penfft_config.o penfft_fine.o powerspectrum.o
