module purge
module load gcc/5.2.0 openmpi/gcc/1.8.3 use.experimental caf/gcc/5.2.0-openmpi
module load fftw/3.3.0-gcc-openmpi

module list

rm -f *.o *.out *.mod *~

FC=caf
XFLAG='-O3 -cpp -march=native -mcmodel=medium'
OFLAG=${XFLAG}' -c'
FFTFLAG='-I'${SCINET_FFTW_INC}' ''-L'${SCINET_FFTW_LIB}' -lfftw3f -lm -ldl'

echo $FFTFLAG

$FC $OFLAG ../parameters.f90

$FC $OFLAG $FFTFLAG -Dpenfft_4x penfft_config.f90

$FC $OFLAG $FFTFLAG penfft_fine.f90

$FC $OFLAG $FFTFLAG powerspectrum.f90

$FC $OFLAG $FFTFLAG displacement.f90

$FC $XFLAG $FFTFLAG displacement.o parameters.o penfft_config.o penfft_fine.o powerspectrum.o

