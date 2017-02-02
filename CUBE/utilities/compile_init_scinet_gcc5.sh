module purge
module load gcc/5.2.0 openmpi/gcc/1.8.3 use.experimental caf/gcc/5.2.0-openmpi
module load fftw/3.3.0-gcc-openmpi
module list

rm -f *.o *.out *.mod *~

FC=caf
XFLAG='-O3 -cpp -march=native -mcmodel=medium -DFFTFINE'
OFLAG=${XFLAG}' -c'
FFTFLAG='-I'${SCINET_FFTW_INC}' ''-L'${SCINET_FFTW_LIB}' -lfftw3f -lm -ldl'

echo FFTFLAG:
echo $FFTFLAG

$FC $OFLAG ../main/parameters.f90
$FC $OFLAG ../main/pencil_fft.f90 $FFTFLAG
$FC $OFLAG initial_conditions.f90 $FFTFLAG

$FC $XFLAG parameters.o pencil_fft.o initial_conditions.o $FFTFLAG
