module purge
module load intel/16.0.3 intelmpi/5.0.3.048 fftw/3.3.4-intel-impi
module list
#/opt/fftw/3.3.5-gcc-5.4.0-openmpi-2.0.0/lib/
#/opt/fftw/3.3.5-gcc-5.4.0-openmpi-2.0.0/include/

rm -f *.o *.out *.mod *~

FC=ifort
XFLAG='-O3 -xHost -fpp -mcmodel=medium -coarray=shared -qopenmp -DFFTFINE'
OFLAG=${XFLAG}' -c'
FFTFLAG='-I'${SCINET_FFTW_INC}' ''-L'${SCINET_FFTW_LIB}' -lfftw3f -lm -ldl'

echo FFTFLAG:
echo $FFTFLAG

$FC $OFLAG ../main/parameters.f90
$FC $OFLAG ../main/pencil_fft.f90 $FFTFLAG
$FC $OFLAG initial_conditions.f90 $FFTFLAG

$FC $XFLAG parameters.o pencil_fft.o initial_conditions.o $FFTFLAG

export FOR_COARRAY_NUM_IMAGES=1
