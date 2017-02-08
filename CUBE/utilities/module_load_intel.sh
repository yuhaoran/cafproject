module purge
module load intel/16.0.3 intelmpi/5.0.3.048 fftw/3.3.4-intel-impi
#module list

export FC='ifort'
export XFLAG='-O3 -fpp -mcmodel=medium -coarray=shared'
# GPC nodes do not support intel AVX instructions "-xHost"
export OFLAG=${XFLAG}' -c'
export FFTFLAG='-I'${SCINET_FFTW_INC}' ''-L'${SCINET_FFTW_LIB}' -lfftw3f -lm -ldl'

export FOR_COARRAY_NUM_IMAGES=1
