export FC='mpif90'
export XFLAG='-O3 -cpp -fcoarray=single -fopenmp'
export OFLAG=${XFLAG}' -c'
export FFTFLAG='-lfftw3f -lm -ldl'
