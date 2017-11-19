$FC $OFLAG ../main/parameters.f90 -DFFTFINE
$FC $OFLAG ../main/pencil_fft.f90 $FFTFLAG
$FC $OFLAG powerspectrum.f90 $FFTFLAG
$FC $OFLAG displacement.f90 $FFTFLAG

$FC $XFLAG parameters.o pencil_fft.o powerspectrum.o displacement.o $FFTFLAG

export FOR_COARRAY_NUM_IMAGES=1
