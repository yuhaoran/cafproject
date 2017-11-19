$FC $OFLAG ../main/parameters.f90 -DFFTFINE
$FC $OFLAG ../main/pencil_fft.f90 $FFTFLAG
$FC $OFLAG initial_conditions.f90 $FFTFLAG

$FC $XFLAG parameters.o pencil_fft.o initial_conditions.o $FFTFLAG

export FOR_COARRAY_NUM_IMAGES=1
