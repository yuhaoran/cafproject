source ../utilities/module_load_intel.sh

rm -f *.o *.out

$FC $XFLAG test.f90 $FFTFLAG
export FOR_COARRAY_NUM_IMAGES=8
