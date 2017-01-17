rm -f *.o *.out

mpif90 -m64 -O3 -cpp -fcoarray=single -mcmodel=medium test.f90 -lfftw3f -lm -ldl
