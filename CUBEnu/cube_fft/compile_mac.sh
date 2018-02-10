rm -f *.o *.out

gfortran -O3 -cpp -fcoarray=single -mcmodel=medium test.f90 -I/usr/local/include/ -L/usr/local/lib/ -lfftw3f -lm -ldl
#gfortran -fopenmp -cpp -fcoarray=single test_parallel.f90 -I/usr/local/include/ -L/usr/local/lib/ -lfftw3f -lm -ldl
