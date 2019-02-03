#module purge
#module load intel/16.0.3 intelmpi/5.0.3.048 fftw/3.3.4-intel-impi
#module list

rm -f *.o *.out

gfortran -O3 -cpp -fopenmp -mcmodel=medium test.f90 -lfftw3f -I/usr/local/include/ -L/usr/local/lib/ -lm -ldl -o fft1.x
gfortran -O3 -cpp -fopenmp -mcmodel=medium test_parallel.f90 -lfftw3_omp -lfftw3f_omp  -I/usr/local/include/ -L/usr/local/lib/ -lm -ldl -o fft32.x
export OMP_NUM_THREADS=32
