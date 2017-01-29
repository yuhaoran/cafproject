rm *.o *.out
mpif90 -O3 -cpp -fcoarray=single   -c ../main/parameters.f90 -o parameters.o

mpif90 -O3 -cpp -fcoarray=single   -Dpenfft_4x -c penfft_config.f90

mpif90 -O3 -cpp -fcoarray=single   -c penfft_fine.f90 -o penfft_fine.o -lfftw3f -I/usr/local/include/ -L/usr/local/lib/ -lm -ldl

mpif90 -O3 -cpp -fcoarray=single  -c cross_power.f90 -lfftw3f -I/usr/local/include/ -L/usr/local/lib/ -lm -ldl

mpif90 -O3 -cpp -fcoarray=single  -c displacement.f90 -o displacement.o  -lfftw3f -I/usr/local/include/ -L/usr/local/lib/ -lm -ldl

mpif90 -O3 -cpp -fcoarray=single displacement.o parameters.o penfft_config.o penfft_fine.o cross_power.o -o a.out  -lfftw3f -I/usr/local/include/ -L/usr/local/lib/ -lm -ldl
