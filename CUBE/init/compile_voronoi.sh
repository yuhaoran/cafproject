rm *.o *.out *.mod
mpif90 -O3 -cpp -fcoarray=single   -c ../parameters.f90
mpif90 -O3 -cpp -fcoarray=single   -Dpenfft_4x -c penfft_config.f90
mpif90 -O3 -cpp -fcoarray=single  -c voronoi.f90 

mpif90 -O3 -cpp -fcoarray=single voronoi.o parameters.o penfft_config.o
