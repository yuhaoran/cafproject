#!/bin/bash
for i in {1..5}; do
sed -i '16cprint*,\"Test'''$i'''\"' test.f90 # replace the 16th line of test.f90
gfortran test.f90 -o test
./test
done
sed -i '16cprint*,\"Test1\"' test.f90 # resetting the 16th line of test.f90
