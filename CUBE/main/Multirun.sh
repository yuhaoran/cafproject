##!/bin/bash
for i in {1..2}; do # compile 2 universe
sed -i '6ccharacter(*),parameter :: opath=\".\/output\/universe'''$i'''\/\"' parameters.f90 # change the output universe id in file parameters.90
sed -i '19s/cafcube\.x/cafcube'''$i'''\.x/g' Makefile # change the 'cafcube.x' as 'cafcuben.x', n=1,2,3... in Makefile line 19
sed -i '21s/cafcube\.x/cafcube'''$i'''\.x/g' Makefile #change the 'cafcube.x' as 'cafcuben.x', n=1,2,3... in Makefile line 21
echo $i
# compile start
cd init
source compile_init.sh
./a.out
cd ..
source compile_caf.sh
# compile end 
cp "cafcube$i.x" "multirun/cafcube$i.x" # copy 'cafcuben.x' to path /multirun 
#./cafcube.x
sed -i '19s/cafcube'''$i'''\.x/cafcube\.x/g' Makefile # reset Makefile
sed -i '21s/cafcube'''$i'''\.x/cafcube\.x/g' Makefile # reset Makefile
done
sed -i '6ccharacter(*),parameter :: opath=\".\/output\/universe1\/\"' parameters.f90 # reset parameters.f90