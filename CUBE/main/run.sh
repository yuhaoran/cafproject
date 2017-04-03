source ../utilities/module_load_intel.sh 
module list

cd ../utilities/
make clean
make
cafrun -np 1 ./ic.x

cd ../main/
make clean
make
#cafrun -np 1 ./cafcube.x
