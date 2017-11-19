source ../utilities/module_load_intel.sh 
module list

cd ../utilities/
make clean
make
cafrun -np 8 ./ic.x

cd ../main/
make clean
make
cafrun -np 8 ./cafcube.x

cd ../utilities/
cafrun -np 8 ./cicpower.x

cd ../main/
