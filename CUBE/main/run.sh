source ../utilities/module_load_intel.sh 
module list

cd ../utilities/
make clean
make
cafrun -np 1 -N 1 ./ic.x

cd ../main/
make clean
make
cafrun -np 1 -N 1 ./cafcube.x

cd ../utilities/
cafrun -np 1 -N 1 ./dsp.x

cd ../main/
