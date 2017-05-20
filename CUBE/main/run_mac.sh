source ../utilities/module_load_mac.sh

#rm -rf ../output/universe1/*

cd ../utilities/
make clean
make
cd ../main/
make clean
make

cd ../utilities/
./ic.x
cd ../main/
./cafcube.x
cd ../utilities/
./dsp.x

cd ../main/
