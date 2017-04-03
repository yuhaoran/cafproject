source ../utilities/module_load_mac.sh

rm -rf ../output/universe1/*

cd ../utilities/
make
./ic.x

cd ../main/
make
./cafcube.x
