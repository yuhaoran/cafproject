source ../utilities/module_load_mac.sh

#rm -rf ../output/universe1/*

cd ../utilities/
make clean
make ic.x
make cicpower.x
./ic.x

# Main N-body code
cd ../main/
make clean
make
./cafcube.x

cd ../utilities/
./cicpower.x


#cd ../utilities/
#./qspace.x
#./corr.x
#./potential.x # get -grad(phi_q)
#./Jmatrix.x
#./detJ.x # get inv(J), and apply -inv(J)grad(phi_q)

cd ../main/
