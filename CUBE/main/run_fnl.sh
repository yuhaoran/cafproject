source ../utilities/module_load_mac.sh

cp parameters.f90_fnl1 parameters.f90
source run_mac.sh > log1

cp parameters.f90_fnl2 parameters.f90
source run_mac.sh > log2

cp parameters.f90_fnl3 parameters.f90
source run_mac.sh > log3

cp parameters.f90_fnl4 parameters.f90
source run_mac.sh > log4

cp parameters.f90_fnl5 parameters.f90
source run_mac.sh > log5
