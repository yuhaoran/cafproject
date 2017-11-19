##!/bin/bash

# Note that this is in Fortran convension
n1=5  # first universe to compile
n2=8  # last universe to compile

# make directories for executables
mkdir -p many/

# check if opath=='../output/universe1/'
echo check if opath=='../output/universe1/'
grep universe ../main/parameters.f90

# change opath name
sed -i 's/universe1/universe'"$n1"'/g' ../main/parameters.f90

for ((i=$n1; i<=$n2; i++))
do
  echo start to compile uninverse$i

  cd ../utilities/
  source module_load_intel.sh
  make clean
  make
  mv ic.x ../batch/many/ic_universe$i.x
  #mv dsp.x ../batch/many/dsp_universe$i.x
  #mv convert.x ../batch/many/convert_universe$i.x
  mv cicpower.x ../batch/many/cicpower_universe$i.x  

  cd ../main/
  #make clean
  #make
  #mv cafcube.x ../batch/many/cube_universe$i.x

  # change file parameter.f90 to next version ('sed' has system dependency)
  sed -i 's/universe'"$i"'/universe'"$((i+1))"'/g' ../main/parameters.f90
  cd ../batch/
done

# reset parameters.f90 ('sed' has system dependency)
sed -i 's/universe'"$((n2+1))"'/universe1/g' ../main/parameters.f90

# check again if opath=='../output/universe1/'
echo check again if opath=='../output/universe1/'
grep universe ../main/parameters.f90

echo Done
