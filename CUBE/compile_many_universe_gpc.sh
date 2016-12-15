##!/bin/bash
mkdir -p many
mkdir -p init/many

grep universe parameters.f90

for i in {1..2}; do # compile many universes
  echo compile uninverse$i

  cd init/
  module purge
  module load gcc/5.2.0 openmpi/gcc/1.8.3 use.experimental caf/gcc/5.2.0-openmpi
  module load fftw/3.3.0-gcc-openmpi
  source compile_init_scinet_gcc5.sh
  mv a.out many/ic_universe$i.x
  source compile_disp_scinet_gcc5.sh
  mv a.out many/dsp_universe$i.x

  cd ..
  module purge
  module load intel/16.0.3 intelmpi/5.0.3.048 fftw/3.3.4-intel-impi
  source compile_caf_scinet_intel.sh
  mv cafcube.x many/cube_universe$i.x

  # change file parameter.f90 to next version ('sed' has system dependency)
  sed -i 's/universe'"$i"'/universe'"$((i+1))"'/g' parameters.f90
done

# reset parameters.f90 ('sed' has system dependency)
sed -i 's/universe'"$((i+1))"'/universe1/g' parameters.f90

grep universe parameters.f90

echo 'done'
