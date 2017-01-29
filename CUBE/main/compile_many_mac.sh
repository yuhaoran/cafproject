##!/bin/bash
mkdir -p many
mkdir -p init/many

grep universe parameters.f90

for i in {1..2}; do # compile many universes
  echo compile uninverse$i

  cd init/
  source compile_init_mac.sh
  mv a.out many/ic$i.x

  cd ..
  source compile_caf_mac.sh
  mv cafcube.x many/cube$i.x

  # change file parameter.f90 to next version ('sed' has system dependency)
  sed -i '' 's/universe'"$i"'/universe'"$((i+1))"'/g' parameters.f90
done

# reset parameters.f90 ('sed' has system dependency)
sed -i '' 's/universe'"$((i+1))"'/universe1/g' parameters.f90

grep universe parameters.f90

echo 'done'
