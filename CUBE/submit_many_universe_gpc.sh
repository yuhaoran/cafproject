##!/bin/bash

for i in {1..2}; do # submit many universes
  echo submit uninverse$i
  grep universe submit_one_universe.qsub
  qsub submit_one_universe.qsub

  sed -i 's/universe'"$i"'/universe'"$((i+1))"'/g' submit_one_universe.qsub

done

# reset submit_one_universe ('sed' has system dependency)

sed -i 's/universe'"$((i+1))"'/universe1/g' submit_one_universe.qsub

echo 'done'
