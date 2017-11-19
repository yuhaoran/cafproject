##!/bin/bash

n1=1  # first universe to submit
n2=8  # last universe to submit

#for i in {1..2}; do # submit many universes

sed -i 's/universe1/universe'"$n1"'/g' submit_one_universe.qsub

for ((i=$n1; i<=$n2; i++))
do
  echo submit uninverse$i
  grep universe submit_one_universe.qsub
  qsub submit_one_universe.qsub

  sed -i 's/universe'"$i"'/universe'"$((i+1))"'/g' submit_one_universe.qsub

done

# reset submit_one_universe ('sed' has system dependency)

sed -i 's/universe'"$((n2+1))"'/universe1/g' submit_one_universe.qsub

echo 'done'
