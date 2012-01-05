#!/bin/bash

echo R0,i
for ((i=50;$i<151;i=$i+1))
do
  scaling=$(echo "scale=2;$i/100" | bc)
  r0=$(bin/solve_ode.x -o params/ode.prm -m params/model.prm --ic-file init/full.init --scaling $scaling  --tmax 1000)
  eq=$(cat try1.final | sed 's/\t/\n/g'  | awk '{sum += $1}END{print sum}')
  echo $r0,$eq
done
