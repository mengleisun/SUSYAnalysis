#!/bin/bash

#For data
#ptbins=( 30 32 34 36 38 40 45 50 55 60 65 70 80 90 120 150 100000)
#ptbins=( 30 32 34 36 38 40 45 50 55 60 65 70 80 90 120 150 180 200 100000)
ptbins=( 200 230 260 100000)
#For MC
#ptbins=( 30 60 100000)

for (( i = 0 ; i < ${#ptbins[@]} ; i++ ))
do
	j=$((i+1))
  root -b -q "FitJetFake.C+(${ptbins[$i]},${ptbins[$j]},1,0,1.4442)"
  #root -b -q "newFit.C+(${ptbins[$i]},${ptbins[$j]})"
done
