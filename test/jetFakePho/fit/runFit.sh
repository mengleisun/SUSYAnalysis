#!/bin/bash
ptbins=( 30 32 34 36 38 40 45 50 55 60 65 70 80 90 120 150 180 220 10000)
for (( i = 0 ; i < ${#ptbins[@]} ; i++ ))
do
	j=$((i+1))
  #root -b -q "FitJetFake.C+(${ptbins[$i]},${ptbins[$j]})"
  root -b -q "FitJetFake.C+(${ptbins[$i]},${ptbins[$j]},1)"
done

#for (( i = 0 ; i < ${#ptbins[@]} ; i++ ))
#do
#	j=$((i+1))
#  #root -b -q "FitJetFake.C+(${ptbins[$i]},${ptbins[$j]})"
#  root -b -q "FitJetFake.C+(${ptbins[$i]},${ptbins[$j]},2)"
#done
