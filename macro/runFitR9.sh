#!/bin/bash
#set ptbins =(30 32 34 36 38 40 44 48 52 56 60 66 72 78 86 90 100 120 150 10000)
ptbins=(10 20 30 40 50 200)
EtaBins=(0 0.8 1.444) 

low=60
high=120
for (( i = 0 ; i <=5 ; i++ ))
do
	j=$((i+1))
  root -b -q "FitR9.C+(0, ${ptbins[$i]},${ptbins[$j]}, -2.5,   -2.0 )" 
  root -b -q "FitR9.C+(0, ${ptbins[$i]},${ptbins[$j]}, -2.0,   -1.566 )" 
  root -b -q "FitR9.C+(0, ${ptbins[$i]},${ptbins[$j]}, -1.566, -1.444 )" 
  root -b -q "FitR9.C+(0, ${ptbins[$i]},${ptbins[$j]}, -1.444, -0.8 )" 
  root -b -q "FitR9.C+(0, ${ptbins[$i]},${ptbins[$j]}, -0.8,   0)" 
  root -b -q "FitR9.C+(0, ${ptbins[$i]},${ptbins[$j]}, 0,  0.8)" 
  root -b -q "FitR9.C+(0, ${ptbins[$i]},${ptbins[$j]}, 0.8, 1.444)" 
  root -b -q "FitR9.C+(0, ${ptbins[$i]},${ptbins[$j]}, 1.444, 1.566)"
	root -b -q "FitR9.C+(0, ${ptbins[$i]},${ptbins[$j]}, 1.566, 2.0)"
	root -b -q "FitR9.C+(0, ${ptbins[$i]},${ptbins[$j]}, 2.0,   2.5)" 
done

for (( m = 0 ; m <=6 ; m++ ))
do
	n=$((m+1))
  root -b -q "FitR9.C+(1, ${ptbins[$m]},${ptbins[$n]}, -2.5,   -2.0 )" 
  root -b -q "FitR9.C+(1, ${ptbins[$m]},${ptbins[$n]}, -2.0,   -1.566 )" 
  root -b -q "FitR9.C+(1, ${ptbins[$m]},${ptbins[$n]}, -1.566, -1.444 )" 
  root -b -q "FitR9.C+(1, ${ptbins[$m]},${ptbins[$n]}, -1.444, -0.8 )" 
  root -b -q "FitR9.C+(1, ${ptbins[$m]},${ptbins[$n]}, -0.8,   0)" 
  root -b -q "FitR9.C+(1, ${ptbins[$m]},${ptbins[$n]}, 0,  0.8)" 
  root -b -q "FitR9.C+(1, ${ptbins[$m]},${ptbins[$n]}, 0.8, 1.444)" 
  root -b -q "FitR9.C+(1, ${ptbins[$m]},${ptbins[$n]}, 1.444, 1.566)"
	root -b -q "FitR9.C+(1, ${ptbins[$m]},${ptbins[$n]}, 1.566, 2.0)"
	root -b -q "FitR9.C+(1, ${ptbins[$m]},${ptbins[$n]}, 2.0,   2.5)" 
done
