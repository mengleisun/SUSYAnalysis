#!/bin/csh
#set ptbins =(30 32 34 36 38 40 44 48 52 56 60 66 72 78 86 90 100 120 150 10000)
set ptbins =(10 20 35 50 90 150 500 10000)
set EtaBins =(0 0.8 1.444) 

set low = 60
set high = 120
set i = 1
set j = 2
while ( $i <= 6 )
  root -b -q "FitKer.C+(0, ${ptbins[$i]},${ptbins[$j]}, -2.5,   -2.0 )" 
  root -b -q "FitKer.C+(0, ${ptbins[$i]},${ptbins[$j]}, -2.0,   -1.566 )" 
  root -b -q "FitKer.C+(0, ${ptbins[$i]},${ptbins[$j]}, -1.566, -1.444 )" 
  root -b -q "FitKer.C+(0, ${ptbins[$i]},${ptbins[$j]}, -1.444, -0.8 )" 
  root -b -q "FitKer.C+(0, ${ptbins[$i]},${ptbins[$j]}, -0.8,   0)" 
  root -b -q "FitKer.C+(0, ${ptbins[$i]},${ptbins[$j]}, 0,  0.8)" 
  root -b -q "FitKer.C+(0, ${ptbins[$i]},${ptbins[$j]}, 0.8, 1.444)" 
  root -b -q "FitKer.C+(0, ${ptbins[$i]},${ptbins[$j]}, 1.444, 1.566)"
	root -b -q "FitKer.C+(0, ${ptbins[$i]},${ptbins[$j]}, 1.566, 2.0)"
	root -b -q "FitKer.C+(0, ${ptbins[$i]},${ptbins[$j]}, 2.0,   2.5)" 
  @ i+= 1
	@ j+= 1
end

set m = 1
set n = 2
while ( $m <= 6 )
  root -b -q "FitKer.C+(1, ${ptbins[$m]},${ptbins[$n]}, -2.5,   -2.0 )" 
  root -b -q "FitKer.C+(1, ${ptbins[$m]},${ptbins[$n]}, -2.0,   -1.566 )" 
  root -b -q "FitKer.C+(1, ${ptbins[$m]},${ptbins[$n]}, -1.566, -1.444 )" 
  root -b -q "FitKer.C+(1, ${ptbins[$m]},${ptbins[$n]}, -1.444, -0.8 )" 
  root -b -q "FitKer.C+(1, ${ptbins[$m]},${ptbins[$n]}, -0.8,   0)" 
  root -b -q "FitKer.C+(1, ${ptbins[$m]},${ptbins[$n]}, 0,  0.8)" 
  root -b -q "FitKer.C+(1, ${ptbins[$m]},${ptbins[$n]}, 0.8, 1.444)" 
  root -b -q "FitKer.C+(1, ${ptbins[$m]},${ptbins[$n]}, 1.444, 1.566)"
	root -b -q "FitKer.C+(1, ${ptbins[$m]},${ptbins[$n]}, 1.566, 2.0)"
	root -b -q "FitKer.C+(1, ${ptbins[$m]},${ptbins[$n]}, 2.0,   2.5)" 
  @ m+= 1
	@ n+= 1
end

