#!/bin/csh
set ptbins =(30 32 34 36 38 40 45 50 55 60 65 70 80 90 120 150 180 10000)
set i = 1
set j = 2
while ( $i <= 17 )
  #root -b -q "FitJetFake.C+(${ptbins[$i]},${ptbins[$j]})"
  root -b -q "FitJetFake.C+(${ptbins[$i]},${ptbins[$j]})"
  @ i+= 1
  @ j+= 1
end

