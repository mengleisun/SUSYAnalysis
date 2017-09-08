#!/bin/csh
#set ptbins =(30 32 34 36 38 40 44 48 52 56 60 66 72 78 86 90 100 120 150 10000)
set ptbins =(35 40 45 50 55 60 65 70 75 80 90 100 120 150 10000)
set EtaBins =(0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0 1.05 1.1 1.15 1.2 1.25 1.3 1.35 1.4 1.45 1.5)
set VertexBins =(1 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 48 10000)

set low = 60
set high = 120
set i = 1
set j = 2
while ( $i <= 14 )
  root -b -q "FitKer.C+(0, 0, ${ptbins[$i]},${ptbins[$j]}, $low, $high)" 
  @ i+= 1
	@ j+= 1
end

set m = 1
set n = 2
while ( $m <= 14 )
  root -b -q "FitKer.C+(0, 1, ${ptbins[$m]},${ptbins[$n]}, $low, $high)" 
  @ m+= 1
	@ n+= 1
end

set a = 1
set b = 2
while ( $a <= 29 )
  root -b -q "FitKer.C+(1, 0, ${EtaBins[$a]},${EtaBins[$b]}, $low, $high)" 
  @ a+= 1
	@ b+= 1
end

set c = 1
set d = 2
while ( $c <= 29 )
  root -b -q "FitKer.C+(1, 1, ${EtaBins[$c]},${EtaBins[$d]}, $low, $high)" 
  @ c+= 1
	@ d+= 1
end

set e = 1
set f = 2
while ( $e <= 21 )  #vtx
  root -b -q "FitKer.C+(2, 0, ${VertexBins[$e]},${VertexBins[$f]},$low, $high)" 
  @ e+= 1
	@ f+= 1
end

set g = 1
set h = 2
while ( $g <= 21 )  #vtx
  root -b -q "FitKer.C+(2, 1, ${VertexBins[$g]},${VertexBins[$h]},$low, $high)" 
  @ g+= 1
	@ h+= 1
end
