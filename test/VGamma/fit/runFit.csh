#!/bin/bash
for i in {0..100}
do
  root -b -q "Fitfractionmg.C+($i,40,70,0,1000,4)"
done
for j in {0..100}
do
  root -b -q "Fitfractioneg.C+($j,40,70,0,1000,4)"
done
#
#set i = 0
#while ( $i < 1000 )
#  root -b -q "Fitfraction.C+($i,50,70,4)"
#  @ i+= 1
#end
#
#set j = 0
#while ( $j < 1000 )
#  root -b -q "Fitfraction.C+($j,70,100,4)" 
#  @ j+= 1
#end
#
#set n = 0
#while ( $n < 100 )
#  root -b -q "Fitfraction.C+($n,100,20000,4)" 
#  @ n+= 1
#end

#set i = 0
#root -b -q "Fitfractionmg.C+(0,40,70,0,1000,4)"
#root -b -q "FitMass.C+(0,40,70,0,1000,4)"
#root -b -q "Fitfractionmg.C+(0,40,70,0,1000,1000)"
#root -b -q "FitMass.C+(0,40,70,0,1000,1000)"
#
root -b -q "Fitfractioneg.C+(0,40,70,0,30,4)"
root -b -q "Fitfractioneg.C+(0,40,70,30,50,4)"
root -b -q "Fitfractioneg.C+(0,40,70,50,70,4)"
root -b -q "Fitfractioneg.C+(0,40,70,70,200,4)"
