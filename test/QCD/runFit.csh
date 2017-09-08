#!/bin/csh
set i = 0
while ( $i < 1 )
  root -b -q "Fitfraction.C+($i,0,70,4)" 
  @ i+= 1
end

set j = 0
while ( $j < 1 )
  root -b -q "Fitfraction.C+($j,70,100,4)" 
  @ j+= 1
end

set n = 0
while ( $n < 1 )
  root -b -q "Fitfraction.C+($n,100,200,4)" 
  @ n+= 1
end

#set i = 0
#root -b -q "Fitfraction.C+($i,0,70,4)" 
#root -b -q "Fitfraction.C+($i,70,100,4)" 
#root -b -q "Fitfraction.C+($i,100,200,4)" 
#root -b -q "Fitfraction.C+($i,0,10000,30)" 
