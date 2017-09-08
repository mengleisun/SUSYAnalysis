#!/bin/bash
#set ptbins =(30 35 40 45 50 55 60 65 70 75 80 90 100 120 150 10000)
ptbins=(30 35 40 45 50 55 60 65 70 75 80 90 100 120 150 180 10000)
#EtaBins=(0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0 1.05 1.1 1.15 1.2 1.25 1.3 1.35 1.4 1.45 1.5)
EtaBins=(1.55 1.60 1.65 1.70 1.75 1.80 1.85 1.90 1.95 2.00 2.05 2.10 2.15 2.20 2.25 2.30 2.35 2.40 2.45 2.5)
VertexBins=(1 4 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 10000)

#for (( i = 0 ; i < ${#ptbins[@]} ; i++ ))
#do
#	j=$((i+1))
#  root -b -q "FitKer.C+(0, 0, ${ptbins[$i]},${ptbins[$j]}, 40,140)" 
#done 
#
#for (( i = 0 ; i < ${#ptbins[@]} ; i++ ))
#do
#	j=$((i+1))
#  root -b -q "FitKer.C+(0, 1, ${ptbins[$i]},${ptbins[$j]}, 40,140)" 
#done 

for (( i = 0 ; i < ${#EtaBins[@]} ; i++ ))
do
	j=$((i+1))
  root -b -q "FitKer.C+(1, 0, ${EtaBins[$i]},${EtaBins[$j]}, 40,140)" 
done 

for (( i = 0 ; i < ${#EtaBins[@]} ; i++ ))
do
	j=$((i+1))
  root -b -q "FitKer.C+(1, 1, ${EtaBins[$i]},${EtaBins[$j]}, 40,140)" 
done 

#for (( i = 0 ; i < ${#VertexBins[@]} ; i++ ))
#do
#	j=$((i+1))
#  root -b -q "FitKer.C+(2, 0, ${VertexBins[$i]},${VertexBins[$j]}, 40,140)" 
#done 
#
#for (( i = 0 ; i < ${#VertexBins[@]} ; i++ ))
#do
#	j=$((i+1))
#  root -b -q "FitKer.C+(2, 1, ${VertexBins[$i]},${VertexBins[$j]}, 40,140)" 
#done 
