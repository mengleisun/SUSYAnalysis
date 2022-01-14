#!/bin/bash
#To study the possible p T dependence of the Vγ scale factors, we perform the template fit in four lepton pT bins:
#Fitfractioneg.C+($i,METlow,METhigh,ptLow,ptHigh,isolation cut)
#for i in {0..1000}
#do
#  root -b -q "Fitfractioneg.C+($i,40,70,0,-1,4)"
#done
for i in {0..1000}
do
  root -b -q "Fitfractioneg.C+($i,40,70,0,1000,4)"
done

#for i in {0..1000}
#do
#  root -b -q "Fitfractioneg.C+($i,40,70,0,50,4)"
#done
#for i in {0..1000}
#do
#  root -b -q "Fitfractioneg.C+($i,40,70,50,70,4)"
#done
#for i in {0..1000}
#do
#  root -b -q "Fitfractioneg.C+($i,40,70,70,100,4)"
#done
#for i in {0..1000}
#do
#  root -b -q "Fitfractioneg.C+($i,40,70,100,1000,4)"
#done
#
#for i in {0..1000}
#do
#  root -b -q "Fitfractionmg.C+($i,40,70,0,50,4)"
#done
#for i in {0..1000}
#do
#  root -b -q "Fitfractionmg.C+($i,40,70,50,70,4)"
#done
#for i in {0..1000}
#do
#  root -b -q "Fitfractionmg.C+($i,40,70,70,100,4)"
#done
#for i in {0..1000}
#do
#  root -b -q "Fitfractionmg.C+($i,40,70,100,1000,4)"
#done
#
