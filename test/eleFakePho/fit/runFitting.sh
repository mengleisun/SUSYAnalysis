#!/bin/bash
rm RunnerScript.txt
ptbins=(30 35 40 45 50 55 60 65 70 75 80 90 100 120 150 180 10000)
Etabins=(0.0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0 1.05 1.1 1.15 1.2 1.25 1.3 1.35 1.4 1.45 1.5)
VertexBins=(1 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 10000)

for (( i = 0 ; i < ${#ptbins[@]}-1 ; i++ ))
do
	j=$((i+1))
  echo  1 0 0 ${ptbins[$i]} ${ptbins[$j]} 60 120 >> RunnerScript.txt
  echo  2 0 0 ${ptbins[$i]} ${ptbins[$j]} 60 120 >> RunnerScript.txt
  echo  3 0 0 ${ptbins[$i]} ${ptbins[$j]} 60 120 >> RunnerScript.txt
  echo  1 0 1 ${ptbins[$i]} ${ptbins[$j]} 60 120 >> RunnerScript.txt
  echo  2 0 1 ${ptbins[$i]} ${ptbins[$j]} 60 120 >> RunnerScript.txt
  echo  3 0 1 ${ptbins[$i]} ${ptbins[$j]} 60 120 >> RunnerScript.txt
done 


for (( i = 0 ; i < ${#Etabins[@]}-1 ; i++ ))
do
	j=$((i+1))
  echo  1 1 0 ${Etabins[$i]} ${Etabins[$j]} 60 120 >> RunnerScript.txt
  echo  2 1 0 ${Etabins[$i]} ${Etabins[$j]} 60 120 >> RunnerScript.txt
  echo  3 1 0 ${Etabins[$i]} ${Etabins[$j]} 60 120 >> RunnerScript.txt
  echo  1 1 1 ${Etabins[$i]} ${Etabins[$j]} 60 120 >> RunnerScript.txt
  echo  2 1 1 ${Etabins[$i]} ${Etabins[$j]} 60 120 >> RunnerScript.txt
  echo  3 1 1 ${Etabins[$i]} ${Etabins[$j]} 60 120 >> RunnerScript.txt
done 

for (( i = 0 ; i < ${#VertexBins[@]}-1 ; i++ ))
do
	j=$((i+1))
  echo  1 2 0 ${VertexBins[$i]} ${VertexBins[$j]} 60 120 >> RunnerScript.txt
  echo  2 2 0 ${VertexBins[$i]} ${VertexBins[$j]} 60 120 >> RunnerScript.txt
  echo  3 2 0 ${VertexBins[$i]} ${VertexBins[$j]} 60 120 >> RunnerScript.txt
  echo  1 2 1 ${VertexBins[$i]} ${VertexBins[$j]} 60 120 >> RunnerScript.txt
  echo  2 2 1 ${VertexBins[$i]} ${VertexBins[$j]} 60 120 >> RunnerScript.txt
  echo  3 2 1 ${VertexBins[$i]} ${VertexBins[$j]} 60 120 >> RunnerScript.txt
done 

g++ -fPIC -g -Wall -Wno-maybe-uninitialized ../../../lib/RooDCBShape.o ../../../lib/RooUserPoly.o ../../../lib/libAnaClasses.so ../../../lib/libRooFitClasses.so -L$ROOTSYS -lRooFit -lHtml -lMinuit -lRooFitCore -lRooStats -lHistFactory `root-config --cflags` FitKer.C -o FitKer.exe `root-config --libs` -L$ROOTSYS/lib -lGui -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic
rm configFit_Bw_ker.txt
echo 'useCMSShape 0' >> configFit_Bw_ker.txt 
echo 'useExpo 0'     >> configFit_Bw_ker.txt
echo 'useKer  1'     >> configFit_Bw_ker.txt
echo 'useDY  0'      >> configFit_Bw_ker.txt
echo 'SaveOutput 1'  >> configFit_Bw_ker.txt

rm configFit_DY_ker.txt
echo 'useCMSShape 0' >> configFit_DY_ker.txt 
echo 'useExpo 0'     >> configFit_DY_ker.txt
echo 'useKer  1'     >> configFit_DY_ker.txt
echo 'useDY  1'      >> configFit_DY_ker.txt
echo 'SaveOutput 1'  >> configFit_DY_ker.txt

rm configFit_Bw_expo.txt
echo 'useCMSShape 0' >> configFit_Bw_expo.txt 
echo 'useExpo 1'     >> configFit_Bw_expo.txt
echo 'useKer  0'     >> configFit_Bw_expo.txt
echo 'useDY  0'      >> configFit_Bw_expo.txt
echo 'SaveOutput 1'  >> configFit_Bw_expo.txt

#./FitKer.exe 1 0 0 40 45 60 120


cd /uscms/homes/t/tmishra/work/CMSSW_10_2_22/src

tar -czvf files.tar.gz SUSYAnalysis/test/eleFakePho/fit/FitKer.C SUSYAnalysis/include/RooCMSShape.h SUSYAnalysis/include/RooDCBShape.h SUSYAnalysis/include/analysis_rawData.h SUSYAnalysis/include/analysis_photon.h SUSYAnalysis/include/analysis_muon.h SUSYAnalysis/include/analysis_ele.h SUSYAnalysis/include/analysis_tools.h SUSYAnalysis/include/analysis_mcData.h SUSYAnalysis/include/tdrstyle.C SUSYAnalysis/lib/libAnaClasses.so SUSYAnalysis/lib/libRooFitClasses.so SUSYAnalysis/test/eleFakePho/fit/make.sh SUSYAnalysis/lib/roofit_Dict_rdict.pcm SUSYAnalysis/lib/analysis_photonDict_rdict.pcm SUSYAnalysis/test/eleFakePho/fit/RunnerScript.txt SUSYAnalysis/test/eleFakePho/fit/command_HTcondor.sh SUSYAnalysis/test/eleFakePho/fit/command_HTcondor.py SUSYAnalysis/test/eleFakePho/fit/Run_Analyser_short.py SUSYAnalysis/include/RooUserPoly.h SUSYAnalysis/test/eleFakePho/fit/FitKer.C SUSYAnalysis/lib/RooDCBShape.o SUSYAnalysis/lib/RooUserPoly.o SUSYAnalysis/test/eleFakePho/fit/make.sh SUSYAnalysis/test/eleFakePho/fit/RooDCBShape.h SUSYAnalysis/test/eleFakePho/fit/RooUserPoly.h SUSYAnalysis/lib/* SUSYAnalysis/include/* SUSYAnalysis/src/*

cp files.tar.gz SUSYAnalysis/test/eleFakePho/fit/
