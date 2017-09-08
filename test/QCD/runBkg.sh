#!/bin/csh
root -b -q	"analysis_qcdBkg.C+(1, 40,70,0,70,4)"
root -b -q	"analysis_sig.C+(1, 40,70,0,70)"

root -b -q	"analysis_qcdBkg.C+(1, 40,70,70,100,4)"
root -b -q	"analysis_sig.C+(1, 40,70,70,100)"

root -b -q	"analysis_qcdBkg.C+(1, 40,70,100,200,4)"
root -b -q	"analysis_sig.C+(1, 40,70,100,200)"

#root -b -q	"analysis_VGBkg.C+(1, 40,70,0,70)"
#root -b -q	"analysis_VGBkg.C+(1, 40,70,70,100)"
#root -b -q	"analysis_VGBkg.C+(1, 40,70,100,200)"

