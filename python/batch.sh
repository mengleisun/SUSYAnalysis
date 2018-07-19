
#!/usr/bin/env python

import os,sys

scriptFile = open('run_analysis_signal.sh', 'w')
scriptFile.write('#!/bin/bash\n')
scriptFile.write('cp -r /afs/cern.ch/work/m/msun/Analysis/SUSYAnalysis/include/ .\n')
scriptFile.write('cp -r /afs/cern.ch/work/m/msun/Analysis/SUSYAnalysis/test/ .\n')
scriptFile.write('source /afs/cern.ch/sw/lcg/contrib/gcc/4.9/x86_64-slc6-gcc49-opt/setup.sh\n') 
scriptFile.write('source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.32/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh\n')
scriptFile.write('cd test\n')
scriptFile.write('root -b -q "analysis_trigger.C+(\\"job_DoubleEG_Run2015C_PR_25ns.root\\")"\n')
scriptFile.write('cp *.root /afs/cern.ch/work/m/msun/Analysis/SUSYAnalysis/test/\n')
