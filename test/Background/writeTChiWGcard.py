import re
import os
import ROOT
import argparse

limdir = '/uscms_data/d3/mengleis/work/CMSSW_7_4_7/src/HiggsAnalysis/CombinedLimit/tchiwg/'

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--MET1', dest='METbin1', type=int)
parser.add_argument('--MET2', dest='METbin2', type=int)
args = parser.parse_args()


met1 = args.METbin1
met2 = args.METbin2


n_channels = 18
n_processes= 6
pro_names = ['TChiWG','elefakepho', 'jetfakepho', 'qcdfakelep','VGamma','rare']
syst_names = ['jes','jer','esf','scale','eleshape','jetshape','xs','lumi']


file_in = ROOT.TFile('./SignalSystematic_' + str(met1) + '_' + str(met2) + '.root', 'read')
susy_in = ROOT.TFile('./signalTree_TChiWG_' + str(met1) + '_' + str(met2) + '.root','read')
h_SUSYmass = susy_in.Get('h_chan1_rate_nom')
h_rates = {}
for hname in pro_names:
    if(hname == 'TChiWG'):
        for ich in range(1,n_channels+1):
            h_rates['h_chan' + str(ich) + '_rate_nom'] = susy_in.Get('h_chan' + str(ich) + '_rate_nom')
            for sys in syst_names:
                h_rates['h_chan' + str(ich) + '_syserr_' + sys] = susy_in.Get('h_chan' + str(ich) + '_syserr_' + sys)
    else:
        h_rates['h_' + hname + '_norm'] = file_in.Get('h_' + hname + '_norm')
        for sys in syst_names:
            h_rates['h_' + hname + '_syserr_' + sys] = file_in.Get('h_' + hname + '_syserr_' + sys)


for xbin in range(1, h_SUSYmass.GetXaxis().GetNbins() + 1):
    if(h_SUSYmass.GetBinContent( xbin) <= 0):
        continue
    mass1 = h_SUSYmass.GetXaxis().GetBinCenter(xbin)

    file_out = open(limdir+'counting_tchiwg_' + str(int(mass1)) + '_' + str(met1) + '_' + str(met2) +  '.txt', 'w')


    file_out.write("imax 18  number of channels\n")
    file_out.write("jmax 5  number of backgrounds\n")
    file_out.write("kmax *  number of nuisance parameters\n")
    file_out.write("------------\n")
    file_out.write("bin             1       2       3			4			5			6			7			8			9       10      11      12	  13		14	  15	  16		17		18\n")
    file_out.write("observation     1       1       1     1     1     1     1     1     1       1       1       1     1     1     1     1     1     1\n")
    file_out.write("------------\n")
		
    file_out.write('bin   ')
    for i in range(1,n_channels+1):
        for j in range(0,n_processes):
            file_out.write('%d      ' %i)
    file_out.write('\n')

    file_out.write('process   ')
    for i in range(0,n_channels):
        for p in pro_names: 
            file_out.write(p+" ")
    file_out.write('\n')

    file_out.write('process   ')
    for i in range(0,n_channels):
        for j in range(0,n_processes):
            file_out.write('%d      ' %j)
    file_out.write('\n')

    file_out.write('rate   ')
    for ich in range(1,n_channels+1):
        for p in pro_names:
            evtrate = 0
            if(p == 'TChiWG'):
                evtrate = h_rates['h_chan' + str(ich) + '_rate_nom'].GetBinContent(xbin)
            else:
                evtrate = h_rates['h_' + p + '_norm'].GetBinContent(ich)
            file_out.write(str(round(evtrate, 3)) + '   ')
    file_out.write('\n')
	  
    file_out.write('------------\n')

    for syst in syst_names:
        file_out.write(syst + '   lnN ')
        for ich in range(1,n_channels+1):
            for p in pro_names:
                syserr = 1;
                if(p == 'TChiWG'):
                    if(h_rates['h_chan' + str(ich) + '_rate_nom'].GetBinContent(xbin)==0):
                        file_out.write(str(round(1, 3)) + '  ')
                    else:
                        syserr = h_rates['h_chan' + str(ich) + '_syserr_' + sys].GetBinContent(xbin)/h_rates['h_chan' + str(ich) + '_rate_nom'].GetBinContent(xbin)
                        if(syserr < 0):
                             file_out.write('-  ')
                        else:
                             syserr = syserr + 1
                             file_out.write(str(round(syserr, 3)) + '  ')
                else:
                    if(h_rates['h_' + p + '_norm'].GetBinContent(ich)==0):
                        file_out.write(str(round(1, 3)) + '  ')
                    else:
                        syserr = h_rates['h_' + p + '_syserr_' + syst].GetBinContent(ich)/h_rates['h_' + p + '_norm'].GetBinContent(ich)
                        if(syserr < 0):
                             file_out.write('-  ')
                        else:
                             syserr = syserr + 1
                             file_out.write(str(round(syserr, 3)) + '  ')
        file_out.write('\n')
 

    file_out.close()


