import re
import os
import ROOT
import argparse

n_channels = 36
syst_names = ['jes','jer','esf','scale','eleshape','jetshape','qcdshape','xs','lumi','isr']

susy_in = ROOT.TFile('../signalTree_TChiWG.root','read')
h_SUSYmass = susy_in.Get('h_chan1_rate_nom')
h_rates = {}
for ich in range(1,n_channels+1):
    h_rates['h_chan' + str(ich) + '_rate_nom'] = susy_in.Get('h_chan' + str(ich) + '_rate_nom')
    for sys in syst_names:
        h_rates['h_chan' + str(ich) + '_syserr_' + sys] = susy_in.Get('h_chan' + str(ich) + '_syserr_' + sys)
h_rates['tchiwg_h_syserr_PU'] = susy_in.Get('tchiwg_h_syserr_PU')

file_template = open('./counting_exp_XXX_YYY.txt', 'r')
lines = [line for line in file_template.readlines()]

low_p = 2 
high_p = 0

for xbin in range(1, h_SUSYmass.GetXaxis().GetNbins() + 1):
    if(h_SUSYmass.GetBinContent( xbin) <= 0):
        continue
    mass1 = h_SUSYmass.GetXaxis().GetBinCenter(xbin)

    for br in range(0,51):
     
        bf = br*0.02 
        brweight = 2*bf*(1-bf)
        #file_out = open('/uscms_data/d3/mengleis/work/CMSSW_7_4_7/src/HiggsAnalysis/CombinedLimit/DataCard/TChiWG_Scan/card/counting_tchiwg_' + str(int(mass1))  + '_BF_' + str(int(br*2)) +  '.txt', 'w')
        file_out = open('/uscmst1b_scratch/lpc1/3DayLifetime/mengleis/TChiWG/counting_tchiwg_' + str(int(mass1))  + '_BF_' + str(int(br*2)) +  '.txt', 'w')
        avg_jes = 0
        avg_jer = 0
        avg_esf = 0
        n_nonzero = 0
        for k in range(1, n_channels+1):
            n_nom = h_rates['h_chan' + str(k) + '_rate_nom'].GetBinContent(xbin)
            e_jes = h_rates['h_chan' + str(k) + '_syserr_jes'].GetBinContent(xbin)
            e_jer = h_rates['h_chan' + str(k) + '_syserr_jer'].GetBinContent(xbin)
            e_esf = h_rates['h_chan' + str(k) + '_syserr_esf'].GetBinContent(xbin)
            if n_nom != 0.0:
                avg_jes +=  e_jes / n_nom
                avg_jer +=  e_jer / n_nom 
                avg_esf +=  e_esf / n_nom
                n_nonzero += 1
        avg_jes = avg_jes/n_nonzero
        avg_jer = avg_jer/n_nonzero
        avg_esf = avg_esf/n_nonzero
    
        for l in lines:
            if (re.search('NSC', l)
                    or re.search('NEVT', l)
                    or re.search('STSC', l)
                    or re.search('STAT', l)
                    or re.search('JESS', l)
                    or re.search('JERS', l)
                    or re.search('ESF', l)
                    or re.search('XSS', l)
                    or re.search('PUS', l)):
                for k in range(1, n_channels+1):
                    n_nom = h_rates['h_chan' + str(k) + '_rate_nom'].GetBinContent(xbin)
                    e_nom = h_rates['h_chan' + str(k) + '_rate_nom'].GetBinError(xbin)
                    e_jes = h_rates['h_chan' + str(k) + '_syserr_jes'].GetBinContent(xbin)
                    e_jer = h_rates['h_chan' + str(k) + '_syserr_jer'].GetBinContent(xbin)
                    e_esf = h_rates['h_chan' + str(k) + '_syserr_esf'].GetBinContent(xbin)
                    unc_PU = h_rates['tchiwg_h_syserr_PU'].GetBinContent(xbin) + 1.0
                    unc_stat = 1.0
                    unc_jes = 1.0 
                    unc_jer = 1.0
                    unc_esf = 1.0
                    unc_xs  = 1.0
                    if n_nom != 0.0:
                        unc_stat = 1.0 +  e_nom / n_nom 
                        unc_jes = 1.0 +  e_jes / n_nom
                        unc_jer = 1.0 +  e_jer / n_nom 
                        unc_esf = 1.0 +  e_esf / n_nom
                    else:
                        unc_stat = 2.0
                        unc_jes = 1.0 + avg_jes 
                        unc_jer = 1.0 + avg_jer 
                        unc_esf = 1.0 + avg_esf 
    
                    if re.search('NSC', l):
                        l = l.replace('NSC' + str(k) + ' ', str(round(n_nom*brweight, 3)))
                    elif re.search('STSC', l):
                        l = l.replace('STSC' + str(k)+ ' ', str(round(unc_stat, 3)))
                    elif re.search('JESS', l):
                        l = l.replace('JESS' + str(k)+ ' ', str(round(unc_jes, 3)))
                    elif re.search('JERS', l):
                        l = l.replace('JERS' + str(k)+ ' ', str(round(unc_jer, 3)))
                    elif re.search('ESF', l):
                        l = l.replace('ESF' + str(k)+ ' ', str(round(unc_esf, 3)))
                    elif re.search('XSS', l):
                        l = l.replace('XSS' + str(k)+ ' ', str(round(unc_xs , 3)))
                    elif re.search('PUS', l):
                        l = l.replace('PUS' + ' ', str(round(unc_PU , 3)))
                        if(n_nom > 0):
                          if(unc_PU < low_p):
                              low_p = unc_PU
                          if(unc_PU > high_p):
                              high_p = unc_PU
    
            file_out.write(l)
    
        file_out.close()

print low_p, high_p
