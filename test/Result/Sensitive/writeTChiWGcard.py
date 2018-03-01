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

for inputchan in range(1, n_channels + 1):
  file_template = open('./counting_exp_XXX_YYY_' + str(inputchan) + '.txt', 'r')
  lines = [line for line in file_template.readlines()]

  for xbin in range(1, h_SUSYmass.GetXaxis().GetNbins() + 1):
      if(h_SUSYmass.GetBinContent( xbin) <= 0):
          continue
      mass1 = h_SUSYmass.GetXaxis().GetBinCenter(xbin)
  
      file_out = open('/uscms_data/d3/mengleis/work/CMSSW_7_4_7/src/HiggsAnalysis/CombinedLimit/DataCard/Sensitive/cards/counting_tchiwg_' + str(int(mass1)) + '_' + str(inputchan) +  '.txt', 'w')
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
                  or re.search('XSS', l)):
              for k in range(1, n_channels+1):
                  n_nom = h_rates['h_chan' + str(k) + '_rate_nom'].GetBinContent(xbin)
                  e_nom = h_rates['h_chan' + str(k) + '_rate_nom'].GetBinError(xbin)
                  e_jes = h_rates['h_chan' + str(k) + '_syserr_jes'].GetBinContent(xbin)
                  e_jer = h_rates['h_chan' + str(k) + '_syserr_jer'].GetBinContent(xbin)
                  e_esf = h_rates['h_chan' + str(k) + '_syserr_esf'].GetBinContent(xbin)
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
                      l = l.replace('NSC' + str(k) + ' ', str(round(n_nom, 3)))
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
  
          file_out.write(l)
  
      file_out.close()
  
