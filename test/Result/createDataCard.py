import re
import os
import ROOT
import argparse

limdir = './'

n_channels = 36
n_processes= 6
pro_names = ['SUSY','elefakepho', 'jetfakepho', 'qcdfakelep','VGamma','rare']
syst_names = ['jes','jer','esf','scale','e_to_pho_syst','j_to_pho_syst','fakelep_shape','xs','lumi','isr']

file_in = ROOT.TFile('./SignalSystematic.root', 'read')
h_rates = {}
for hname in pro_names:
    h_rates['h_' + hname + '_norm'] = file_in.Get('h_' + hname + '_norm')
    for sys in syst_names:
        h_rates['h_' + hname + '_syserr_' + sys] = file_in.Get('h_' + hname + '_syserr_' + sys)
h_rates['h_elefakepho_controlsample'] = file_in.Get('h_elefakepho_controlsample')
h_rates['h_elefakepho_transferfactor'] = file_in.Get('h_elefakepho_transferfactor')
h_rates['h_jetfakepho_controlsample'] = file_in.Get('h_jetfakepho_controlsample')
h_rates['h_jetfakepho_transferfactor'] = file_in.Get('h_jetfakepho_transferfactor')
h_rates['h_qcdfakelep_controlsample'] = file_in.Get('h_qcdfakelep_controlsample')
h_rates['h_qcdfakelep_transferfactor'] = file_in.Get('h_qcdfakelep_transferfactor')

logfile_out= open( 't5wg.log','w')
file_out = open( limdir+'counting_exp_XXX_YYY.txt', 'w')

file_out.write("imax 36 number of channels\n")
file_out.write("jmax 5  number of backgrounds\n")
file_out.write("kmax *  number of nuisance parameters\n")
file_out.write("------------\n")
file_out.write("bin             1  2  3  4  5  6  7  8  9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36   \n")
file_out.write("observation    309  494  85  32  64  45  1  1  5  12  23  20  4  12  7  1  1  0  153  276  67   32   46   32   1  1  4  10  21  14  6  9  4  0  1  3  \n")

file_out.write("------------\n")

file_out.write('{:26s}'.format('bin'))
for i in range(1,n_channels+1):
    for j in range(0,n_processes):
        file_out.write('{:12d} '.format(i))
file_out.write('\n')

file_out.write('{:26s}'.format('process'))
for i in range(0,n_channels):
    for p in pro_names: 
        file_out.write('{:>12s} '.format(p))
file_out.write('\n')

file_out.write('{:26s}'.format('process'))
for i in range(0,n_channels):
    for j in range(0,n_processes):
        file_out.write('{:12d} '.format(j))
file_out.write('\n')

file_out.write('{:26s}'.format('rate'))
for ich in range(1,n_channels+1):
    signal_rate = 0
    bkg_rate = 0
    ratio = 0
    for p in pro_names:
        evtrate = 0
        if(p == 'SUSY'):
            file_out.write('{:>12s} '.format('NSC'+str(ich)))
        else:
            evtrate = h_rates['h_' + p + '_norm'].GetBinContent(ich)
            bkg_rate += evtrate
            file_out.write('{:12.3f} '.format(evtrate))
            logfile_out.write('{:15} {:2d} {:10.5f}\n'.format(p, ich, evtrate))
    ratio = signal_rate/bkg_rate
file_out.write('\n')

file_out.write('------------\n')

low_p = {}
high_p = {}
for p in pro_names:
  low_p[p] = 2
  high_p[p] = 0

for syst in syst_names:
    file_out.write('{:15s} {:3s} {:6s}'.format(syst,'lnN',''))
    for ich in range(1,n_channels+1):
        for p in pro_names:
            if(p == 'SUSY'):
                if(syst == 'jes'):
                    file_out.write('{:>12s} '.format('JESS'+ str(ich)))
                elif(syst == 'jer'):
                    file_out.write('{:>12s} '.format('JERS'+ str(ich)))
                elif(syst == 'esf'):
                    file_out.write('{:>12s} '.format('ESF'+ str(ich)))
                elif(syst == 'lumi'):
                    file_out.write('{:12.3f} '.format(1.026))
                else: 
                    file_out.write('{:>12s} '.format('-'))
            else:
                syserr = 0
                if(h_rates['h_' + p + '_norm'].GetBinContent(ich)==0):
                    if(p == 'elefakepho' and syst == 'e_to_pho_syst'):
                        syserr = h_rates['h_elefakepho_transferfactor'].GetBinError( ich) + 1
                        file_out.write('{:12.3f} '.format(syserr))
                    elif(p == 'jetfakepho' and syst == 'j_to_pho_syst'):
                        syserr = h_rates['h_jetfakepho_transferfactor'].GetBinError( ich) + 1
                        file_out.write('{:12.3f} '.format(syserr))
                    elif(p == 'qcdfakelep' and syst == 'scale'):
                        syserr = 1 - h_rates['h_qcdfakelep_transferfactor'].GetBinError(ich)
                        file_out.write('{:12.3f} '.format(syserr))
                    else:
                        file_out.write('{:>12s} '.format('-'))
                else:
                    syserr = h_rates['h_' + p + '_syserr_' + syst].GetBinContent(ich)/h_rates['h_' + p + '_norm'].GetBinContent(ich)
                    if(syserr < 0):
                         file_out.write('{:>12s} '.format('-'))
                    else:
                         if(p == "qcdfakelep" and syst == "scale"):
                             syserr = 1 - syserr
                         else: 
                             syserr = syserr + 1
                         file_out.write('{:12.3f} '.format(syserr))
    file_out.write('\n')


file_out.write('{:15s} {:3s} {:6s}'.format('PU','lnN',''))
for ich in range(1,n_channels+1):
    for p in pro_names:
        if(p == 'SUSY'):
            file_out.write('{:>12s} '.format('PUS'))
        else:
            file_out.write('{:>12s} '.format('-'))
file_out.write('\n')

for ich in range(1,n_channels+1):
    file_out.write('{:15s} {:3s} {:6s}'.format('SUSY_stat'+str(ich),'lnN',''))
    for k in range(1,n_channels+1):
        if( k == ich):
            file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} '.format('STSC'+ str(ich),'-', '-', '-', '-', '-'))
        else:
            file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} '.format('-', '-', '-', '-', '-', '-'))
    file_out.write('\n')

for ich in range(1,n_channels+1):
    nevt = h_rates['h_elefakepho_controlsample'].GetBinContent(ich)
    fakerate = h_rates['h_elefakepho_transferfactor'].GetBinContent(ich)
    file_out.write('{:15s} {:3s} {:6d}'.format('e_to_pho_stat'+str(ich),'gmN',int(nevt)))
    for k in range(1,n_channels+1):
        if( k == ich):
            file_out.write('{:>12s} {:12.5f} {:>12s} {:>12s} {:>12s} {:>12s} '.format('-', fakerate, '-', '-', '-', '-'))
            logfile_out.write('{:15} {:2d} {:10.5f}\n'.format('elefakepho', ich, nevt*fakerate))
        else:
            file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} '.format('-', '-', '-', '-', '-', '-'))
    file_out.write('\n')
for ich in range(1,n_channels+1):
    nevt = h_rates['h_jetfakepho_controlsample'].GetBinContent(ich)
    fakerate = h_rates['h_jetfakepho_transferfactor'].GetBinContent(ich)
    file_out.write('{:15s} {:3s} {:6d}'.format('j_to_pho_stat'+str(ich),'gmN',int(nevt)))
    for k in range(1,n_channels+1):
        if( k == ich):
            file_out.write('{:>12s} {:>12s} {:12.5f} {:>12s} {:>12s} {:>12s} '.format('-', '-', fakerate, '-', '-', '-'))
            logfile_out.write('{:15} {:2d} {:10.5f}\n'.format('jetfakepho', ich, nevt*fakerate))
        else:
            file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} '.format('-', '-', '-', '-', '-', '-'))
    file_out.write('\n')
for ich in range(1,n_channels+1):
    nevt = h_rates['h_qcdfakelep_controlsample'].GetBinContent(ich)
    fakerate = h_rates['h_qcdfakelep_transferfactor'].GetBinContent(ich)
    file_out.write('{:15s} {:3s} {:6d}'.format('j_to_lep_stat'+str(ich),'gmN',int(nevt)))
    for k in range(1,n_channels+1):
        if( k == ich):
            file_out.write('{:>12s} {:>12s} {:>12s} {:12.5f} {:>12s} {:>12s} '.format('-', '-', '-', fakerate, '-', '-'))
            logfile_out.write('{:15} {:2d} {:10.5f}\n'.format('jetfakelep', ich, nevt*fakerate))
        else:
            file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} '.format('-', '-', '-', '-', '-', '-'))
    file_out.write('\n')

for ich in range(1,n_channels+1):
    file_out.write('{:15s} {:3s} {:6s}'.format('rare_stat'+str(ich),'lnN',''))
    staterror= 1.0 + h_rates['h_rare_norm'].GetBinError(ich)/h_rates['h_rare_norm'].GetBinContent(ich)
    for k in range(1,n_channels+1):
        if( k == ich):
            file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:12.3f} '.format('-', '-', '-', '-', '-',staterror))
        else:
            file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} '.format('-', '-', '-', '-', '-', '-'))
    file_out.write('\n')

for ich in range(1,n_channels+1):
    file_out.write('{:15s} {:3s} {:6s}'.format('VG_stat'+str(ich),'lnN',''))
    staterror= 1.0 + h_rates['h_VGamma_norm'].GetBinError(ich)/h_rates['h_VGamma_norm'].GetBinContent(ich)
    for k in range(1,n_channels+1):
        if( k == ich):
            file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:12.3f} {:>12s} '.format('-', '-', '-', '-',staterror, '-'))
        else:
            file_out.write('{:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} '.format('-', '-', '-', '-', '-', '-'))
    file_out.write('\n')
file_out.close()
logfile_out.close()

