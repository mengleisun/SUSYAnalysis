#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TF1.h"
#include "TMath.h"
#include "../include/tdrstyle.C"

void plotMC(){//main  
TFile *file[6];

file[0]= TFile::Open("../data/plot_T5Wg_mGl-800to1000.root");  
file[1]= TFile::Open("../data/plot_T5Wg_mGl-1050to1100.root");
file[2]= TFile::Open("../data/plot_T5Wg_mGl-1250to1300.root");
file[3]= TFile::Open("../data/plot_T5Wg_mGl-1350to1400.root");
file[4]= TFile::Open("../data/plot_T5Wg_mGl-1450to1500.root");
file[5]= TFile::Open("../data/plot_T5Wg_mGl-1550.root");

TH2F *p_SUSYMass[6];
TH1F *p_mcphotonET[6];
TH1F *p_mceleET[6];
TH1F *p_phoEBR9_true[6]; 
TH1F *p_phoEER9_true[6];
TH1F *p_phoEBHoverE_true[6];
TH1F *p_phoEEHoverE_true[6];
TH1F *p_phoEBsigma_true[6];
TH1F *p_phoEEsigma_true[6];
TH1F *p_phoEBChIso_true[6];
TH1F *p_phoEEChIso_true[6];
TH1F *p_phoEBNeuIso_true[6];
TH1F *p_phoEENeuIso_true[6];
TH1F *p_phoEBPhoIso_true[6];
TH1F *p_phoEEPhoIso_true[6];

TH1F *p_eleEBR9_true[6];
TH1F *p_eleEER9_true[6];

TH1F *p_eleEBR9_idselect[6];
TH1F *p_eleEER9_idselect[6];

TH1F *p_Mt[6];
TH1F *p_PhoELeDeltaR[6];
TH1F *p_PhoEleMass[6];
TH1F *p_selectPhoEleMass[6];

for(unsigned iF(0); iF < 6; iF++){
  p_SUSYMass[iF] = (TH2F*)file[iF]->Get("Mass");

  p_mcphotonET[iF] = (TH1F*)file[iF]->Get("photonET");
  p_mceleET[iF] = (TH1F*)file[iF]->Get("eleET");

  p_phoEBR9_true[iF] = (TH1F*)file[iF]->Get("p_phoEBR9");
  p_phoEBR9_true[iF]->GetXaxis()->SetRangeUser(0.5,1);
  p_phoEER9_true[iF] = (TH1F*)file[iF]->Get("p_phoEER9");
  p_phoEER9_true[iF]->GetXaxis()->SetRangeUser(0.5,1);
  p_phoEBHoverE_true[iF] = (TH1F*)file[iF]->Get("p_phoEBHoverE");
  p_phoEEHoverE_true[iF] = (TH1F*)file[iF]->Get("p_phoEEHoverE");
  p_phoEBsigma_true[iF] = (TH1F*)file[iF]->Get("p_phoEBsigma");
  p_phoEEsigma_true[iF] = (TH1F*)file[iF]->Get("p_phoEEsigma");
  p_phoEBChIso_true[iF] = (TH1F*)file[iF]->Get("p_phoEBChIso");
  p_phoEEChIso_true[iF] = (TH1F*)file[iF]->Get("p_phoEEChIso");
  p_phoEBNeuIso_true[iF] = (TH1F*)file[iF]->Get("p_phoEBNeuIso");
  p_phoEENeuIso_true[iF] = (TH1F*)file[iF]->Get("p_phoEENeuIso");
  p_phoEBPhoIso_true[iF] = (TH1F*)file[iF]->Get("p_phoEBPhoIso");
  p_phoEEPhoIso_true[iF] = (TH1F*)file[iF]->Get("p_phoEEPhoIso");

  p_eleEBR9_true[iF] = (TH1F*)file[iF]->Get("p_eleEBR9");
  p_eleEBR9_true[iF]->GetXaxis()->SetRangeUser(0.5,1);
  p_eleEER9_true[iF] = (TH1F*)file[iF]->Get("p_eleEER9");
  p_eleEER9_true[iF]->GetXaxis()->SetRangeUser(0.5,1);

  p_eleEBR9_idselect[iF] = (TH1F*)file[iF]->Get("p_eleEBR9_idselect");
  p_eleEBR9_idselect[iF]->GetXaxis()->SetRangeUser(0.5,1);
  p_eleEER9_idselect[iF] = (TH1F*)file[iF]->Get("p_eleEER9_idselect");
  p_eleEER9_idselect[iF]->GetXaxis()->SetRangeUser(0.5,1);

  p_Mt[iF] = (TH1F*)file[iF]->Get("Mt");
  p_PhoELeDeltaR[iF] = (TH1F*)file[iF]->Get("p_PhoELeDeltaR");
  p_PhoEleMass[iF] = (TH1F*)file[iF]->Get("p_PhoEleMass");
  p_selectPhoEleMass[iF] = (TH1F*)file[iF]->Get("p_selectPhoEleMass");
}

p_SUSYMass[0]->Add(p_SUSYMass[1]); 
p_SUSYMass[0]->Add(p_SUSYMass[2]);
p_SUSYMass[0]->Add(p_SUSYMass[3]);
p_SUSYMass[0]->Add(p_SUSYMass[4]);
p_SUSYMass[0]->Add(p_SUSYMass[5]);

setTDRStyle();    
TCanvas *hist_MC[6][19];
std::ostringstream canvas; 

for(unsigned iF(0); iF<6; iF++){
  for(unsigned iC(0); iC<19; iC++){
    canvas.str("");
    canvas << "energy" << iF << "_canvas" << iC;
    hist_MC[iF][iC] = new TCanvas(canvas.str().c_str(),canvas.str().c_str(),600,600);
  }
}

TCanvas *canMass = new TCanvas("SUSY Mass","SUSY Mass",600,600);
TCanvas *canEt = new TCanvas("Photon Et","Photon Et",600,600);
canEt->cd();
int norm(1);
int LineColor[] = {1,2,3,4,6,41};
string FileList[] = {"M_{#tilde{g}}=800 to 1000","M_{#tilde{g}}=1050 to 1100", "M_{#tilde{g}}=1250 to 1300","M_{#tilde{g}}=1350 to 1400","M_{#tilde{g}}=1450 to 1500","M_{#tilde{g}}=1550"};
TLegend *legEt =  new TLegend(0.4,0.55,0.86,0.85);
std::ostringstream histname;
for(unsigned ii(0); ii<6; ii++){
  histname.str("");
  histname << FileList[ii];
  norm = p_mcphotonET[ii]->GetEntries();
  p_mcphotonET[ii]->Scale(1.0/norm);
  p_mcphotonET[ii]->SetLineColor(LineColor[ii]);
  legEt->AddEntry(p_mcphotonET[ii],histname.str().c_str());
  if(ii==0)p_mcphotonET[ii]->Draw();
  else p_mcphotonET[ii]->Draw("same");
}
legEt->Draw("same");

canMass->cd();
Int_t PaletteColors[] = {kBlue+1, kBlue, kBlue-7, kCyan,kGreen, kGreen-6}; // #colors >= #levels - 1
TLine *MassGridX[32];
TLine *MassGridY[16];
for(unsigned ix(0); ix < 32; ix++)MassGridX[ix]= new TLine(ix*50,800,ix*50,1600);
for(unsigned iy(0); iy < 16; iy++)MassGridY[iy]= new TLine(0,iy*50+800,1600,iy*50+800);
gStyle->SetPalette((sizeof(PaletteColors)/sizeof(Int_t)), PaletteColors);
p_SUSYMass[0]->GetYaxis()->SetTitleOffset(1.7);
p_SUSYMass[0]->Draw("colz");
for(unsigned ix(0); ix < 32; ix++)MassGridX[ix]->Draw("same");
for(unsigned iy(0); iy < 16; iy++)MassGridY[iy]->Draw("same");

TCanvas *canPt = new TCanvas("Electron Pt","Electron Pt",600,600);
canPt->cd();
for(unsigned ii(0); ii<6; ii++){
  norm = p_mceleET[ii]->GetEntries();
  p_mceleET[ii]->Scale(1.0/norm);
  p_mceleET[ii]->SetLineColor(LineColor[ii]);
  if(ii==0)p_mceleET[ii]->Draw();
  else p_mceleET[ii]->Draw("same");
}
legEt->Draw("same");


TCanvas *canMT = new TCanvas("MT","MT",600,600);
canMT->cd();
for(unsigned ii(0); ii<6; ii++){
  norm = p_Mt[ii]->GetEntries();
  p_Mt[ii]->Scale(1.0/norm);
  p_Mt[ii]->SetLineColor(LineColor[ii]);
  if(ii==0)p_Mt[ii]->Draw();
  else p_Mt[ii]->Draw("same");
}
legEt->Draw("same");


for(unsigned iF(0); iF<6; iF++){
  TAxis *axis = p_selectPhoEleMass[iF]->GetXaxis();
  Int_t bmin = axis->FindBin(0.0); 
  Int_t bmax = axis->FindBin(95.0); 
  double integral = p_selectPhoEleMass[iF]->Integral(bmin,bmax);
  std::cout << integral << " events less then 95GeV out of " << p_selectPhoEleMass[iF]->GetEntries() << std::endl;

  hist_MC[iF][0]->cd();
  p_selectPhoEleMass[iF]->Draw();
  TLine *masscut = new TLine(95,0,95,1000);
  masscut->SetLineColor(kRed);
  masscut->Draw("same");

  hist_MC[iF][1]->cd();
  p_mceleET[iF]->Draw();
  hist_MC[iF][2]->cd();
  p_phoEBR9_true[iF]->Draw();
  hist_MC[iF][3]->cd();
  p_phoEER9_true[iF]->Draw();
  hist_MC[iF][4]->cd();
  gPad->SetLogy();
  p_phoEBHoverE_true[iF]->Draw();
  hist_MC[iF][5]->cd();
  gPad->SetLogy();
  p_phoEEHoverE_true[iF]->Draw();
  hist_MC[iF][6]->cd();
  p_phoEBsigma_true[iF]->Draw();
  hist_MC[iF][7]->cd();
  p_phoEEsigma_true[iF]->Draw();
  hist_MC[iF][8]->cd();
  gPad->SetLogy();
  p_phoEBChIso_true[iF]->Draw();
  hist_MC[iF][9]->cd();
  gPad->SetLogy();
  p_phoEEChIso_true[iF]->Draw();
  hist_MC[iF][10]->cd();
  gPad->SetLogy();
  p_phoEBNeuIso_true[iF]->Draw();
  hist_MC[iF][11]->cd();
  gPad->SetLogy();
  p_phoEENeuIso_true[iF]->Draw();
  hist_MC[iF][12]->cd();
  gPad->SetLogy();
  p_phoEBPhoIso_true[iF]->Draw();
  hist_MC[iF][13]->cd();
  gPad->SetLogy();
  p_phoEEPhoIso_true[iF]->Draw();
  hist_MC[iF][14]->cd();
  p_eleEBR9_true[iF]->SetLineColor(kBlue);
  p_eleEBR9_idselect[iF]->SetLineColor(kRed);
  TLegend *legEBR9 =  new TLegend(0.4,0.65,0.86,0.85);
  legEBR9->AddEntry(p_eleEBR9_true[iF],"all e#pm R9");
  legEBR9->AddEntry(p_eleEBR9_idselect[iF],"pt>25GeV,medium WP e#pm,R9"); 
  p_eleEBR9_true[iF]->Draw();
  p_eleEBR9_idselect[iF]->Draw("same");
  legEBR9->Draw("same");

  hist_MC[iF][15]->cd();
  p_eleEER9_true[iF]->SetLineColor(kBlue);
  p_eleEER9_idselect[iF]->SetLineColor(kRed);
  TLegend *legEER9 =  new TLegend(0.4,0.65,0.86,0.85);
  legEER9->AddEntry(p_eleEER9_true[iF],"all e#pm R9");
  legEER9->AddEntry(p_eleEER9_idselect[iF],"pt>25GeV,medium WP e#pm, R9"); 
  p_eleEER9_true[iF]->Draw();
  p_eleEER9_idselect[iF]->Draw("same");
  legEER9->Draw("same");

  hist_MC[iF][16]->cd();
  p_Mt[iF]->Draw();
  hist_MC[iF][17]->cd();
  p_PhoELeDeltaR[iF]->Draw();
  hist_MC[iF][18]->cd();
  p_PhoEleMass[iF]->Draw();
}


for(unsigned iF(0); iF<6; iF++){
  for(unsigned iC(0); iC<19; iC++){
    canvas.str("");
    canvas << "energy" << iF << "_canvas" << iC << ".pdf";
    hist_MC[iF][iC]->SaveAs(canvas.str().c_str());
  }
}

canEt->SaveAs("MC_photonEt.pdf");
canMass->SaveAs("SUSY_Mass.pdf");
canPt->SaveAs("MC_elePt.pdf");
canMT->SaveAs("MC_MT.pdf");
}


