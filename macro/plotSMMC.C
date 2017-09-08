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

void plotSMMC(){//main  
TFile *file[2];

file[0]= TFile::Open("../data/SM_MC_WGToLNuG.root");
file[1]= TFile::Open("../data/T5Wg_mGl-1150to1200.root");

TH1F *p_mcphotonET[2];
TH1F *p_mceleET[2];
TH1F *p_Mt[2];
TH1F *p_PhoELeDeltaR[2];
TH1F *p_PhoEleMass[2];
TH1F *p_selectPhoEleMass[2];

for(unsigned iF(0); iF < 2; iF++){
  p_mcphotonET[iF] = (TH1F*)file[iF]->Get("photonET");
  p_mcphotonET[iF]->Scale(1.0/p_mcphotonET[iF]->GetEntries());
  p_mceleET[iF] = (TH1F*)file[iF]->Get("eleET");
  p_mceleET[iF]->Scale(1.0/p_mceleET[iF]->GetEntries());
  p_mceleET[iF]->SetTitle("signal e (#tilde{#chi}^{#pm}#rightarrowW^{#pm}#tilde{#chi}^{0}, W^{#pm}#rightarrowe) p_{T};p_{T} (GeV)");
  p_Mt[iF] = (TH1F*)file[iF]->Get("Mt");
  p_Mt[iF]->Scale(1.0/p_Mt[iF]->GetEntries());
  p_PhoELeDeltaR[iF] = (TH1F*)file[iF]->Get("p_PhoELeDeltaR");
  p_PhoELeDeltaR[iF]->Scale(1.0/p_PhoELeDeltaR[iF]->GetEntries());
  p_PhoEleMass[iF] = (TH1F*)file[iF]->Get("p_PhoEleMass");
  p_PhoEleMass[iF]->Scale(1.0/p_PhoEleMass[iF]->GetEntries());
  p_selectPhoEleMass[iF] = (TH1F*)file[iF]->Get("p_selectPhoEleMass");
  p_selectPhoEleMass[iF]->Scale(1.0/p_selectPhoEleMass[iF]->GetEntries());
}

setTDRStyle();    
TCanvas *hist_MC[3];
std::ostringstream canvas; 

for(unsigned iF(0); iF<3; iF++){
    canvas.str("");
    canvas << "canvas" << iF;
    hist_MC[iF] = new TCanvas(canvas.str().c_str(),canvas.str().c_str(),600,600);
}

TCanvas *canEt = new TCanvas("Photon Et","Photon Et",600,600);
canEt->cd();
gPad->SetLogy();
int LineColor[] = {2,3,4,6,41};
string FileList[] = {"WGToLNuG","T5WG(M_{#tilde{g}}=1150GeV)"};
TLegend *legEt =  new TLegend(0.45,0.8,0.9,0.90);
std::ostringstream histname;
for(unsigned ii(0); ii<2; ii++){
  histname.str("");
  histname << FileList[ii];
  p_mcphotonET[ii]->SetLineColor(LineColor[ii]);
  if(ii>0)p_mcphotonET[ii]->SetFillColor(LineColor[ii]);
  p_mcphotonET[ii]->SetFillStyle(3444);
  legEt->AddEntry(p_mcphotonET[ii],histname.str().c_str());
  if(ii==0)p_mcphotonET[ii]->Draw();
  else p_mcphotonET[ii]->Draw("same");
}
legEt->Draw("same");

TCanvas *canPt = new TCanvas("Electron Pt","Electron Pt",600,600);
canPt->cd();
gPad->SetLogy();
for(unsigned ii(0); ii<2; ii++){
  p_mceleET[ii]->SetLineColor(LineColor[ii]);
  if(ii>0)p_mceleET[ii]->SetFillColor(LineColor[ii]);
  p_mceleET[1]->SetFillStyle(3444);
  p_mceleET[ii]->GetXaxis()->SetRangeUser(25,400);
  if(ii==0)p_mceleET[ii]->Draw();
  else p_mceleET[ii]->Draw("same");
}
legEt->Draw("same");


TCanvas *canMT = new TCanvas("MT","MT",600,600);
canMT->cd();
gPad->SetLogy();
for(unsigned ii(0); ii<2; ii++){
  p_Mt[ii]->SetLineColor(LineColor[ii]);
  if(ii>0)p_Mt[ii]->SetFillColor(LineColor[ii]);
   p_Mt[1]->SetFillStyle(3444);
  if(ii==0)p_Mt[ii]->Draw();
  else p_Mt[ii]->Draw("same");
}
legEt->Draw("same");

  hist_MC[1]->cd();
gPad->SetLogy();
  p_PhoELeDeltaR[0]->GetXaxis()->SetRangeUser(0,6);
  p_PhoELeDeltaR[0]->SetLineColor(LineColor[0]);
  p_PhoELeDeltaR[0]->Draw();
  p_PhoELeDeltaR[1]->SetLineColor(LineColor[1]);
  p_PhoELeDeltaR[1]->SetFillColor(LineColor[1]);
  p_PhoELeDeltaR[1]->SetFillStyle(3444);
  p_PhoELeDeltaR[1]->Draw("same");
  legEt->Draw("same");
  hist_MC[2]->cd();
gPad->SetLogy();
  p_PhoEleMass[0]->SetLineColor(LineColor[0]);
  p_PhoEleMass[0]->Draw();
  p_PhoEleMass[1]->SetLineColor(LineColor[1]);
  p_PhoEleMass[1]->SetFillColor(LineColor[1]);
  p_PhoEleMass[1]->SetFillStyle(3444);
  p_PhoEleMass[1]->Draw("same");
  legEt->Draw("same");


for(unsigned iF(0); iF<2; iF++){
    canvas.str("");
    canvas << "canvas" << iF << ".pdf";
    hist_MC[iF]->SaveAs(canvas.str().c_str());
    canvas.str("");
    canvas << "canvas" << iF << ".png";
    hist_MC[iF]->SaveAs(canvas.str().c_str());
}

canEt->SaveAs("SMMC_photonEt.pdf");
canPt->SaveAs("SMMC_elePt.pdf");
canMT->SaveAs("SMMC_MT.pdf");
canEt->SaveAs("SMMC_photonEt.png");
canPt->SaveAs("SMMC_elePt.png");
canMT->SaveAs("SMMC_MT.png");
}


