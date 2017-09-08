#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include "TH1F.h"
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

void plotmg(){//main  

TFile *file = TFile::Open("../data/plot_mg.root");

setTDRStyle();    
TCanvas *hist_MC[18];
std::ostringstream canvas; 

for(unsigned iC(0); iC<18; iC++){
   canvas.str("");
   canvas << "canvas" << iC;
   hist_MC[iC] = new TCanvas(canvas.str().c_str(),canvas.str().c_str(),600,600);
}

  TH1F *p_photonET = (TH1F*)file->Get("photonET");
  TH1F *p_muET = (TH1F*)file->Get("eleET");
  TH1F *p_Mt = (TH1F*)file->Get("Mt");
  TH1F *p_MET = (TH1F*)file->Get("MET");

  TH1F *p_phoEBR9_Nm1 = (TH1F*)file->Get("p_phoEBR9_Nm1");
  TH1F *p_phoEBHoverE_Nm1 = (TH1F*)file->Get("p_phoEBHoverE_Nm1");
  TH1F *p_phoEBsigma_Nm1 = (TH1F*)file->Get("p_phoEBsigma_Nm1");
  TH1F *p_phoEBChIso_Nm1 = (TH1F*)file->Get("p_phoEBChIso_Nm1");
  TH1F *p_phoEBNeuIso_Nm1 = (TH1F*)file->Get("p_phoEBNeuIso_Nm1");
  TH1F *p_phoEBPhoIso_Nm1 = (TH1F*)file->Get("p_phoEBPhoIso_Nm1");
  TH1F *p_phoEER9_Nm1 = (TH1F*)file->Get("p_phoEER9_Nm1");
  TH1F *p_phoEEHoverE_Nm1 = (TH1F*)file->Get("p_phoEEHoverE_Nm1");
  TH1F *p_phoEEsigma_Nm1 = (TH1F*)file->Get("p_phoEEsigma_Nm1");
  TH1F *p_phoEEChIso_Nm1 = (TH1F*)file->Get("p_phoEEChIso_Nm1");
  TH1F *p_phoEENeuIso_Nm1 = (TH1F*)file->Get("p_phoEENeuIso_Nm1");
  TH1F *p_phoEEPhoIso_Nm1 = (TH1F*)file->Get("p_phoEEPhoIso_Nm1");
  TH1F *p_PhoELeDeltaR = (TH1F*)file->Get("p_PhoELeDeltaR");
  TH1F *p_PhoEleMass = (TH1F*)file->Get("p_PhoEleMass");

hist_MC[0]->cd();
p_photonET->GetXaxis()->SetRangeUser(35,200);
p_photonET->Draw();
hist_MC[1]->cd();
p_muET->GetXaxis()->SetRangeUser(25,100);
p_muET->Draw();
hist_MC[2]->cd();
p_Mt->GetXaxis()->SetRangeUser(0,200);
p_Mt->Draw();
hist_MC[3]->cd();
p_MET->Draw();
hist_MC[4]->cd();
p_phoEBR9_Nm1->Draw();
TLine *EBR9cut = new TLine(0.5,0,0.5,8000);
EBR9cut->Draw("same");
hist_MC[5]->cd();
gPad->SetLogy();
p_phoEBHoverE_Nm1->Draw();
TLine *EBHEcut = new TLine(0.05,0,0.05,10000);
EBHEcut->Draw("same");
hist_MC[6]->cd();
p_phoEBsigma_Nm1->Draw();
TLine *EBSigmacut = new TLine(0.0102,0,0.0102,12000);
EBSigmacut->Draw("same");
hist_MC[7]->cd();
gPad->SetLogy();
p_phoEBChIso_Nm1->Draw();
TLine *EBChIsocut = new TLine(3.32,0,3.32,10000);
EBChIsocut->Draw("same");
hist_MC[8]->cd();
gPad->SetLogy();
p_phoEBNeuIso_Nm1->Draw();
hist_MC[9]->cd();
gPad->SetLogy();
p_phoEBPhoIso_Nm1->Draw();
hist_MC[10]->cd();
p_phoEER9_Nm1->Draw();
TLine *EER9cut = new TLine(0.8,0,0.8,8000);
EER9cut->Draw("same");
hist_MC[11]->cd();
gPad->SetLogy();
p_phoEEHoverE_Nm1->Draw();
TLine *EEHEcut = new TLine(0.05,0,0.05,10000);
EEHEcut->Draw("same");
hist_MC[12]->cd();
p_phoEEsigma_Nm1->Draw();
TLine *EESigmacut = new TLine(0.0274,0,0.0274,12000);
EESigmacut->Draw("same");
hist_MC[13]->cd();
gPad->SetLogy();
p_phoEEChIso_Nm1->Draw();
TLine *EEChIsocut = new TLine(1.97,0,1.97,10000);
EEChIsocut->Draw("same");
hist_MC[14]->cd();
gPad->SetLogy();
p_phoEENeuIso_Nm1->Draw();
hist_MC[15]->cd();
gPad->SetLogy();
p_phoEEPhoIso_Nm1->Draw();
hist_MC[16]->cd();
p_PhoELeDeltaR->GetXaxis()->SetRangeUser(0,5);
p_PhoELeDeltaR->Draw();
hist_MC[17]->cd();
p_PhoEleMass->Draw();


for(unsigned iC(0); iC<18; iC++){
   canvas.str("");
   canvas << "mg_plot" << iC << ".png";
   hist_MC[iC]->SaveAs(canvas.str().c_str());
}
}


