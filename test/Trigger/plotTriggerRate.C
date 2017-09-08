#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include "TH1F.h"
#include "TGraph.h"
#include "TCanvas.h"

void plotTriggerRate(){//main 
 
Double_t lumis[1160];
Double_t rates1[1160];
Double_t rates2[1160];

int nLumi(0);
float rate_Diphoton(0), rate_MuonEG(0), instLumi(0);
std::ifstream fRate("HLTRate_run283865.txt");
if(fRate.is_open()){
  for(int i(15); i< 1175; i++){
    fRate >> nLumi >> rate_Diphoton >> rate_MuonEG >> instLumi;
    lumis[i-15] = instLumi*1e30;
    rates1[i-15] = rate_Diphoton;
    rates2[i-15] = rate_MuonEG;
  }
}

TCanvas *can1=new TCanvas("can1","",600,600);
can1->cd();
TGraph *Diphoton = new TGraph(1175-15, lumis, rates1);
Diphoton->SetMarkerColor(kBlue);
Diphoton->SetTitle("Rate  of HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90");
Diphoton->GetXaxis()->SetTitle("Lumi");
Diphoton->GetYaxis()->SetTitle("rate(Hz)");
Diphoton->Fit("pol1");
Diphoton->Draw("A*");

TCanvas *can2=new TCanvas("can2","",600,600);
can2->cd();
TGraph *MuEG = new TGraph(1175-15, lumis, rates2);
MuEG->SetMarkerColor(kBlue);
MuEG->SetTitle("Rate  of HLT_Mu17_Photon30_CaloIdL_*");
MuEG->GetXaxis()->SetTitle("Lumi");
MuEG->GetYaxis()->SetTitle("rate(Hz)");
MuEG->Fit("pol1");
MuEG->Draw("A*");
}
