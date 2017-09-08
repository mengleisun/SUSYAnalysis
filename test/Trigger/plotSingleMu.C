#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TChain.h"
#include "TSystem.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLine.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TLorentzVector.h"
#include "TF1.h"

#include "../../include/analysis_rawData.h"
#include "../../include/analysis_photon.h"
#include "../../include/analysis_muon.h"
#include "../../include/analysis_ele.h"
#include "../../include/analysis_tools.h"

Double_t turnon_func(Double_t *x, Double_t *par)
{
  double halfpoint = par[0];
  double slope1 = par[1];
  double plateau = par[2];
  double constant = par[3];
  double expoffset = par[4];
  //double offset = par[3];
  //double plateau = 1.0;
  double offset = 0;

  double pt = TMath::Max(x[0],0.000001);

   double arg = 0;
   //cout << pt <<", "<< halfpoint <<", " << slope <<endl;
   arg = (pt - halfpoint)/(TMath::Sqrt(pt)*slope1);
   double fitval = (offset+0.5*plateau*(1+TMath::Erf(arg)))*(1-exp(-constant*(pt-expoffset)));
   return fitval;
}

void plotSingleMu(){//main  

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

  TFile *file=TFile::Open("/uscms_data/d3/mengleis/plot_SingleMuonTrigger.root");
  TTree *es = (TTree*)file->Get("mgTree");
  int   mg_run(0);
  float mg_muPt(0);
  float mg_muEta(0);
  bool  isTight;
  float mg_muMiniIso(0);
  float mg_dR(0); 
  int   mg_phofireL1;
  int   mg_mufireL1;
  bool  mg_mgfireHLT;
  bool  passL1;
  bool  passIsoMuHLT;
  bool  passTrkMuHLT;
 
  es->SetBranchAddress("runN",      &mg_run);
  es->SetBranchAddress("muPt",      &mg_muPt);
  es->SetBranchAddress("muEta",     &mg_muEta);
  es->SetBranchAddress("isTight",   &isTight);
  es->SetBranchAddress("muMiniIso", &mg_muMiniIso);
  es->SetBranchAddress("dR",        &mg_dR);
  es->SetBranchAddress("phofireL1", &mg_phofireL1);
  es->SetBranchAddress("mufireL1",  &mg_mufireL1); 
  es->SetBranchAddress("mgfireHLT", &mg_mgfireHLT);
  es->SetBranchAddress("passL1",    &passL1);
  es->SetBranchAddress("passIsoMuHLT",   &passIsoMuHLT);
  es->SetBranchAddress("passTrkMuHLT",   &passTrkMuHLT);

  const unsigned nEvts = es->GetEntries(); 
  std::cout << "total=" << es->GetEntries() << std::endl;
 
  float xaxis[] = {0,10,12,14,16,18,20,22,24,26,28,30,34,38,42,46,50,70,100};
  TProfile *p_Isoeff = new TProfile("p_Isoeff","",100,0,100);
  TProfile *p_eff= new TProfile("p_eff",";#mu p_{T} (GeV)",100,0,100);

  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    es->GetEntry(ievt);
    if(fabs(mg_muEta) > 1.5)continue;
    if(mg_muMiniIso > 0.2)continue;

    if(passL1){
      if(passIsoMuHLT)p_Isoeff->Fill(mg_muPt, 1);
      else p_Isoeff->Fill(mg_muPt, 0);

      if(passIsoMuHLT || passTrkMuHLT)p_eff->Fill(mg_muPt, 1);
      else p_eff->Fill(mg_muPt, 0);
    }
  }

	TLegend *leg=new TLegend(0.5, 0.2, 0.9,0.4);
	leg->AddEntry(p_eff, "IsoMu24 || IsoTkMu24");
	leg->AddEntry(p_Isoeff, "IsoMu24");
  TCanvas *can=new TCanvas("can","",600,600);
  gStyle->SetOptStat(0);
  can->cd();
  p_eff->Draw();
  p_Isoeff->SetLineColor(kRed);
  p_Isoeff->Draw("same");
	leg->Draw("same");
}


