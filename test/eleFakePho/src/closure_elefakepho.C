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
#include "TLorentzVector.h"

#include "../../../include/analysis_tools.h"

double num[]={1094.43 ,1643.62 ,2421.58 ,3682.88 ,816.148 ,359.317 ,243.047 ,48.3986};
double den[]={86611 ,145658 ,212125 ,428289 ,129435 ,50437 ,39536.2 ,23002.2};

void closure_elefakepho(){//main  
	 
  float DYJet_XS = 5943.2;
  float DYJet_nevt = 28747969;
  float TT_XS = 670.3; 
  float TT_nevt = 37459078;
  float WW_XS = 63.21;
  float WW_nevt = 993214; 
  
  TChain *proxytree = new TChain("signalTree");
  proxytree->Add("/uscms_data/d3/mengleis/SUSYAnalysisData/data/proxyTree_egMC.root");

  const unsigned nEvts = proxytree->GetEntries(); 

  std::ostringstream outputname;
  outputname << "closure_DYTT.root";
  TFile *outputfile = TFile::Open(outputname.str().c_str(),"RECREATE");
  outputfile->cd();

//************ Signal Tree **********************//
  float proxyphoEt(0);
  float proxyphoEta(0);
  float proxyphoPhi(0);
  float proxylepPt(0);
  float proxylepEta(0);
  float proxylepPhi(0);
  float proxysigMT(0);
  float proxysigMET(0);
  float proxysigMETPhi(0);
  float proxydPhiLepMET(0);
  int   proxynVertex(0);
  float proxydRPhoLep(0);
  float proxyHT(0);
  float proxynJet(0);
  
  proxytree->SetBranchAddress("phoEt",     &proxyphoEt);
  proxytree->SetBranchAddress("phoEta",    &proxyphoEta);
  proxytree->SetBranchAddress("phoPhi",    &proxyphoPhi);
  proxytree->SetBranchAddress("lepPt",     &proxylepPt);
  proxytree->SetBranchAddress("lepEta",    &proxylepEta);
  proxytree->SetBranchAddress("lepPhi",    &proxylepPhi);
  proxytree->SetBranchAddress("sigMT",     &proxysigMT);
  proxytree->SetBranchAddress("sigMET",    &proxysigMET);
  proxytree->SetBranchAddress("sigMETPhi", &proxysigMETPhi);
  proxytree->SetBranchAddress("dPhiLepMET",&proxydPhiLepMET);
  proxytree->SetBranchAddress("nVertex",   &proxynVertex);
  proxytree->SetBranchAddress("dRPhoLep",  &proxydRPhoLep);
  proxytree->SetBranchAddress("HT",        &proxyHT);
  proxytree->SetBranchAddress("nJet",      &proxynJet);
  
  TChain *proxymctree = new TChain("MCTree");
  proxymctree->Add("/uscms_data/d3/mengleis/SUSYAnalysisData/data/proxyTree_egMC.root");
  int proxymcType = 0;
  proxymctree->SetBranchAddress("mcType",          &proxymcType);
//*********** histo list **********************//
  TH1F *p_allPhoEt = new TH1F("p_allPhoEt","#gamma E_{T}; E_{T} (GeV)",500,0,1500);
  TH1F *p_allPhoEta = new TH1F("p_allPhoEta","#gamma #eta; #eta;",60,-3,3);
  TH1F *p_allMt = new TH1F("p_allMt","M_{T}; M_{T} (GeV);",200,0,400); 
  TH1F *p_alldPhiEleMET = new TH1F("p_alldPhiEleMET","dPhiEleMET",70,-3.5,3.5); 

  TH1F *p_elebkgPhoEt = new TH1F("p_elebkgPhoEt","#gamma E_{T}; E_{T} (GeV)",500,0,1500);
  TH1F *p_elebkgPhoEta = new TH1F("p_elebkgPhoEta","#gamma #eta; #eta;",60,-3,3);
  TH1F *p_elebkgMt = new TH1F("p_elebkgMt","M_{T}; M_{T} (GeV);",200,0,400); 
  TH1F *p_eledPhiEleMET = new TH1F("p_eledPhiEleMET","dPhiEleMET",70,-3.5,3.5); 
 
  TH1F *p_elebkgPhoEtTT = new TH1F("p_elebkgPhoEtTT","#gamma E_{T}; E_{T} (GeV)",500,0,1500);
  TH1F *p_elebkgPhoEtaTT = new TH1F("p_elebkgPhoEtaTT","#gamma #eta; #eta;",60,-3,3);
  TH1F *p_elebkgMtTT = new TH1F("p_elebkgMtTT","M_{T}; M_{T} (GeV);",200,0,400); 
  TH1F *p_eledPhiEleMETTT = new TH1F("p_eledPhiEleMETTT","dPhiEleMET",70,-3.5,3.5); 

  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
	proxytree->GetEntry(ievt);
    proxymctree->GetEntry(ievt);

    double eventweight = 1;
    if(proxymcType==1)eventweight = 100000.0*DYJet_XS/DYJet_nevt;
    else if(proxymcType==2)eventweight = 100000.0*TT_XS/TT_nevt;
    else if(proxymcType==4)eventweight = 100000.0*WW_XS/WW_nevt;

    double fakerate = pow(0.0293*proxyphoEt+3.37,-3.04);
    if(proxyphoEt>110.0)fakerate = 0.5*(fakerate + (243.047+48.3986)/(243.047+48.3986+39536.2+23002.2)); 
	double w_ele = fakerate/(1-fakerate)*eventweight;
	float deltaPhi = proxylepPhi - proxysigMETPhi;
	if(fabs(deltaPhi) > TMath::Pi()){
	  if(deltaPhi > 0)deltaPhi = -1.0*(TMath::TwoPi() - fabs(deltaPhi));
	  else deltaPhi = TMath::TwoPi() - fabs(deltaPhi);
	}

    if(proxymcType==1){
	  p_elebkgPhoEt->Fill(proxyphoEt,w_ele);
	  p_elebkgPhoEta->Fill(proxyphoEta, w_ele);
	  p_elebkgMt->Fill(proxysigMT, w_ele);
	  p_eledPhiEleMET->Fill(deltaPhi, w_ele);
    }
    else if(proxymcType==2 || proxymcType==4){
	  p_elebkgPhoEtTT->Fill(proxyphoEt,w_ele);
	  p_elebkgPhoEtaTT->Fill(proxyphoEta, w_ele);
	  p_elebkgMtTT->Fill(proxysigMT, w_ele);
	  p_eledPhiEleMETTT->Fill(deltaPhi, w_ele);
    }
  }        
 
//************ Signal Tree **********************//
  TChain *sigtree = new TChain("signalTree");
  sigtree->Add("/uscms_data/d3/mengleis/SUSYAnalysisData/data/resTree_egMC.root");

  float sigphoEt(0);
  float sigphoEta(0);
  float sigphoPhi(0);
  float siglepPt(0);
  float siglepEta(0);
  float siglepPhi(0);
  float sigsigMT(0);
  float sigsigMET(0);
  float sigsigMETPhi(0);
  float sigdPhiLepMET(0);
  int   signVertex(0);
  float sigdRPhoLep(0);
  float sigHT(0);
  float signJet(0);
  
  sigtree->SetBranchAddress("phoEt",     &sigphoEt);
  sigtree->SetBranchAddress("phoEta",    &sigphoEta);
  sigtree->SetBranchAddress("phoPhi",    &sigphoPhi);
  sigtree->SetBranchAddress("lepPt",     &siglepPt);
  sigtree->SetBranchAddress("lepEta",    &siglepEta);
  sigtree->SetBranchAddress("lepPhi",    &siglepPhi);
  sigtree->SetBranchAddress("sigMT",     &sigsigMT);
  sigtree->SetBranchAddress("sigMET",    &sigsigMET);
  sigtree->SetBranchAddress("sigMETPhi", &sigsigMETPhi);
  sigtree->SetBranchAddress("dPhiLepMET",&sigdPhiLepMET);
  sigtree->SetBranchAddress("nVertex",   &signVertex);
  sigtree->SetBranchAddress("dRPhoLep",  &sigdRPhoLep);
  sigtree->SetBranchAddress("HT",        &sigHT);
  sigtree->SetBranchAddress("nJet",      &signJet);

  TChain *mctree = new TChain("MCTree");
  mctree->Add("/uscms_data/d3/mengleis/SUSYAnalysisData/data/resTree_egMC.root");
  int mcType = 0;
  std::vector<int> *mcPID=0;
  std::vector<float> *mcEta=0;
  std::vector<float> *mcPhi=0;
  std::vector<float> *mcPt=0;
  std::vector<int> *mcMomPID=0;
  mctree->SetBranchAddress("mcType",          &mcType);
  mctree->SetBranchAddress("mcPID",           &mcPID);
  mctree->SetBranchAddress("mcEta",           &mcEta);
  mctree->SetBranchAddress("mcPhi",           &mcPhi);
  mctree->SetBranchAddress("mcPt",            &mcPt);
  mctree->SetBranchAddress("mcMomPID",        &mcMomPID);

  TH1F *p_dRMC = new TH1F("p_dRMC","p_dRMC",30,0,0.3);
  TH1F *p_dEMC = new TH1F("p_dEMC","p_dEMC",100,0,1);
  TH1F *p_nMatch = new TH1F("p_nMatch","p_nMatch",10,0,10);
  for (unsigned ievt(0); ievt<sigtree->GetEntries(); ++ievt){//loop on entries
	sigtree->GetEntry(ievt);
    mctree->GetEntry(ievt);	

    double eventweight = 1;
    if(mcType==1)eventweight = 100000.0*DYJet_XS/DYJet_nevt;
    else if(mcType==2)eventweight = 100000.0*TT_XS/TT_nevt; 
    else if(mcType==4)eventweight = 100000.0*WW_XS/WW_nevt;

	double mindR(0.3),deltaE(1);
	unsigned matchIndex(0);
    int nMatch(0);
    int nMatchPho(0);
	for(unsigned iMC(0); iMC < mcPID->size(); iMC++){
      if((*mcPt)[iMC] < 10)continue;
	  double dR = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], sigphoEta, sigphoPhi);
	  double dE = fabs((*mcPt)[iMC] - sigphoEt)/sigphoEt;
	  if(dR < mindR && dE < 0.5){mindR=dR; matchIndex=iMC;deltaE = dE;}
	}
	if(mindR < 0.15 && deltaE < 0.5){
	  if(((*mcPID)[matchIndex] == 11 || (*mcPID)[matchIndex] == -11)){
        float deltaPhi = siglepPt - sigsigMETPhi;
		if(fabs(deltaPhi) > TMath::Pi()){
		  if(deltaPhi > 0)deltaPhi = -1.0*(TMath::TwoPi() - fabs(deltaPhi));
		  else deltaPhi = TMath::TwoPi() - fabs(deltaPhi);
		}

		p_allPhoEt->Fill(sigphoEt,eventweight);
		p_allPhoEta->Fill(sigphoEta, eventweight);
		p_allMt->Fill(sigsigMT, eventweight);
		p_alldPhiEleMET->Fill(deltaPhi, eventweight);
      }
    }

  }        
 
outputfile->Write();

gStyle->SetOptStat(0);
Double_t plotPtBins[]={35,39,42,45,48,51,54,57,60,63,66,69,72,75,78,81,90,100,110,120,130,140,200};
TCanvas *c_pt = new TCanvas("Photon_Pt", "Photon P_{T}",600,600);
c_pt->cd();
gPad->SetLogy();
TH1 *new_allPhoEt = p_allPhoEt->Rebin(2);
TH1 *new_elebkgPhoEt = p_elebkgPhoEt->Rebin(2);
TH1 *new_elebkgPhoEtTT = p_elebkgPhoEtTT->Rebin(2);
new_allPhoEt->GetXaxis()->SetRangeUser(35,200);
new_allPhoEt->Draw();
new_allPhoEt->SetLineColor(kBlack);
new_allPhoEt->SetMarkerStyle(20);
new_elebkgPhoEt->SetFillStyle(1001);
new_elebkgPhoEt->SetLineColor(kRed);
new_elebkgPhoEt->SetFillColor(kRed);
new_elebkgPhoEt->Draw("hist same");
new_elebkgPhoEtTT->Add(new_elebkgPhoEt);
new_elebkgPhoEtTT->SetLineColor(kGreen);
new_elebkgPhoEtTT->SetFillStyle(1001);
new_elebkgPhoEtTT->SetFillColor(kGreen);
new_elebkgPhoEtTT->Draw("hist same");
new_elebkgPhoEt->Draw("hist same");
TLegend *leg =  new TLegend(0.6,0.7,0.9,0.9);
leg->SetFillStyle(0);
gStyle->SetLegendBorderSize(1);
gStyle->SetLegendFillColor(0);
leg->AddEntry(new_allPhoEt,"observed");
leg->AddEntry(new_elebkgPhoEtTT,"EWK (tt, WW)");
leg->AddEntry(new_elebkgPhoEt,"DY");
leg->Draw("same");
new_allPhoEt->Draw("same");
c_pt->SaveAs("EfakePho_closure.pdf");

TCanvas *c_mt = new TCanvas("MT", "MT",600,600);
c_mt->cd();
gPad->SetLogy();
TH1 *new_allMt = p_allMt->Rebin(2);
TH1 *new_elebkgMt = p_elebkgMt->Rebin(2);
TH1 *new_elebkgMtTT = p_elebkgMtTT->Rebin(2);
new_allMt->GetXaxis()->SetRangeUser(35,200);
new_allMt->Draw();
new_allMt->SetLineColor(kBlack);
new_allMt->SetMarkerStyle(20);
new_elebkgMt->SetFillStyle(1001);
new_elebkgMt->SetLineColor(kRed);
new_elebkgMt->SetFillColor(kRed);
new_elebkgMt->Draw("hist same");
new_elebkgMtTT->Add(new_elebkgMt);
new_elebkgMtTT->SetLineColor(kGreen);
new_elebkgMtTT->SetFillStyle(1001);
new_elebkgMtTT->SetFillColor(kGreen);
new_elebkgMtTT->Draw("hist same");
new_elebkgMt->Draw("hist same");
TLegend *leg2 =  new TLegend(0.6,0.7,0.9,0.9);
leg2->SetFillStyle(0);
gStyle->SetLegendBorderSize(1);
gStyle->SetLegendFillColor(0);
leg2->AddEntry(new_allMt,"observed");
leg2->AddEntry(new_elebkgMtTT,"EWK (tt, WW)");
leg2->AddEntry(new_elebkgMt,"DY");
leg2->Draw("same");
new_allMt->Draw("same");
c_mt->SaveAs("EfakeMT_closure.pdf");
}


