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
#include "TVector2.h"

#include "../include/analysis_rawData.h"
#include "../include/analysis_photon.h"
#include "../include/analysis_muon.h"
#include "../include/analysis_ele.h"
#include "../include/analysis_mcData.h"
#include "../include/analysis_tools.h"
#include "../include/analysis_jet.h"

void test_leak(){//main  

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

  TChain* es = new TChain("SUSYtree");
  es->Add("/uscms_data/d3/mengleis/Sep1/resTree_TChiWG.root");

	float Mgluino(0);
  float Mchargino(0);
  float Mneutralino(0);
  float mcPhotonEt(0);
  float mcPhotonEta(0);
  float mcPhotonPhi(0);
  float recoPhotonEt(0);
  float recoPhotonEta(0);
  float recoPhotonPhi(0);
  float PhoR9(0);
  float PhoHoverE(0);
  float PhoSigma(0);
  float PhoChIso(0);
  float PhoNeuIso(0);
  float PhoPhoIso(0);
  bool  PhoPassID;
  float mcElePt(0);
  float mcEleEta(0);
  float mcElePhi(0);
  float recoEleEt(0);
  float recoEleEta(0);
  float recoElePhi(0);
  float EleR9(0);
  float EleHoverE(0);
  float EleSigma(0);
  float EleChIso(0);
  float EleNeuIso(0);
  float ElePhoIso(0);
  float EleMiniIso(0);
  float EledEtaIn(0);
  float EledPhiIn(0);
  float EleD0(0);
  float EleDz(0);
  float EleooEmooP(0);
  bool  ElePassID;
  float dRPhoEle(0);
  float Invmass(0);
  int   nJet(0);
  float HT(0);
	float sigMET(0);
  double MT_(0), ThreeBodyMass_(0);
	int   nVertex(0);

  es->SetBranchAddress("Mgluino",        &Mgluino);
  es->SetBranchAddress("Mchargino",      &Mchargino);
  es->SetBranchAddress("Mneutralino",    &Mneutralino);
	es->SetBranchAddress("nVertex",        &nVertex);
  es->SetBranchAddress("MT",&MT_);
  es->SetBranchAddress("ThreeBodyMass",&ThreeBodyMass_);
  es->SetBranchAddress("mcPhotonEt"     ,&mcPhotonEt); 
  es->SetBranchAddress("mcPhotonEta"    ,&mcPhotonEta);
  es->SetBranchAddress("mcPhotonPhi"    ,&mcPhotonPhi);
  es->SetBranchAddress("recoPhotonEt",   &recoPhotonEt);
  es->SetBranchAddress("recoPhotonEta",  &recoPhotonEta);
  es->SetBranchAddress("recoPhotonPhi",  &recoPhotonPhi);
  es->SetBranchAddress("PhoR9"          ,&PhoR9);
  es->SetBranchAddress("PhoHoverE"      ,&PhoHoverE);
  es->SetBranchAddress("PhoSigma"       ,&PhoSigma);
  es->SetBranchAddress("PhoChIso"       ,&PhoChIso);
  es->SetBranchAddress("PhoNeuIso"      ,&PhoNeuIso);
  es->SetBranchAddress("PhoPhoIso"      ,&PhoPhoIso);
  es->SetBranchAddress("PhoPassID"      ,&PhoPassID);
  es->SetBranchAddress("mcElePt"        ,&mcElePt);
  es->SetBranchAddress("mcEleEta"       ,&mcEleEta);
  es->SetBranchAddress("mcElePhi"       ,&mcElePhi);
  es->SetBranchAddress("recoEleEt",     &recoEleEt);
  es->SetBranchAddress("recoEleEta",    &recoEleEta);
  es->SetBranchAddress("recoElePhi",    &recoElePhi);
  es->SetBranchAddress("EleR9"          ,&EleR9);
  es->SetBranchAddress("EleHoverE"      ,&EleHoverE);
  es->SetBranchAddress("EleSigma"       ,&EleSigma);
  es->SetBranchAddress("EleChIso"       ,&EleChIso);
  es->SetBranchAddress("EleNeuIso"      ,&EleNeuIso);
  es->SetBranchAddress("ElePhoIso"      ,&ElePhoIso);
  es->SetBranchAddress("EleMiniIso"     ,&EleMiniIso);
  es->SetBranchAddress("EledEtaIn"      ,&EledEtaIn);
  es->SetBranchAddress("EledPhiIn"      ,&EledPhiIn);
  es->SetBranchAddress("EleD0"          ,&EleD0);
  es->SetBranchAddress("EleDz"          ,&EleDz);
  es->SetBranchAddress("EleooEmooP"     ,&EleooEmooP);
  es->SetBranchAddress("ElePassID"      ,&ElePassID);
  es->SetBranchAddress("dRPhoEle"       ,&dRPhoEle);
  es->SetBranchAddress("Invmass"        ,&Invmass);    
  es->SetBranchAddress("nJet"           ,&nJet);
  es->SetBranchAddress("HT"             ,&HT);
	es->SetBranchAddress("sigMET",         &sigMET);

	TH1D *h1=new TH1D("h1","",1000,0,2000);
	TH1D *h2=new TH1D("h2","",1000,0,2000);
	for(unsigned ievt(0); ievt < es->GetEntries(); ievt++){
		es->GetEntry(ievt);
		if(Mgluino > 1000)continue;
		if(recoPhotonEt > 35 && fabs(recoPhotonEta) < 1.4442 )h1->Fill(recoPhotonEt);
		if(recoPhotonEt > 35 && fabs(recoPhotonEta) < 1.4442 && PhoPassID > 0)h2->Fill(recoPhotonEt);
	}
	h2->Divide(h1);
	h2->Draw(); 


}
