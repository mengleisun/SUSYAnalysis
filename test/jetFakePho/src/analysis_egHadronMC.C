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

#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_tools.h"
#include "../../../include/analysis_mcData.h"

void analysis_egHadronMC(){//main  

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

  char outputname[100] = "/uscms_data/d3/mengleis/Sep1/plot_hadron_GJet.root";
  ofstream logfile;
  logfile.open("/uscms_data/d3/mengleis/Sep1/plot_hadron_GJet.log"); 

  logfile << "analysis_hadron()" << std::endl;

  RunType datatype(MC); 

  TChain* es = new TChain("ggNtuplizer/EventTree");
  es->Add("root://cmseos.fnal.gov///store/group/lpcsusystealth/ggNtuple_leppho/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1.root");
  logfile << "GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1.root" << std::endl;

  TFile *outputfile = TFile::Open(outputname,"NEW");
  outputfile->cd();

  int mcType = MCType::WGJetInclusive;
  if(datatype == MC && mcType == MCType::NOMC){std::cout << "wrong MC type" << std::endl; throw;} 
  logfile << "mcType" << mcType << std::endl;
//************ Signal Tree **********************//
  TTree *egtree = new TTree("egTree","egTree");
  float eg_phoEt(0);
  float eg_phoEta(0);
  float eg_phoPhi(0);
  float eg_phoSigma(0);
  float eg_phoChIso(0);
  float eg_sigMT(0);
  float eg_sigMET(0);
  float eg_sigMETPhi(0);
  float eg_dPhiLepMET(0);
  int   eg_nVertex(0);
  float eg_HT(0);
  float eg_nJet(0);
  std::vector<int> eg_mcPID;
  std::vector<float> eg_mcEta;
  std::vector<float> eg_mcPhi;
  std::vector<float> eg_mcPt;
  std::vector<int> eg_mcMomPID;
  
  egtree->Branch("phoEt",     &eg_phoEt);
  egtree->Branch("phoEta",    &eg_phoEta);
  egtree->Branch("phoPhi",    &eg_phoPhi);
  egtree->Branch("phoSigma",  &eg_phoSigma);
  egtree->Branch("phoChIso",  &eg_phoChIso);
  egtree->Branch("sigMT",     &eg_sigMT);
  egtree->Branch("sigMET",    &eg_sigMET);
  egtree->Branch("sigMETPhi", &eg_sigMETPhi);
  egtree->Branch("dPhiLepMET",&eg_dPhiLepMET);
  egtree->Branch("nVertex",   &eg_nVertex);
  egtree->Branch("HT",        &eg_HT);
  egtree->Branch("nJet",      &eg_nJet);
  egtree->Branch("mcPID",     &eg_mcPID);
  egtree->Branch("mcEta",     &eg_mcEta);
  egtree->Branch("mcPhi",     &eg_mcPhi);
  egtree->Branch("mcPt",      &eg_mcPt);
  egtree->Branch("mcMomPID",  &eg_mcMomPID);
  egtree->Branch("mcType",    &mcType);

  TTree *hadtree = new TTree("hadTree","hadTree");
  float had_phoEt(0);
  float had_phoEta(0);
  float had_phoPhi(0);
  float had_phoSigma(0);
  float had_phoChIso(0);
  float had_sigMT(0);
  float had_sigMET(0);
  float had_sigMETPhi(0);
  float had_dPhiLepMET(0);
  int   had_nVertex(0);
  float had_HT(0);
  float had_nJet(0);
  
  hadtree->Branch("phoEt",     &had_phoEt);
  hadtree->Branch("phoEta",    &had_phoEta);
  hadtree->Branch("phoPhi",    &had_phoPhi);
  hadtree->Branch("phoSigma",  &had_phoSigma);
  hadtree->Branch("phoChIso",  &had_phoChIso);
  hadtree->Branch("sigMT",     &had_sigMT);
  hadtree->Branch("sigMET",    &had_sigMET);
  hadtree->Branch("sigMETPhi", &had_sigMETPhi);
  hadtree->Branch("dPhiLepMET",&had_dPhiLepMET);
  hadtree->Branch("nVertex",   &had_nVertex);
  hadtree->Branch("HT",        &had_HT);
  hadtree->Branch("nJet",      &had_nJet);

  rawData raw(es, datatype);
  std::vector<mcData>  MCData;
  std::vector<recoPhoton> Photon;
  std::vector<recoMuon>   Muon;
  std::vector<recoEle>   Ele;
  float MET(0);
  float METPhi(0);
  int nVtx(0);
  int jetNumber(0);
  int METFilter(0);
  logfile << "RunType: " << datatype << std::endl;

  const unsigned nEvts = es->GetEntries(); 
  std::cout << "total " << nEvts << std::endl;
  logfile <<   "total " << nEvts << std::endl;

  for(unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
  
    if(ievt%10000==0)std::cout << " -- Processing event " << ievt << std::endl;

	raw.GetData(es, ievt);
	MCData.clear();
	Photon.clear();
	Muon.clear();
	Ele.clear();
	if(datatype == MC)for(int iMC(0); iMC < raw.nMC; iMC++){MCData.push_back(mcData(raw, iMC));}
	for(int iPho(0); iPho < raw.nPho; iPho++){Photon.push_back(recoPhoton(raw, iPho));}
	for(int iMu(0); iMu < raw.nMu; iMu++){Muon.push_back(recoMuon(raw, iMu));}
	for(int iEle(0); iEle < raw.nEle; iEle++){Ele.push_back(recoEle(raw, iEle));}
	MET = raw.pfMET;
	METPhi = raw.pfMETPhi;
	METFilter = raw.metFilters;
	nVtx = raw.nVtx;
	jetNumber = raw.nJet;

	if(!raw.passHLT())continue;
	if(raw.nPho <1)continue;
    if(MET>70)continue;
   
	//bool hasEle(false);
	//std::vector<recoEle>::iterator signalEle = Ele.begin();
	//for(std::vector<recoEle>::iterator itEle = Ele.begin(); itEle != Ele.end(); itEle++){
	//  if(itEle->passSignalSelection()){
	//	if(!hasEle){
	//	  hasEle=true; 
	//	  signalEle = itEle;
	//	}
	//  }
	//}
      
    //if(!hasEle)continue;
    if(Photon.size() < 2)continue;
    if(Photon[0].getCalibEt() < 25 || Photon[1].getCalibEt() < 25)continue;
    //Get reco-level photons
	for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
      bool isegCandidate(false), ishadCandidate(false);
      bool hasegCandidate(false),  hashadCandidate(false);
	  if(itpho->getCalibEt() < 25 || !itpho->passHoverE(1)|| !itpho->passNeuIso(1) || !itpho->passPhoIso(1))continue;
      if(itpho->getChIso() > 20.0)continue;
      //double DeltaPhoEle = DeltaR(itpho->getEta(), itpho->getPhi(), signalEle->getEta(), signalEle->getPhi());
      //if(DeltaPhoEle < 0.8)continue;
	  bool PixelVeto = itpho->PixelSeed()==0? true: false;
	  bool GSFveto(true);
	  bool FSRVeto(true);
	  for(std::vector<recoEle>::iterator ie = Ele.begin(); ie != Ele.end(); ie++){
		 if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.02)GSFveto = false;
		 if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.3 && ie->getCalibPt()>2.0)FSRVeto=false;
	  }
	  for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++)
		 if(DeltaR(itpho->getEta(), itpho->getPhi(), im->getEta(), im->getPhi()) < 0.3 && im->getPt()>2.0)FSRVeto=false;
	  if(GSFveto && PixelVeto && FSRVeto){
        isegCandidate = true;
        ishadCandidate = true;
	  }
      if(!itpho->fireDoubleTrg(5) && !itpho->fireDoubleTrg(6))isegCandidate=false;

	  //float deltaPhi = DeltaPhi(signalEle->getPhi(), METPhi);
	  //float MT = sqrt(2*MET*signalEle->getPt()*(1-std::cos(deltaPhi)));
      if(isegCandidate && !hasegCandidate){
		eg_phoEt = itpho->getCalibEt();
		eg_phoEta = itpho->getEta();
		eg_phoPhi = itpho->getPhi();
		eg_phoSigma = itpho->getSigma();
		eg_phoChIso = itpho->getChIso();
		//eg_sigMT = MT;
		eg_sigMET = MET;
		eg_sigMETPhi = METPhi;
		//eg_dPhiLepMET = deltaPhi;
		eg_nVertex = nVtx;
		eg_HT = 0;
		eg_nJet = jetNumber;
		eg_mcPID.clear();
		eg_mcEta.clear();
		eg_mcPhi.clear();
		eg_mcPt.clear();
		eg_mcMomPID.clear();
		if(datatype == MC){
		  for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
			if(itMC->getEt() < 1.0)continue;
			float mcdR = DeltaR(itpho->getEta(), itpho->getPhi(), itMC->getEta(), itMC->getPhi());
			if(mcdR < 0.3){
			  eg_mcPID.push_back(itMC->getPID());
			  eg_mcMomPID.push_back(itMC->getMomPID());
			  eg_mcEta.push_back(itMC->getEta());      
			  eg_mcPhi.push_back(itMC->getPhi());
			  eg_mcPt.push_back(itMC->getEt());
			}
		  }
		}
        egtree->Fill();
        hasegCandidate = true;
      }
      if(ishadCandidate && !hashadCandidate){
		had_phoEt = itpho->getCalibEt();
		had_phoEta = itpho->getEta();
		had_phoPhi = itpho->getPhi();
		had_phoSigma = itpho->getSigma();
		had_phoChIso = itpho->getChIso();
	//	had_sigMT = MT;
		had_sigMET = MET;
		had_sigMETPhi = METPhi;
//		had_dPhiLepMET = deltaPhi;
		had_nVertex = nVtx;
		had_HT = 0;
		had_nJet = jetNumber;
        hadtree->Fill();
        hashadCandidate = true;
      }
	}

  } 
outputfile->Write();
}


