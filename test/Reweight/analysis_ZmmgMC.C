#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<ctime>

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
#include "TFileCollection.h"

#include "../../include/analysis_rawData.h"
#include "../../include/analysis_photon.h"
#include "../../include/analysis_muon.h"
#include "../../include/analysis_ele.h"
#include "../../include/analysis_mcData.h"
#include "../../include/analysis_tools.h"
#include "../../include/analysis_jet.h"

bool passMETFilter(int filter){
  bool passfilter(true);
  for(int im(1); im <= 8; im++)
    if(((filter >> im)&1)!=0){passfilter = false; return passfilter;}

//  if(((filter >> 9)&1)!=1){passfilter = false; return passfilter;}
//  if(((filter >> 10)&1)!=1){passfilter = false; return passfilter;}

  return passfilter;
}


void analysis_ZmmgMC(){//main  

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

  char outputname[100] = "/uscms_data/d3/mengleis/test/resTree_ZGamma_ZG130.root";
  ofstream logfile;
  logfile.open("/uscms_data/d3/mengleis/test/resTree_ZGamma_ZG130.log"); 

  logfile << "analysis_mg()" << std::endl;

  RunType datatype(MC); 
  TChain* es = new TChain("ggNtuplizer/EventTree");
//	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_MuonEG_FebReminiAOD_BCD.root");
//	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_MuonEG_FebReminiAOD_EFG.root");
//	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_MuonEG_FebReminiAOD_H.root");
//
	//es->Add("root://cmseos.fnal.gov//store/user/msun/MCSummer16/ZGTo2LG_RunIISummer16MiniAODv2-TrancheIV_v6-v1.root");

//	es->Add("root://cmseos.fnal.gov//store/user/msun/MCSummer16/DYJetsToLL_M-50_nlo.root");
//	es->Add("root://cmseos.fnal.gov//store/user/msun/MCSummer16/TTGJets_RunIISummer16MiniAODv2-TrancheIV_v6_ext1.root");
//	es->Add("root://cmseos.fnal.gov//store/user/msun/MCSummer16/WWG_RunIISummer16MiniAODv2-TrancheIV_v6_ext1.root");
	//
	//need to rerun
	//es->Add("root://cmseos.fnal.gov//store/user/msun/MCSummer16/WZG_RunIISummer16MiniAODv2-TrancheIV_v6.root");
	es->Add("root://cmseos.fnal.gov//store/user/msun/MCSummer16/ZGTo2LG_PtG-130_RunIISummer16MiniAODv2-TrancheIV_v6-v1.root");

  const unsigned nEvts = es->GetEntries(); 
	//float MCweight = 4895*35.8*1000.0/nEvts;	
	//float MCweight = 3.697*35.8*1000.0/nEvts;	
	//float MCweight = 0.21*35.8*1000.0/nEvts;	
	//float MCweight = 0.04*35.8*1000.0/nEvts;//WZG
	float MCweight = 177.8*35.8*1000.0/nEvts;	
	//float MCweight = 0.14*35.8*1000.0/nEvts;//WZG
  logfile << "Total event: " << nEvts << std::endl;
  std::cout << "Total event: " << nEvts << std::endl;
  logfile << "Output file: " << outputname << std::endl;

	int nTotal(0),npassHLT(0), npassPho(0), npassLep(0), npassdR(0), npassZ(0), npassMETFilter(0);

  TFile *outputfile = TFile::Open(outputname,"RECREATE");
  outputfile->cd();

	TH1F *p_nMu = new TH1F("p_nMu","p_nMu",10,0,10);
//************ Signal Tree **********************//
  TTree *sigtree = new TTree("signalTree","signalTree");
  int   run(0);
  Long64_t  event(0);
  int   lumis(0);
  float phoEt(0);
  float phoEta(0);
  float phoPhi(0);
  float lepPt(0);
  float lepEta(0);
  float lepPhi(0);
  float sigMT(0);
  float sigMET(0);
  float sigMETPhi(0);
  float dPhiLepMET(0);
	float threeMass(0);
  int   nVertex(0);
  float dRPhoLep(0);
  float HT(0);
  float nJet(0);

  sigtree->Branch("run",       &run);
  sigtree->Branch("event",     &event);
  sigtree->Branch("lumis",     &lumis);
  sigtree->Branch("phoEt",     &phoEt);
  sigtree->Branch("phoEta",    &phoEta);
  sigtree->Branch("phoPhi",    &phoPhi);
  sigtree->Branch("lepPt",     &lepPt);
  sigtree->Branch("lepEta",    &lepEta);
  sigtree->Branch("lepPhi",    &lepPhi);
  sigtree->Branch("sigMT",     &sigMT);
  sigtree->Branch("sigMET",    &sigMET);
  sigtree->Branch("sigMETPhi", &sigMETPhi);
  sigtree->Branch("dPhiLepMET",&dPhiLepMET);
	sigtree->Branch("threeMass", &threeMass);
  sigtree->Branch("nVertex",   &nVertex);
  sigtree->Branch("dRPhoLep",  &dRPhoLep);
  sigtree->Branch("HT",        &HT);
  sigtree->Branch("nJet",      &nJet);

//************ Signal Tree **********************//
  TTree *proxytree = new TTree("proxyTree","proxyTree");
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
	float proxythreeMass(0);
  int   proxynVertex(0);
  float proxydRPhoLep(0);
  float proxyHT(0);
  float proxynJet(0);
	float proxysigMETJESup(0);
	float proxysigMETJESdo(0);
	float proxysigMETJERup(0);
	float proxysigMETJERdo(0);
	float proxysigMTJESup(0);
	float proxysigMTJESdo(0);
	float proxysigMTJERup(0);
	float proxysigMTJERdo(0);
	float proxyHTJESup(0);
	float proxyHTJESdo(0);
  
  proxytree->Branch("phoEt",     &proxyphoEt);
  proxytree->Branch("phoEta",    &proxyphoEta);
  proxytree->Branch("phoPhi",    &proxyphoPhi);
  proxytree->Branch("lepPt",     &proxylepPt);
  proxytree->Branch("lepEta",    &proxylepEta);
  proxytree->Branch("lepPhi",    &proxylepPhi);
  proxytree->Branch("sigMT",     &proxysigMT);
  proxytree->Branch("sigMET",    &proxysigMET);
  proxytree->Branch("sigMETPhi", &proxysigMETPhi);
  proxytree->Branch("dPhiLepMET",&proxydPhiLepMET);
	proxytree->Branch("threeMass", &proxythreeMass);
  proxytree->Branch("nVertex",   &proxynVertex);
  proxytree->Branch("dRPhoLep",  &proxydRPhoLep);
  proxytree->Branch("HT",        &proxyHT);
  proxytree->Branch("nJet",      &proxynJet);
	proxytree->Branch("sigMETJESup",     &proxysigMETJESup);
	proxytree->Branch("sigMETJESdo",     &proxysigMETJESdo);
	proxytree->Branch("sigMETJERup",     &proxysigMETJERup);
	proxytree->Branch("sigMETJERdo",     &proxysigMETJERdo);
	proxytree->Branch("sigMTJESup",      &proxysigMTJESup);
	proxytree->Branch("sigMTJESdo",      &proxysigMTJESdo);
	proxytree->Branch("sigMTJERup",      &proxysigMTJERup);
	proxytree->Branch("sigMTJERdo",      &proxysigMTJERdo);
	proxytree->Branch("HTJESup",     &proxyHTJESup);
	proxytree->Branch("HTJESdo",     &proxyHTJESdo);

//************ Signal Tree **********************//
  TTree *jettree = new TTree("jetTree","jetTree");
  float jetphoEt(0);
  float jetphoEta(0);
  float jetphoPhi(0);
  float jetlepPt(0);
  float jetlepEta(0);
  float jetlepPhi(0);
  float jetsigMT(0);
  float jetsigMET(0);
  float jetsigMETPhi(0);
  float jetdPhiLepMET(0);
	float jetthreeMass(0);
  int   jetnVertex(0);
  float jetdRPhoLep(0);
  float jetHT(0);
  float jetnJet(0);
	float jetsigMETJESup(0);
	float jetsigMETJESdo(0);
	float jetsigMETJERup(0);
	float jetsigMETJERdo(0);
	float jetsigMTJESup(0);
	float jetsigMTJESdo(0);
	float jetsigMTJERup(0);
	float jetsigMTJERdo(0);
	float jetHTJESup(0);
	float jetHTJESdo(0);
  
  jettree->Branch("phoEt",     &jetphoEt);
  jettree->Branch("phoEta",    &jetphoEta);
  jettree->Branch("phoPhi",    &jetphoPhi);
  jettree->Branch("lepPt",     &jetlepPt);
  jettree->Branch("lepEta",    &jetlepEta);
  jettree->Branch("lepPhi",    &jetlepPhi);
  jettree->Branch("sigMT",     &jetsigMT);
  jettree->Branch("sigMET",    &jetsigMET);
  jettree->Branch("sigMETPhi", &jetsigMETPhi);
  jettree->Branch("dPhiLepMET",&jetdPhiLepMET);
	jettree->Branch("threeMass", &jetthreeMass);
  jettree->Branch("nVertex",   &jetnVertex);
  jettree->Branch("dRPhoLep",  &jetdRPhoLep);
  jettree->Branch("HT",        &jetHT);
  jettree->Branch("nJet",      &jetnJet);
	jettree->Branch("sigMETJESup",     &jetsigMETJESup);
	jettree->Branch("sigMETJESdo",     &jetsigMETJESdo);
	jettree->Branch("sigMETJERup",     &jetsigMETJERup);
	jettree->Branch("sigMETJERdo",     &jetsigMETJERdo);
	jettree->Branch("sigMTJESup",      &jetsigMTJESup);
	jettree->Branch("sigMTJESdo",      &jetsigMTJESdo);
	jettree->Branch("sigMTJERup",      &jetsigMTJERup);
	jettree->Branch("sigMTJERdo",      &jetsigMTJERdo);
	jettree->Branch("HTJESup",     &jetHTJESup);
	jettree->Branch("HTJESdo",     &jetHTJESdo);
//*********** fake lepton *********************//
  TTree *fakeLeptree = new TTree("fakeLepTree","fakeLepTree");
  float fakeLepphoEt(0);
  float fakeLepphoEta(0);
  float fakeLepphoPhi(0);
  float fakeLepPt(0);
  float fakeLepEta(0);
  float fakeLepPhi(0);
	float fakeLepMiniIso(0);
	int   fakeLepIsMedium(0);
  float fakeLepsigMT(0);
  float fakeLepsigMET(0);
  float fakeLepsigMETPhi(0);
  float fakeLepdPhiLepMET(0);
	float fakethreeMass(0);
  int   fakeLepnVertex(0);
  float fakeLepdRPhoLep(0);
  float fakeLepHT(0);
  float fakeLepnJet(0);
	float fakeLepsigMETJESup(0);
	float fakeLepsigMETJESdo(0);
	float fakeLepsigMETJERup(0);
	float fakeLepsigMETJERdo(0);
	float fakeLepsigMTJESup(0);
	float fakeLepsigMTJESdo(0);
	float fakeLepsigMTJERup(0);
	float fakeLepsigMTJERdo(0);
	float fakeLepHTJESup(0);
	float fakeLepHTJESdo(0);
  
  fakeLeptree->Branch("phoEt",     &fakeLepphoEt);
  fakeLeptree->Branch("phoEta",    &fakeLepphoEta);
  fakeLeptree->Branch("phoPhi",    &fakeLepphoPhi);
  fakeLeptree->Branch("lepPt",     &fakeLepPt);
  fakeLeptree->Branch("lepEta",    &fakeLepEta);
  fakeLeptree->Branch("lepPhi",    &fakeLepPhi);
  fakeLeptree->Branch("fakeLepMiniIso",&fakeLepMiniIso);
  fakeLeptree->Branch("fakeLepIsMedium",&fakeLepIsMedium);
  fakeLeptree->Branch("sigMT",     &fakeLepsigMT);
  fakeLeptree->Branch("sigMET",    &fakeLepsigMET);
  fakeLeptree->Branch("sigMETPhi", &fakeLepsigMETPhi);
  fakeLeptree->Branch("dPhiLepMET",&fakeLepdPhiLepMET);
	fakeLeptree->Branch("threeMass", &fakethreeMass);
  fakeLeptree->Branch("nVertex",   &fakeLepnVertex);
  fakeLeptree->Branch("dRPhoLep",  &fakeLepdRPhoLep);
  fakeLeptree->Branch("HT",        &fakeLepHT);
  fakeLeptree->Branch("nJet",      &fakeLepnJet);
	fakeLeptree->Branch("sigMETJESup",     &fakeLepsigMETJESup);
	fakeLeptree->Branch("sigMETJESdo",     &fakeLepsigMETJESdo);
	fakeLeptree->Branch("sigMETJERup",     &fakeLepsigMETJERup);
	fakeLeptree->Branch("sigMETJERdo",     &fakeLepsigMETJERdo);
	fakeLeptree->Branch("sigMTJESup",      &fakeLepsigMTJESup);
	fakeLeptree->Branch("sigMTJESdo",      &fakeLepsigMTJESdo);
	fakeLeptree->Branch("sigMTJERup",      &fakeLepsigMTJERup);
	fakeLeptree->Branch("sigMTJERdo",      &fakeLepsigMTJERdo);
	fakeLeptree->Branch("HTJESup",     &fakeLepHTJESup);
	fakeLeptree->Branch("HTJESdo",     &fakeLepHTJESdo);

//************ ZGamma Tree **********************//
  TTree *Ztree = new TTree("ZTree","ZTree");
  float ZphoEt(0);
  float ZphoEta(0);
  float ZphoPhi(0);
	float ZphoChIso(0);
	float ZphoNeuIso(0);
  float ZlepPt(0);
  float ZlepEta(0);
  float ZlepPhi(0);
	float ZtrailPt(0);
	float ZtrailEta(0);
	float ZtrailPhi(0);
  float ZsigMT(0);
  float ZsigMET(0);
  float ZsigMETPhi(0);
  float ZdPhiLepMET(0);
	float ZthreeMass(0);
	float ZcalibthreeMass(0);
	float ZdilepMass(0);
  int   ZnVertex(0);
  float ZdRPhoLep(0);
	float ZdRPhoTrail(0);
  float ZHT(0);
  float ZnJet(0);
	int   ZnMu(0);
	int   ZnGoodMu(0);
  std::vector<int> Z_mcPID;
  std::vector<float> Z_mcEta;
  std::vector<float> Z_mcPhi;
  std::vector<float> Z_mcPt;
  std::vector<int> Z_mcMomPID;
  std::vector<int> Z_mcGMomPID;
	std::vector<int> Z_mcStatus;
 
	Ztree->Branch("MCweight",  &MCweight); 
  Ztree->Branch("phoEt",     &ZphoEt);
  Ztree->Branch("phoEta",    &ZphoEta);
  Ztree->Branch("phoPhi",    &ZphoPhi);
	Ztree->Branch("phoChIso",  &ZphoChIso);
	Ztree->Branch("phoNeuIso", &ZphoNeuIso);
  Ztree->Branch("lepPt",     &ZlepPt);
  Ztree->Branch("lepEta",    &ZlepEta);
  Ztree->Branch("lepPhi",    &ZlepPhi);
	Ztree->Branch("trailPt",   &ZtrailPt);
	Ztree->Branch("trailEta",  &ZtrailEta);
	Ztree->Branch("trailPhi",  &ZtrailPhi);
  Ztree->Branch("sigMT",     &ZsigMT);
  Ztree->Branch("sigMET",    &ZsigMET);
  Ztree->Branch("sigMETPhi", &ZsigMETPhi);
  Ztree->Branch("dPhiLepMET",&ZdPhiLepMET);
	Ztree->Branch("threeMass", &ZthreeMass);
	Ztree->Branch("calibthreeMass", &ZcalibthreeMass);
	Ztree->Branch("dilepMass", &ZdilepMass);
  Ztree->Branch("nVertex",   &ZnVertex);
  Ztree->Branch("dRPhoLep",  &ZdRPhoLep);
	Ztree->Branch("dRPhoTrail",&ZdRPhoTrail);
  Ztree->Branch("HT",        &ZHT);
  Ztree->Branch("nJet",      &ZnJet);
	Ztree->Branch("nMu",       &ZnMu);
	Ztree->Branch("nGoodMu",   &ZnGoodMu);
  Ztree->Branch("mcPID",    &Z_mcPID);
  Ztree->Branch("mcEta",    &Z_mcEta);
  Ztree->Branch("mcPhi",    &Z_mcPhi);
  Ztree->Branch("mcPt",     &Z_mcPt);
  Ztree->Branch("mcMomPID", &Z_mcMomPID);
  Ztree->Branch("mcGMomPID",&Z_mcGMomPID);
	Ztree->Branch("mcStatus", &Z_mcStatus);
//************ ZJet Tree **********************//
  TTree *ZJettree = new TTree("ZJetTree","ZJetTree");
  float ZJetphoEt(0);
  float ZJetphoEta(0);
  float ZJetphoPhi(0);
  float ZJetlepPt(0);
  float ZJetlepEta(0);
  float ZJetlepPhi(0);
  float ZJetsigMT(0);
  float ZJetsigMET(0);
  float ZJetsigMETPhi(0);
  float ZJetdPhiLepMET(0);
	float ZJetthreeMass(0);
	float ZJetdilepMass(0);
	float ZJetphoSigma(0);
	float ZJetphoChIso(0);
  int   ZJetnVertex(0);
  float ZJetdRPhoLep(0);
  float ZJetHT(0);
  float ZJetnJet(0);
  
  ZJettree->Branch("phoEt",     &ZJetphoEt);
  ZJettree->Branch("phoEta",    &ZJetphoEta);
  ZJettree->Branch("phoPhi",    &ZJetphoPhi);
  ZJettree->Branch("lepPt",     &ZJetlepPt);
  ZJettree->Branch("lepEta",    &ZJetlepEta);
  ZJettree->Branch("lepPhi",    &ZJetlepPhi);
  ZJettree->Branch("sigMT",     &ZJetsigMT);
  ZJettree->Branch("sigMET",    &ZJetsigMET);
  ZJettree->Branch("sigMETPhi", &ZJetsigMETPhi);
  ZJettree->Branch("dPhiLepMET",&ZJetdPhiLepMET);
	ZJettree->Branch("threeMass", &ZJetthreeMass);
	ZJettree->Branch("dilepMass", &ZJetdilepMass);
	ZJettree->Branch("phoSigma",  &ZJetphoSigma);
	ZJettree->Branch("phoChIso",  &ZJetphoChIso);
  ZJettree->Branch("nVertex",   &ZJetnVertex);
  ZJettree->Branch("dRPhoLep",  &ZJetdRPhoLep);
  ZJettree->Branch("HT",        &ZJetHT);
  ZJettree->Branch("nJet",      &ZJetnJet);

//*********** histo list **********************//
  TH1F *p_eventcount = new TH1F("p_eventcount","p_eventcount",7,0,7);

  rawData raw(es, datatype);
  std::vector<mcData>  MCData;
  std::vector<recoPhoton> Photon;
  std::vector<recoMuon>   Muon;
  std::vector<recoEle>   Ele;
  std::vector<recoJet>   JetCollection;
  float MET(0);
  float METPhi(0);
	float MET_T1JERUp(0);
	float MET_T1JERDo(0);
	float MET_T1JESUp(0);
	float MET_T1JESDo(0);	
	float	METPhi_T1JESUp(0);
	float	METPhi_T1JESDo(0);
	float	METPhi_T1UESUp(0);
	float	METPhi_T1UESDo(0);
  int nVtx(0);
  int jetNumber(0);
  int METFilter(0);
  logfile << "RunType: " << datatype << std::endl;

  std::cout << "Total evetns : " << nEvts << std::endl;
  logfile << "Total evetns : " << nEvts << std::endl;
	for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries

		if (ievt%100000==0) std::cout << " -- Processing event " << ievt << std::endl;
		if (ievt%100000==0) logfile  << " -- Processing event " << ievt << std::endl;

			raw.GetData(es, ievt);
			MCData.clear();
			Photon.clear();
			Muon.clear();
			Ele.clear();
			JetCollection.clear();
			if(datatype == MC)for(int iMC(0); iMC < raw.nMC; iMC++){MCData.push_back(mcData(raw, iMC));}
			for(int iPho(0); iPho < raw.nPho; iPho++){Photon.push_back(recoPhoton(raw, iPho));}
			for(int iMu(0); iMu < raw.nMu; iMu++){Muon.push_back(recoMuon(raw, iMu));}
			for(int iEle(0); iEle < raw.nEle; iEle++){Ele.push_back(recoEle(raw, iEle));}
			for(int iJet(0); iJet < raw.nJet; iJet++){JetCollection.push_back(recoJet(raw, iJet));}
			MET = raw.pfMET;
			METPhi = raw.pfMETPhi;
			METFilter = raw.metFilters;
			MET_T1JERUp = raw.pfMET_T1JERUp;
			MET_T1JERDo = raw.pfMET_T1JERDo;
			MET_T1JESUp = raw.pfMET_T1JESUp;
			MET_T1JESDo = raw.pfMET_T1JESDo;
			METPhi_T1JESUp = raw.pfMETPhi_T1JESUp;
			METPhi_T1JESDo = raw.pfMETPhi_T1JESDo;
			METPhi_T1UESUp = raw.pfMETPhi_T1UESUp;
			METPhi_T1UESDo = raw.pfMETPhi_T1UESDo;
			nVtx = raw.nVtx;
			run=raw.run;
			event=raw.event;
			lumis=raw.lumis;

			nTotal+=1;
      if(((raw.HLTEleMuX >> 8) &1) ==0 && ((raw.HLTEleMuX >> 41) &1) ==0)continue;   // Mu+gamma HLT
      //if(((raw.HLTEleMuX >> 19) &1)==0 && ((raw.HLTEleMuX >> 20) &1)==0)continue;  // SingleMuon HLT
			if(raw.nGoodVtx < 1)continue;
			npassHLT+=1;

			int nGoodMu(0);
			for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++){
			  if(im->isMedium()){nGoodMu+=1;}
			}	
			p_nMu->Fill(nGoodMu);

			if(raw.nMu < 2 || raw.nPho <1)continue;


			bool hasPho(false);
			std::vector<recoPhoton>::iterator signalPho = Photon.begin();
			std::vector< std::vector<recoPhoton>::iterator >  proxyPhoCollection;
			proxyPhoCollection.clear();
			std::vector< std::vector<recoPhoton>::iterator >  jetPhoCollection;
			jetPhoCollection.clear();
			for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
				if(itpho->getR9() < 0.5)continue;
				if(itpho->fireDoubleTrg(28) || itpho->fireDoubleTrg(29) || itpho->fireDoubleTrg(30)){
					if(!itpho->passBasicSelection())continue;
					bool passSigma = itpho->passSigma(1);
					bool passChIso = itpho->passChIso(1);
					bool PixelVeto = itpho->PixelSeed()==0? true: false;
					bool GSFveto(true);
					bool FSRVeto(true);
					for(std::vector<recoEle>::iterator ie = Ele.begin(); ie != Ele.end(); ie++){
						if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) <= 0.03)GSFveto = false;
						if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.3)FSRVeto=false;
					}
					for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++)
						if(DeltaR(itpho->getEta(), itpho->getPhi(), im->getEta(), im->getPhi()) < 0.3 && im->getEt()>2.0)FSRVeto=false;

					if(itpho->getChIso()<20 && itpho->getSigma()< 0.02 && itpho->isEB()){
						//if(itpho->isEB())  // Loose the proxy definition
							if(GSFveto && PixelVeto && FSRVeto)jetPhoCollection.push_back(itpho);
					}
					if(!itpho->passSignalSelection())continue;
					if(GSFveto && PixelVeto && FSRVeto){
						if(!hasPho){
							hasPho=true;
							npassPho +=1;
							signalPho = itpho;
						}
					}

					if((!PixelVeto || !GSFveto)){
						if(itpho->isEB())proxyPhoCollection.push_back(itpho);
					}
				}
			}

			bool hasLep(false);
			std::vector<recoMuon>::iterator signalLep = Muon.begin();
			std::vector< std::vector<recoMuon>::iterator > proxyLepCollection;
			proxyLepCollection.clear();
			std::vector< std::vector<recoMuon>::iterator > fakeLepCollection;
			fakeLepCollection.clear();
			for(std::vector<recoMuon>::iterator itMu = Muon.begin(); itMu != Muon.end(); itMu++){
				if(itMu->getPt() < 25)continue;
				if(itMu->fireSingleTrg(2) || itMu->fireSingleTrg(21) || itMu->fireSingleTrg(22)){  //Mu+gamma HLT
        //if(itMu->fireSingleTrg(1) || itMu->fireSingleTrg(19)){
					if(itMu->isFakeProxy())fakeLepCollection.push_back(itMu);
					if(itMu->passSignalSelection()){
						proxyLepCollection.push_back(itMu);
						if(!hasLep && hasPho){
							hasLep=true; 
							npassLep +=1;
							signalLep = itMu;
						}
					}
				}
			}

			if(hasPho && hasLep && signalPho->isEB()){
				double dRlepphoton = DeltaR(signalPho->getEta(), signalPho->getPhi(), signalLep->getEta(), signalLep->getPhi());
				if(dRlepphoton > 0.3){
					npassdR+=1;
					if(passMETFilter(METFilter)){ 
						npassMETFilter +=1;
						if(fabs((signalPho->getCalibP4()+signalLep->getP4()).M() - 91.188) > 10.0)npassZ+=1;

						float deltaPhi = DeltaPhi(signalLep->getPhi(), METPhi);
						float MT = sqrt(2*MET*signalLep->getPt()*(1-std::cos(deltaPhi)));
						float ThreeBodyMass = sqrt(2*MET*(signalPho->getP4()+ signalLep->getP4()).Pt()*(1-std::cos(DeltaR(0, (signalPho->getP4()+signalLep->getP4()).Phi(), 0, METPhi))));

						phoEt = signalPho->getCalibEt();
						phoEta= signalPho->getEta();
						phoPhi= signalPho->getPhi();
						lepPt = signalLep->getPt();
						lepEta= signalLep->getEta();
						lepPhi= signalLep->getPhi();
						sigMT = MT;
						sigMET= MET;
						sigMETPhi = METPhi;
						dPhiLepMET = deltaPhi;
						threeMass = ThreeBodyMass;
						nVertex = nVtx;
						dRPhoLep= dRlepphoton;

						nJet = 0;
						HT = 0;
						for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
							if(!itJet->passSignalSelection())continue;
							if(DeltaR(itJet->getEta(), itJet->getPhi(), signalPho->getEta(),signalPho->getPhi()) <= 0.4)continue;	
							if(DeltaR(itJet->getEta(), itJet->getPhi(), signalLep->getEta(),signalLep->getPhi()) <= 0.4)continue;
							nJet += 1;
							HT += itJet->getPt();
						}	
						sigtree->Fill();

					}//MET Filter
				}//dR Filter
			}//Candidate Filter
		 
		for(unsigned ip(0); ip < proxyPhoCollection.size(); ip++){
		for(unsigned ie(0); ie < proxyLepCollection.size(); ie++){
			std::vector<recoPhoton>::iterator proxyPho = proxyPhoCollection[ip];
			std::vector<recoMuon>::iterator proxyMuon = proxyLepCollection[ie];
			double dRlepphoton = DeltaR(proxyPho->getEta(), proxyPho->getPhi(), proxyMuon->getEta(), proxyMuon->getPhi());
			if(dRlepphoton>0.3){
			if(passMETFilter(METFilter)){

				float proxy_deltaPhi = DeltaPhi(proxyMuon->getPhi(), METPhi);
				float proxy_MT = sqrt(2*MET*proxyMuon->getPt()*(1-std::cos(proxy_deltaPhi)));
				float proxy_ThreeBodyMass = sqrt(2*MET*(proxyPho->getP4()+ proxyMuon->getP4()).Pt()*(1-std::cos(DeltaR(0, (proxyPho->getP4()+proxyMuon->getP4()).Phi(), 0, METPhi))));
				proxyphoEt = proxyPho->getCalibEt();
				proxyphoEta= proxyPho->getEta();
				proxyphoPhi= proxyPho->getPhi();
				proxylepPt = proxyMuon->getPt();
				proxylepEta= proxyMuon->getEta();
				proxylepPhi= proxyMuon->getPhi();
				proxysigMT = proxy_MT;
				proxysigMET= MET;
				proxysigMETPhi = METPhi;
				proxydPhiLepMET = proxy_deltaPhi; 
				proxythreeMass = proxy_ThreeBodyMass;
				proxynVertex = nVtx; 
				proxydRPhoLep= dRlepphoton;
				proxysigMETJESup = MET_T1JESUp;
				proxysigMETJESdo = MET_T1JESDo;
				proxysigMETJERup = MET_T1JERUp;
				proxysigMETJERdo = MET_T1JERDo;
				float proxydPhiLepMETJESup = DeltaPhi(proxyMuon->getPhi(), METPhi_T1JESUp);
				float proxydPhiLepMETJESdo = DeltaPhi(proxyMuon->getPhi(), METPhi_T1JESDo);
				float proxydPhiLepMETJERup = proxy_deltaPhi;
				float proxydPhiLepMETJERdo = proxy_deltaPhi;
				proxysigMTJESup = sqrt(2*MET_T1JESUp*proxyMuon->getPt()*(1-std::cos(proxydPhiLepMETJESup)));
				proxysigMTJESdo = sqrt(2*MET_T1JESDo*proxyMuon->getPt()*(1-std::cos(proxydPhiLepMETJESdo)));
				proxysigMTJERup = sqrt(2*MET_T1JERUp*proxyMuon->getPt()*(1-std::cos(proxydPhiLepMETJERup)));
				proxysigMTJERdo = sqrt(2*MET_T1JERDo*proxyMuon->getPt()*(1-std::cos(proxydPhiLepMETJERdo)));

				proxynJet = 0;
				proxyHT = 0;
				proxyHTJESup = 0;
				proxyHTJESdo = 0;
				for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
					if(!itJet->passSignalSelection())continue;
					if(DeltaR(itJet->getEta(), itJet->getPhi(), proxyPho->getEta(),proxyPho->getPhi()) <= 0.4)continue;	
					if(DeltaR(itJet->getEta(), itJet->getPhi(), proxyMuon->getEta(),proxyMuon->getPhi()) <= 0.4)continue;
					proxynJet += 1;
					proxyHT += itJet->getPt();
					proxyHTJESup += itJet->getPt()*(1+itJet->getPtUnc());
					proxyHTJESdo += itJet->getPt()*(1-itJet->getPtUnc());
				}

				proxytree->Fill();

			}//MET Filter
			}//dR filter
		}// loop on ele collection
		} // loop on pho collection

		if(!hasPho){ 
		for(unsigned ip(0); ip < jetPhoCollection.size(); ip++){
		for(unsigned ie(0); ie < proxyLepCollection.size(); ie++){
			std::vector<recoPhoton>::iterator jetPho = jetPhoCollection[ip];
			std::vector<recoMuon>::iterator jetMuon = proxyLepCollection[ie];
			double dRlepphoton = DeltaR(jetPho->getEta(), jetPho->getPhi(), jetMuon->getEta(), jetMuon->getPhi());
			if(dRlepphoton>0.3){
			if(passMETFilter(METFilter)){

				float jet_deltaPhi = DeltaPhi(jetMuon->getPhi(), METPhi);
				float jet_MT = sqrt(2*MET*jetMuon->getPt()*(1-std::cos(jet_deltaPhi)));
				float jet_ThreeBodyMass = sqrt(2*MET*(jetPho->getP4()+ jetMuon->getP4()).Pt()*(1-std::cos(DeltaR(0, (jetPho->getP4()+jetMuon->getP4()).Phi(), 0, METPhi))));
				jetphoEt = jetPho->getCalibEt();
				jetphoEta= jetPho->getEta();
				jetphoPhi= jetPho->getPhi();
				jetlepPt = jetMuon->getPt();
				jetlepEta= jetMuon->getEta();
				jetlepPhi= jetMuon->getPhi();
				jetsigMT = jet_MT;
				jetsigMET= MET;
				jetsigMETPhi = METPhi;
				jetdPhiLepMET = jet_deltaPhi; 
				jetthreeMass = jet_ThreeBodyMass;
				jetnVertex = nVtx; 
				jetdRPhoLep= dRlepphoton;
				jetsigMETJESup = MET_T1JESUp;
				jetsigMETJESdo = MET_T1JESDo;
				jetsigMETJERup = MET_T1JERUp;
				jetsigMETJERdo = MET_T1JERDo;
				float jetdPhiLepMETJESup = DeltaPhi(jetMuon->getPhi(), METPhi_T1JESUp);
				float jetdPhiLepMETJESdo = DeltaPhi(jetMuon->getPhi(), METPhi_T1JESDo);
				float jetdPhiLepMETJERup = jet_deltaPhi;
				float jetdPhiLepMETJERdo = jet_deltaPhi;
				jetsigMTJESup = sqrt(2*MET_T1JESUp*jetMuon->getPt()*(1-std::cos(jetdPhiLepMETJESup)));
				jetsigMTJESdo = sqrt(2*MET_T1JESDo*jetMuon->getPt()*(1-std::cos(jetdPhiLepMETJESdo)));
				jetsigMTJERup = sqrt(2*MET_T1JERUp*jetMuon->getPt()*(1-std::cos(jetdPhiLepMETJERup)));
				jetsigMTJERdo = sqrt(2*MET_T1JERDo*jetMuon->getPt()*(1-std::cos(jetdPhiLepMETJERdo)));

				jetnJet = 0;
				jetHT = 0;
				jetHTJESup = 0;
				jetHTJESdo = 0;
				for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
					if(!itJet->passSignalSelection())continue;
					if(DeltaR(itJet->getEta(), itJet->getPhi(), jetPho->getEta(),jetPho->getPhi()) <= 0.4)continue;	
					if(DeltaR(itJet->getEta(), itJet->getPhi(), jetMuon->getEta(),jetMuon->getPhi()) <= 0.4)continue;
					jetnJet += 1;
					jetHT += itJet->getPt();
					jetHTJESup += itJet->getPt()*(1+itJet->getPtUnc());
					jetHTJESdo += itJet->getPt()*(1-itJet->getPtUnc());
				}
				jettree->Fill();

			}//MET Filter
			}//dR filter
		}// loop on ele collection
		} // loop on pho collection
		}

		if(hasPho && !hasLep && signalPho->isEB()){
			std::vector<recoPhoton>::iterator fakeLepPho = signalPho;
			for(unsigned ip(0); ip < fakeLepCollection.size(); ip++){
				std::vector<recoMuon>::iterator fakeMu = fakeLepCollection[ip];
				double dRlepphoton = DeltaR(fakeLepPho->getEta(), fakeLepPho->getPhi(), fakeMu->getEta(), fakeMu->getPhi());
				if(dRlepphoton>0.3){
					if(passMETFilter(METFilter)){

						float fakeLep_deltaPhi = DeltaPhi(fakeMu->getPhi(), METPhi);
						float fakeLep_MT = sqrt(2*MET*fakeMu->getPt()*(1-std::cos(fakeLep_deltaPhi)));
						float fake_ThreeBodyMass = sqrt(2*MET*(fakeLepPho->getP4()+ fakeMu->getP4()).Pt()*(1-std::cos(DeltaR(0, (fakeLepPho->getP4()+fakeMu->getP4()).Phi(), 0, METPhi))));
						fakeLepphoEt = fakeLepPho->getCalibEt();
						fakeLepphoEta= fakeLepPho->getEta();
						fakeLepphoPhi= fakeLepPho->getPhi();
						fakeLepPt = fakeMu->getPt();
						fakeLepEta= fakeMu->getEta();
						fakeLepPhi= fakeMu->getPhi();
						fakeLepsigMT = fakeLep_MT;
						fakeLepsigMET= MET;
						fakeLepsigMETPhi = METPhi;
						fakeLepdPhiLepMET = fakeLep_deltaPhi; 
						fakethreeMass = fake_ThreeBodyMass;
						fakeLepnVertex = nVtx; 
						fakeLepdRPhoLep= dRlepphoton;
						fakeLepsigMETJESup = MET_T1JESUp;
						fakeLepsigMETJESdo = MET_T1JESDo;
						fakeLepsigMETJERup = MET_T1JERUp;
						fakeLepsigMETJERdo = MET_T1JERDo;
						float fakeLepdPhiLepMETJESup = DeltaPhi(fakeMu->getPhi(), METPhi_T1JESUp);
						float fakeLepdPhiLepMETJESdo = DeltaPhi(fakeMu->getPhi(), METPhi_T1JESDo);
						float fakeLepdPhiLepMETJERup = fakeLep_deltaPhi;
						float fakeLepdPhiLepMETJERdo = fakeLep_deltaPhi;
						fakeLepsigMTJESup = sqrt(2*MET_T1JESUp*fakeMu->getPt()*(1-std::cos(fakeLepdPhiLepMETJESup)));
						fakeLepsigMTJESdo = sqrt(2*MET_T1JESDo*fakeMu->getPt()*(1-std::cos(fakeLepdPhiLepMETJESdo)));
						fakeLepsigMTJERup = sqrt(2*MET_T1JERUp*fakeMu->getPt()*(1-std::cos(fakeLepdPhiLepMETJERup)));
						fakeLepsigMTJERdo = sqrt(2*MET_T1JERDo*fakeMu->getPt()*(1-std::cos(fakeLepdPhiLepMETJERdo)));

						fakeLepMiniIso = fakeMu->getMiniIso();
						if(fakeMu->isMedium())fakeLepIsMedium = 1;
						else fakeLepIsMedium = 0; 
						//fakeLepMiniIso = fakeMu->getRelIso();
						fakeLepnJet = 0;
						fakeLepHT = 0;
						fakeLepHTJESup = 0;
						fakeLepHTJESdo = 0;
						for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
							if(!itJet->passSignalSelection())continue;
							if(DeltaR(itJet->getEta(), itJet->getPhi(), fakeLepPho->getEta(), fakeLepPho->getPhi()) <= 0.4)continue;	
							if(DeltaR(itJet->getEta(), itJet->getPhi(), fakeMu->getEta(),  fakeMu->getPhi()) <= 0.4)continue;
							fakeLepnJet += 1;
							fakeLepHT += itJet->getPt();
							fakeLepHTJESup += itJet->getPt()*(1+itJet->getPtUnc());
							fakeLepHTJESdo += itJet->getPt()*(1-itJet->getPtUnc());
						}	
						fakeLeptree->Fill();
					}//MET Filter
					}//dR filter
				} // loop on pho collection
			}

			/*********  ZG tree************/ 
			if(hasPho && hasLep && signalPho->isEB() && raw.nMu >= 2 ){
				double dRlepphoton = DeltaR(signalPho->getEta(), signalPho->getPhi(), signalLep->getEta(), signalLep->getPhi());
				if(dRlepphoton > 0.3){
					if(passMETFilter(METFilter)){ 
						
						bool foundZG(false); std::vector<recoMuon>::iterator trailLep=Muon.begin();	
						for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++){
							if(im == signalLep)continue;
							if(!im->isMedium())continue;
							double llmass = (im->getP4()+signalLep->getP4()).M();
							double llgmass = (im->getP4()+signalLep->getP4()+signalPho->getP4()).M();
							if( (80 < llmass && llmass < 100) || ( 80 < llgmass && llgmass < 100)){
								foundZG = true;
								trailLep = im;
							}
						}	
					
						if(foundZG){	
							float deltaPhi = DeltaPhi(signalLep->getPhi(), METPhi);
							float MT = sqrt(2*MET*signalLep->getPt()*(1-std::cos(deltaPhi)));
							float ThreeBodyMass = (trailLep->getP4()+signalLep->getP4()+signalPho->getP4()).M(); 

							ZcalibthreeMass = (trailLep->getP4()+signalLep->getP4()+signalPho->getCalibP4()).M();
							ZphoEt = signalPho->getCalibEt();
							ZphoEta= signalPho->getEta();
							ZphoPhi= signalPho->getPhi();
							ZphoChIso= signalPho->getChIso();
							ZphoNeuIso=signalPho->getNeuIso();
							ZlepPt = signalLep->getPt();
							ZlepEta= signalLep->getEta();
							ZlepPhi= signalLep->getPhi();
							ZtrailPt = trailLep->getPt();
							ZtrailEta = trailLep->getEta();
							ZtrailPhi = trailLep->getPhi();
							ZsigMT = MT;
							ZsigMET= MET;
							ZsigMETPhi = METPhi;
							ZdPhiLepMET = deltaPhi;
							ZthreeMass = ThreeBodyMass;
							ZdilepMass = (trailLep->getP4()+signalLep->getP4()).M();
							ZnVertex = nVtx;
							ZdRPhoLep= dRlepphoton;
							ZdRPhoTrail=DeltaR(signalPho->getEta(), signalPho->getPhi(), trailLep->getEta(), trailLep->getPhi());
							ZnMu = raw.nMu;
							ZnGoodMu = nGoodMu;

							ZnJet = 0;
							ZHT = 0;
							for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
								if(!itJet->passSignalSelection())continue;
								if(DeltaR(itJet->getEta(), itJet->getPhi(), signalPho->getEta(),signalPho->getPhi()) <= 0.4)continue;	
								if(DeltaR(itJet->getEta(), itJet->getPhi(), signalLep->getEta(),signalLep->getPhi()) <= 0.4)continue;
								ZnJet += 1;
								ZHT += itJet->getPt();
							}

						  Z_mcPID.clear();
						  Z_mcEta.clear();
						  Z_mcPhi.clear();
						  Z_mcPt.clear();
						  Z_mcMomPID.clear();
						  Z_mcGMomPID.clear();
							Z_mcStatus.clear();
						  for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
						 	  if(itMC->getEt() < 1.0)continue;
							  float mcdRmu = DeltaR(signalLep->getEta(), signalLep->getPhi(), itMC->getEta(), itMC->getPhi());
						 	  float mcdR = DeltaR(signalPho->getEta(), signalPho->getPhi(), itMC->getEta(), itMC->getPhi());
						 	  if(mcdR < 0.3 || mcdRmu < 0.3){
						 		  Z_mcPID.push_back(itMC->getPID());
						 		  Z_mcMomPID.push_back(itMC->getMomPID());
						 		  Z_mcGMomPID.push_back(itMC->getGMomPID());
						 		  Z_mcEta.push_back(itMC->getEta());
						 		  Z_mcPhi.push_back(itMC->getPhi());
						 		  Z_mcPt.push_back(itMC->getEt());
									Z_mcStatus.push_back(itMC->getStatus());
						 	  }
						  }
							
							Ztree->Fill();
						}	

					}//MET Filter
				}//dR Filter
			}//Candidate Filter


		if(raw.nMu >= 2 ){ 
			for(unsigned ip(0); ip < jetPhoCollection.size(); ip++){
			for(unsigned ie(0); ie < proxyLepCollection.size(); ie++){
				std::vector<recoPhoton>::iterator jetPho = jetPhoCollection[ip];
				std::vector<recoMuon>::iterator jetMuon = proxyLepCollection[ie];
				double dRlepphoton = DeltaR(jetPho->getEta(), jetPho->getPhi(), jetMuon->getEta(), jetMuon->getPhi());
				if(dRlepphoton>0.3){
				if(passMETFilter(METFilter)){

					bool foundZG(false); std::vector<recoMuon>::iterator trailLep=Muon.begin();	
					for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++){
						if(im == jetMuon)continue;
						double llmass = (im->getP4()+jetMuon->getP4()).M();
						double llgmass = (im->getP4()+jetMuon->getP4()+jetPho->getP4()).M();
						if( (80 < llmass && llmass < 100) || ( 80 < llgmass && llgmass < 100)){
							foundZG = true;
							trailLep = im;
						}
					}	

					if(foundZG){
						float ZJet_deltaPhi = DeltaPhi(jetMuon->getPhi(), METPhi);
						float ZJet_MT = sqrt(2*MET*jetMuon->getPt()*(1-std::cos(ZJet_deltaPhi)));
						float ZJet_ThreeBodyMass = (trailLep->getP4()+jetMuon->getP4()+jetPho->getP4()).M(); 
						ZJetphoEt = jetPho->getCalibEt();
						ZJetphoEta= jetPho->getEta();
						ZJetphoPhi= jetPho->getPhi();
						ZJetphoSigma = jetPho->getSigma();
						ZJetphoChIso = jetPho->getChIso();
						ZJetlepPt = jetMuon->getPt();
						ZJetlepEta= jetMuon->getEta();
						ZJetlepPhi= jetMuon->getPhi();
						ZJetsigMT = ZJet_MT;
						ZJetsigMET= MET;
						ZJetsigMETPhi = METPhi;
						ZJetdPhiLepMET = ZJet_deltaPhi; 
						ZJetthreeMass = ZJet_ThreeBodyMass;
						ZJetdilepMass = (trailLep->getP4()+jetMuon->getP4()).M();
						ZJetnVertex = nVtx; 
						ZJetdRPhoLep= dRlepphoton;

						ZJetnJet = 0;
						ZJetHT = 0;
						for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
							if(!itJet->passSignalSelection())continue;
							if(DeltaR(itJet->getEta(), itJet->getPhi(), jetPho->getEta(),jetPho->getPhi()) <= 0.4)continue;	
							if(DeltaR(itJet->getEta(), itJet->getPhi(), jetMuon->getEta(),jetMuon->getPhi()) <= 0.4)continue;
							ZJetnJet += 1;
							ZJetHT += itJet->getPt();
						}
						ZJettree->Fill();
					}

			}//MET Filter
			}//dR filter
		}// loop on ele collection
		} // loop on pho collection
		}
	}//loop on  events

p_eventcount->Fill("Total",nTotal);
p_eventcount->Fill("passHLT",npassHLT);
p_eventcount->Fill("passPho",npassPho);
p_eventcount->Fill("passMuon",npassLep);
p_eventcount->Fill("passdR",npassdR);
p_eventcount->Fill("passMETFilter",npassMETFilter);
p_eventcount->Fill("passZ",npassZ);

outputfile->Write();
logfile.close();
}


