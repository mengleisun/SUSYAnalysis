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

#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_mcData.h"
#include "../../../include/analysis_tools.h"
#include "../../../include/analysis_jet.h"

void analysis_mgMC(){//main  

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

  char outputname[100] = "/uscms_data/d3/mengleis/FullStatusOct/resTree_mgsignal_WZ.root";
  ofstream logfile;
  logfile.open("/uscms_data/d3/mengleis/FullStatusOct/resTree_mgsignal_WZ.log"); 

  logfile << "analysis_mg()" << std::endl;
  logfile << "miniIso; one lepton for fakephoton background" << std::endl;
	
  RunType datatype(MCMuonEG);
	bool  isMC(false);
	if(datatype == MC || datatype == MCDoubleEG || datatype == MCMuonEG||  datatype == MCSingleElectron || datatype == MCSingleMuon||  datatype == MCDoubleMuon || datatype == MCMET)isMC=true;
  TChain* es = new TChain("ggNtuplizer/EventTree");
	//es->Add("root://cmseos.fnal.gov///store/group/lpcsusystealth/ggNtuple_leppho/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1.root");
//	es->Add("root://cmseos.fnal.gov//store/user/msun/MCSummer16/DYJetsToLL_M-50_nlo.root");
//	es->Add("root://cmseos.fnal.gov//store/user/msun/MCSummer16/TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8.root");
//	es->Add("root://cmseos.fnal.gov//store/user/msun/MCSummer16/TTJets_madgraphMLM_RunIISummer16MiniAODv2-TrancheIV_v6.root");
//	es->Add("root://cmseos.fnal.gov//store/user/msun/MCSummer16/WW_TuneCUETP8M1_13TeV-pythia8.root");
	es->Add("root://cmseos.fnal.gov//store/user/msun/MCSummer16/WZ_TuneCUETP8M1_13TeV-pythia8.root");
//	es->Add("root://cmseos.fnal.gov//store/user/msun/MCSummer16/WJetsToLNu_RunIISummer16MiniAODv2-TrancheIV_v6-ext2-v1.root");
//	logfile << "root://cmseos.fnal.gov//store/user/msun/MCSummer16/WW_TuneCUETP8M1_13TeV-pythia8.root" << std::endl;

  const unsigned nEvts = es->GetEntries(); 
  logfile << "Total event: " << nEvts << std::endl;
  std::cout << "Total event: " << nEvts << std::endl;
  logfile << "Output file: " << outputname << std::endl;

	int nTotal(0),npassHLT(0), npassPho(0), npassLep(0), npassdR(0), npassZ(0), npassMETFilter(0);

  TFile *outputfile = TFile::Open(outputname,"RECREATE");
  outputfile->cd();

  int mcType = MCType::WZ;
  float crosssection = MC_XS[mcType];
  float ntotalevent = es->GetEntries();

	int nBJet(0);
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
  int   nVertex(0);
  float dRPhoLep(0);
  float HT(0);
  int   nJet(0);
	float trailPt(0);
	float trailEta(0);
	float trailPhi(0);
  std::vector<int>   mcPID;
  std::vector<float> mcEta;
  std::vector<float> mcPhi;
  std::vector<float> mcPt;
  std::vector<int>   mcMomPID;
  std::vector<int>   mcGMomPID;

	sigtree->Branch("crosssection",&crosssection);
	sigtree->Branch("ntotalevent", &ntotalevent);
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
  sigtree->Branch("nVertex",   &nVertex);
  sigtree->Branch("dRPhoLep",  &dRPhoLep);
  sigtree->Branch("HT",        &HT);
  sigtree->Branch("nJet",      &nJet);
  sigtree->Branch("nBJet",     &nBJet);
	sigtree->Branch("trailPt",   &trailPt);
	sigtree->Branch("trailEta",  &trailEta);
	sigtree->Branch("trailPhi",  &trailPhi);
	if(isMC){
  	sigtree->Branch("mcPID",     &mcPID);
  	sigtree->Branch("mcEta",     &mcEta);
  	sigtree->Branch("mcPhi",     &mcPhi);
  	sigtree->Branch("mcPt",      &mcPt);
  	sigtree->Branch("mcMomPID",  &mcMomPID);
  	sigtree->Branch("mcGMomPID", &mcGMomPID);
	}

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
  int   proxynVertex(0);
  float proxydRPhoLep(0);
  float proxyHT(0);
  int   proxynJet(0);
	float proxytrailPt(0);
	float proxytrailEta(0);
	float proxytrailPhi(0);
  
	proxytree->Branch("crosssection",&crosssection);
	proxytree->Branch("ntotalevent", &ntotalevent);
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
  proxytree->Branch("nVertex",   &proxynVertex);
  proxytree->Branch("dRPhoLep",  &proxydRPhoLep);
  proxytree->Branch("HT",        &proxyHT);
  proxytree->Branch("nJet",      &proxynJet);
  proxytree->Branch("nBJet",     &nBJet);
	proxytree->Branch("trailPt",   &proxytrailPt);
	proxytree->Branch("trailEta",  &proxytrailEta);
	proxytree->Branch("trailPhi",  &proxytrailPhi);

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
  int   jetnVertex(0);
  float jetdRPhoLep(0);
  float jetHT(0);
  int   jetnJet(0);
	float jettrailPt(0);
	float jettrailEta(0);
	float jettrailPhi(0);
  
	jettree->Branch("crosssection",&crosssection);
	jettree->Branch("ntotalevent", &ntotalevent);
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
  jettree->Branch("nVertex",   &jetnVertex);
  jettree->Branch("dRPhoLep",  &jetdRPhoLep);
  jettree->Branch("HT",        &jetHT);
  jettree->Branch("nJet",      &jetnJet);
  jettree->Branch("nBJet",     &nBJet);
	jettree->Branch("trailPt",   &jettrailPt);
	jettree->Branch("trailEta",  &jettrailEta);
	jettree->Branch("trailPhi",  &jettrailPhi);
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
  int   fakeLepnVertex(0);
  float fakeLepdRPhoLep(0);
  float fakeLepHT(0);
  int   fakeLepnJet(0);
  
	fakeLeptree->Branch("crosssection",&crosssection);
	fakeLeptree->Branch("ntotalevent", &ntotalevent);
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
  fakeLeptree->Branch("nVertex",   &fakeLepnVertex);
  fakeLeptree->Branch("dRPhoLep",  &fakeLepdRPhoLep);
  fakeLeptree->Branch("HT",        &fakeLepHT);
  fakeLeptree->Branch("nJet",      &fakeLepnJet);
  fakeLeptree->Branch("nBJet",     &nBJet);

//*************** for jet-photon fake rate ***********************//
	TTree *hadrontree = new TTree("hadronTree","hadronTree");
	float hadron_phoEt(0);
	float hadron_phoEta(0);
	float hadron_phoPhi(0);
	float hadron_phoSigma(0);
	float hadron_phoChIso(0);
  std::vector<float> hadron_eleproxyEt;
	std::vector<float> hadron_eleproxyEta;
	std::vector<float> hadron_eleproxyPhi;
	std::vector<float> hadron_eleproxySigma;
	std::vector<float> hadron_eleproxyChIso;
	std::vector<int> hadron_eleproxynVertex;
	float hadron_sigMT(0);
	float hadron_sigMET(0);
	float hadron_sigMETPhi(0);
	float hadron_dPhiLepMET(0);
	int   hadron_nVertex(0);
	float hadron_HT(0);
	float hadron_nJet(0);
	std::vector<int> hadron_mcPID;
	std::vector<float> hadron_mcEta;
	std::vector<float> hadron_mcPhi;
	std::vector<float> hadron_mcPt;
	std::vector<int> hadron_mcMomPID;
	std::vector<int> hadron_mcGMomPID;

	hadrontree->Branch("crosssection",&crosssection);
	hadrontree->Branch("ntotalevent", &ntotalevent);
	hadrontree->Branch("phoEt",     &hadron_phoEt);
	hadrontree->Branch("phoEta",    &hadron_phoEta);
	hadrontree->Branch("phoPhi",    &hadron_phoPhi);
	hadrontree->Branch("phoSigma",  &hadron_phoSigma);
	hadrontree->Branch("phoChIso",  &hadron_phoChIso);
	hadrontree->Branch("eleproxyEt",&hadron_eleproxyEt);
	hadrontree->Branch("eleproxyEta",&hadron_eleproxyEta);
	hadrontree->Branch("eleproxyPhi",&hadron_eleproxyPhi);
	hadrontree->Branch("eleproxySigma", &hadron_eleproxySigma);
	hadrontree->Branch("eleproxyChIso", &hadron_eleproxyChIso);
	hadrontree->Branch("eleproxynVertex",&hadron_eleproxynVertex);
	hadrontree->Branch("sigMT",     &hadron_sigMT);
	hadrontree->Branch("sigMET",    &hadron_sigMET);
	hadrontree->Branch("sigMETPhi", &hadron_sigMETPhi);
	hadrontree->Branch("dPhiLepMET",&hadron_dPhiLepMET);
	hadrontree->Branch("nVertex",   &hadron_nVertex);
	hadrontree->Branch("HT",        &hadron_HT);
	hadrontree->Branch("nJet",      &hadron_nJet);
	if(isMC){
		hadrontree->Branch("mcPID",     &hadron_mcPID);
		hadrontree->Branch("mcEta",     &hadron_mcEta);
		hadrontree->Branch("mcPhi",     &hadron_mcPhi);
		hadrontree->Branch("mcPt",      &hadron_mcPt);
		hadrontree->Branch("mcMomPID",  &hadron_mcMomPID);
		hadrontree->Branch("mcGMomPID", &hadron_mcGMomPID);
	}
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
			if(isMC)for(int iMC(0); iMC < raw.nMC; iMC++){MCData.push_back(mcData(raw, iMC));}
			for(int iPho(0); iPho < raw.nPho; iPho++){Photon.push_back(recoPhoton(raw, iPho));}
			for(int iMu(0); iMu < raw.nMu; iMu++){Muon.push_back(recoMuon(raw, iMu));}
			for(int iEle(0); iEle < raw.nEle; iEle++){Ele.push_back(recoEle(raw, iEle));}
			for(int iJet(0); iJet < raw.nJet; iJet++){JetCollection.push_back(recoJet(raw, iJet));}
			MET = raw.pfMET;
			METPhi = raw.pfMETPhi;
			METFilter = raw.metFilters;
			nVtx = raw.nVtx;
			run=raw.run;
			event=raw.event;
			lumis=raw.lumis;

			nTotal+=1;
			if(!raw.passHLT())continue;
			if(raw.nGoodVtx < 1)continue;
			npassHLT+=1;

			if(raw.nMu < 1 || raw.nPho <1)continue;

			nBJet = 0;
			for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
				if(itJet->getPt() < 20)continue;
				if(itJet->isBJet())nBJet+=1;
			}
			/******************************************************************************************************************************************************************************/
			/***********************************                                  Select Photon                                              **********************************************/
			bool hasPho(false);
			std::vector<recoPhoton>::iterator signalPho = Photon.begin();
			std::vector< std::vector<recoPhoton>::iterator >  proxyPhoCollection;
			proxyPhoCollection.clear();
			std::vector< std::vector<recoPhoton>::iterator >  jetPhoCollection;
			jetPhoCollection.clear();
			// ******* for jet-photon fake rate***************//
			bool hasHadronPho(false);
			std::vector<recoPhoton>::iterator hadronPho = Photon.begin();
			std::vector< std::vector<recoPhoton>::iterator >  hadeleproxyPhoCollection;
			hadeleproxyPhoCollection.clear();

			for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
				if(itpho->getR9() < 0.5)continue;
				if(!itpho->passHLTSelection())continue;
				if(!itpho->passBasicSelection())continue;
				bool passSigma = itpho->passSigma(1);
				bool passChIso = itpho->passChIso(1);
				bool PixelVeto = itpho->PixelSeed()==0? true: false;
				bool GSFveto(true);
				bool FSRVeto(true);
				bool muonFSRVeto(true);
				for(std::vector<recoEle>::iterator ie = Ele.begin(); ie != Ele.end(); ie++){
					if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) <= 0.02)GSFveto = false;
					if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.3)FSRVeto=false;
				}
				for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++)
					if(DeltaR(itpho->getEta(), itpho->getPhi(), im->getEta(), im->getPhi()) < 0.3 && im->getEt()>2.0){FSRVeto=false; muonFSRVeto=false;}

				// ******** very loose, before sigma and isolation cut.  For jet-photon fake rate, and hadron proxy ************//	
				if(GSFveto && PixelVeto && FSRVeto){
					if(!hasHadronPho){
						hasHadronPho = true;
						hadronPho = itpho;
					}

					if(!passSigma || !passChIso){
						if( (itpho->getSigma()< 0.02 && itpho->isEB()) || (itpho->getSigma()< 0.04 && itpho->isEE()) )jetPhoCollection.push_back(itpho);
					}
				}
				else if( (!GSFveto || !PixelVeto) )hadeleproxyPhoCollection.push_back(itpho);

        // ****************  standard ID ************************************//
				if(!itpho->passSignalSelection())continue;
				if(GSFveto && PixelVeto && FSRVeto){
					if(!hasPho){
						hasPho=true;
						npassPho +=1;
						signalPho = itpho;
					}
				}

				if((!PixelVeto || !GSFveto) && muonFSRVeto){
						proxyPhoCollection.push_back(itpho);
				}
			}

			bool hasLep(false);
			std::vector<recoMuon>::iterator signalLep = Muon.begin();
			bool hasTrail(false);
			std::vector<recoMuon>::iterator trailLep = Muon.begin();
			std::vector< std::vector<recoMuon>::iterator > proxyLepCollection;
			proxyLepCollection.clear();
			std::vector< std::vector<recoMuon>::iterator > fakeLepCollection;
			fakeLepCollection.clear();
			std::vector< std::vector<recoMuon>::iterator > miniisoLep;
			miniisoLep.clear();
			for(std::vector<recoMuon>::iterator itMu = Muon.begin(); itMu != Muon.end(); itMu++){
				if(itMu->isMedium() && itMu->getPt() > 15 && itMu->getMiniIso() < 0.2)miniisoLep.push_back(itMu);
				if(itMu->getPt() < 25)continue;
				if(!itMu->passHLTSelection())continue;
				if(itMu->isFakeProxy())fakeLepCollection.push_back(itMu);
				if(itMu->passSignalSelection()){
					if(proxyLepCollection.size() == 0)proxyLepCollection.push_back(itMu);
					if(hasLep && !hasTrail){
						hasTrail = true;
						trailLep = itMu;
					}
					if(!hasLep){
						hasLep=true; 
						npassLep +=1;
						signalLep = itMu;
					}
				}
			}

			if(hasPho && hasLep){
				double dRlepphoton = DeltaR(signalPho->getEta(), signalPho->getPhi(), signalLep->getEta(), signalLep->getPhi());
				if(dRlepphoton > 0.8){
					npassdR+=1;
					if(raw.passMETFilter(METFilter)){ 
						npassMETFilter +=1;
						if(fabs((signalPho->getCalibP4()+signalLep->getP4()).M() - 91.188) > 10.0)npassZ+=1;

						float deltaPhi = DeltaPhi(signalLep->getPhi(), METPhi);
						float MT = sqrt(2*MET*signalLep->getPt()*(1-std::cos(deltaPhi)));

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
						nVertex = nVtx;
						dRPhoLep= dRlepphoton;
						if(hasTrail){
							trailPt = trailLep->getPt();
							trailEta = trailLep->getEta();
							trailPhi = trailLep->getPhi();
						}
						else{
							trailPt  = 0; 
							trailEta = 0;
							trailPhi = 0;
						}

						nJet = 0;
						HT = 0;
						for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
							if(!itJet->passSignalSelection())continue;
							if(DeltaR(itJet->getEta(), itJet->getPhi(), signalPho->getEta(),signalPho->getPhi()) <= 0.4)continue;	
							if(DeltaR(itJet->getEta(), itJet->getPhi(), signalLep->getEta(),signalLep->getPhi()) <= 0.4)continue;
							nJet += 1;
							HT += itJet->getPt();
						}

            mcPID.clear();
            mcEta.clear();
            mcPhi.clear();
            mcPt.clear();
            mcMomPID.clear();
            mcGMomPID.clear();
						if(isMC){
           	 for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
           	   if(itMC->getEt() < 1.0)continue;
           	   float mcdR = DeltaR(signalPho->getEta(), signalPho->getPhi(), itMC->getEta(), itMC->getPhi());
           	   if(mcdR < 0.3){
           	     mcPID.push_back(itMC->getPID());
           	     mcMomPID.push_back(itMC->getMomPID());
           	     mcGMomPID.push_back(itMC->getGMomPID());
           	     mcEta.push_back(itMC->getEta());
           	     mcPhi.push_back(itMC->getPhi());
           	     mcPt.push_back(itMC->getEt());
           	   }
           	 }
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
			if(dRlepphoton>0.8){
			if(raw.passMETFilter(METFilter)){

				float proxy_deltaPhi = DeltaPhi(proxyMuon->getPhi(), METPhi);
				float proxy_MT = sqrt(2*MET*proxyMuon->getPt()*(1-std::cos(proxy_deltaPhi)));
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
				proxynVertex = nVtx; 
				proxydRPhoLep= dRlepphoton;
				if(hasTrail){
					proxytrailPt = trailLep->getPt();
					proxytrailEta = trailLep->getEta();
					proxytrailPhi = trailLep->getPhi();
				}
				else{
					proxytrailPt  = 0; 
					proxytrailEta = 0;
					proxytrailPhi = 0;
				}

				proxynJet = 0;
				proxyHT = 0;
				for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
					if(!itJet->passSignalSelection())continue;
					if(DeltaR(itJet->getEta(), itJet->getPhi(), proxyPho->getEta(),proxyPho->getPhi()) <= 0.4)continue;	
					if(DeltaR(itJet->getEta(), itJet->getPhi(), proxyMuon->getEta(),proxyMuon->getPhi()) <= 0.4)continue;
					proxynJet += 1;
					proxyHT += itJet->getPt();
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
			if(dRlepphoton>0.8){
			if(raw.passMETFilter(METFilter)){

				float jet_deltaPhi = DeltaPhi(jetMuon->getPhi(), METPhi);
				float jet_MT = sqrt(2*MET*jetMuon->getPt()*(1-std::cos(jet_deltaPhi)));
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
				jetnVertex = nVtx; 
				jetdRPhoLep= dRlepphoton;
				if(hasTrail){
					jettrailPt = trailLep->getPt();
					jettrailEta = trailLep->getEta();
					jettrailPhi = trailLep->getPhi();
				}
				else{
					jettrailPt  = 0; 
					jettrailEta = 0;
					jettrailPhi = 0;
				}

				jetnJet = 0;
				jetHT = 0;
				for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
					if(!itJet->passSignalSelection())continue;
					if(DeltaR(itJet->getEta(), itJet->getPhi(), jetPho->getEta(),jetPho->getPhi()) <= 0.4)continue;	
					if(DeltaR(itJet->getEta(), itJet->getPhi(), jetMuon->getEta(),jetMuon->getPhi()) <= 0.4)continue;
					jetnJet += 1;
					jetHT += itJet->getPt();
				}
				jettree->Fill();

			}//MET Filter
			}//dR filter
		}// loop on ele collection
		} // loop on pho collection
		}

		if(hasPho && !hasLep){
			std::vector<recoPhoton>::iterator fakeLepPho = signalPho;
			for(unsigned ip(0); ip < fakeLepCollection.size(); ip++){
				std::vector<recoMuon>::iterator fakeMu = fakeLepCollection[ip];
				double dRlepphoton = DeltaR(fakeLepPho->getEta(), fakeLepPho->getPhi(), fakeMu->getEta(), fakeMu->getPhi());
				if(dRlepphoton>0.8){
					if(raw.passMETFilter(METFilter)){

						float fakeLep_deltaPhi = DeltaPhi(fakeMu->getPhi(), METPhi);
						float fakeLep_MT = sqrt(2*MET*fakeMu->getPt()*(1-std::cos(fakeLep_deltaPhi)));
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
						fakeLepnVertex = nVtx; 
						fakeLepdRPhoLep= dRlepphoton;

						fakeLepMiniIso = fakeMu->getMiniIso();
						if(fakeMu->isMedium())fakeLepIsMedium = 1;
						else fakeLepIsMedium = 0; 
						//fakeLepMiniIso = fakeMu->getRelIso();
						fakeLepnJet = 0;
						fakeLepHT = 0;
						for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
							if(!itJet->passSignalSelection())continue;
							if(DeltaR(itJet->getEta(), itJet->getPhi(), fakeLepPho->getEta(), fakeLepPho->getPhi()) <= 0.4)continue;	
							if(DeltaR(itJet->getEta(), itJet->getPhi(), fakeMu->getEta(),  fakeMu->getPhi()) <= 0.4)continue;
							fakeLepnJet += 1;
							fakeLepHT += itJet->getPt();
						}	
						fakeLeptree->Fill();
					}//MET Filter
					}//dR filter
				} // loop on pho collection
			}

		if(hasHadronPho || hadeleproxyPhoCollection.size() > 0){
			hadron_phoEt = 0;
			hadron_phoEta = 0;
			hadron_phoPhi = 0;
			hadron_phoSigma = 0;
			hadron_phoChIso = 0;
			hadron_sigMT = 0;
			hadron_sigMET = 0;
			hadron_sigMETPhi = 0;
			hadron_dPhiLepMET = 0;
			hadron_nVertex = 0;
			hadron_HT = 0;
			hadron_nJet = 0;
  		hadron_eleproxyEt.clear();
			hadron_eleproxyEta.clear();
			hadron_eleproxyPhi.clear();
			hadron_eleproxySigma.clear();
			hadron_eleproxyChIso.clear();
			hadron_eleproxynVertex.clear();

			if(hasLep){
				double DeltaPhoLep = DeltaR(hadronPho->getEta(), hadronPho->getPhi(), signalLep->getEta(), signalLep->getPhi());
				double DoubleMass  = (hadronPho->getP4()+signalLep->getP4()).M();
				if(DeltaPhoLep > 0.8){
				if(raw.passMETFilter(METFilter)){
					float deltaPhi = DeltaPhi(signalLep->getPhi(), METPhi);
					float MT = sqrt(2*MET*signalLep->getPt()*(1-std::cos(deltaPhi)));
					hadron_phoEt = hadronPho->getCalibEt();
					hadron_phoEta = hadronPho->getEta();
					hadron_phoPhi = hadronPho->getPhi();
					hadron_phoSigma = hadronPho->getSigma();
					hadron_phoChIso = hadronPho->getChIso();
					hadron_sigMT = MT;
					hadron_sigMET = MET;
					hadron_sigMETPhi = METPhi;
					hadron_dPhiLepMET = deltaPhi;
					hadron_nVertex = nVtx;
					hadron_HT = 0;
					hadron_nJet = jetNumber;
				}}
			}
			for(unsigned ip(0); ip < hadeleproxyPhoCollection.size(); ip++){
				for(unsigned ie(0); ie < proxyLepCollection.size(); ie++){
					std::vector<recoPhoton>::iterator proxyPho = hadeleproxyPhoCollection[ip];
					std::vector<recoMuon>::iterator proxyMuon = proxyLepCollection[ie];
					double dRlepphoton = DeltaR(proxyPho->getEta(), proxyPho->getPhi(), proxyMuon->getEta(), proxyMuon->getPhi());
					if(dRlepphoton>0.8 && raw.passMETFilter(METFilter)){
						hadron_eleproxyEt.push_back(proxyPho->getCalibEt());
						hadron_eleproxyEta.push_back(proxyPho->getEta());
						hadron_eleproxyPhi.push_back(proxyPho->getPhi());
						hadron_eleproxySigma.push_back(proxyPho->getSigma());
						hadron_eleproxyChIso.push_back(proxyPho->getChIso());
						hadron_eleproxynVertex.push_back(nVtx); 
					}//dR filter
				}// loop on ele collection
			} // loop on pho collection

      hadron_mcPID.clear();
      hadron_mcEta.clear();
      hadron_mcPhi.clear();
      hadron_mcPt.clear();
      hadron_mcMomPID.clear();
      hadron_mcGMomPID.clear();
			if(isMC){
       for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
         if(itMC->getPt() < 10.0)continue;
         mcPID.push_back(itMC->getPID());
         mcMomPID.push_back(itMC->getMomPID());
         mcGMomPID.push_back(itMC->getGMomPID());
         mcEta.push_back(itMC->getEta());
         mcPhi.push_back(itMC->getPhi());
         mcPt.push_back(itMC->getEt());
       }
			}

			hadrontree->Fill();
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


