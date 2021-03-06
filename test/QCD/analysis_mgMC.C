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

void analysis_mgMC(){//main  

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

  char outputname[100] = "/uscms_data/d3/mengleis/FullStatusOct/fakelep_mgsignal_QCD.root";
  ofstream logfile;
  logfile.open("/uscms_data/d3/mengleis/FullStatusOct/fakelep_mgsignal_QCD.log"); 

  logfile << "analysis_mg()" << std::endl;
  logfile << "miniIso; one lepton for fakephoton background" << std::endl;
	
  RunType datatype(MCMuonEG);
	bool  isMC(false);
	if(datatype == MC || datatype == MCDoubleEG || datatype == MCMuonEG||  datatype == MCSingleElectron || datatype == MCSingleMuon||  datatype == MCDoubleMuon || datatype == MCMET)isMC=true;
  TChain* es = new TChain("ggNtuplizer/EventTree");
	//es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1.root");
	es->Add("root://cmseos.fnal.gov//store/user/msun/MCSummer16/QCD_Pt-20toInf_MuEnrichedPt15_pythia8.root");

  const unsigned nEvts = es->GetEntries(); 
  logfile << "Total event: " << nEvts << std::endl;
  std::cout << "Total event: " << nEvts << std::endl;
  logfile << "Output file: " << outputname << std::endl;

	int nTotal(0),npassHLT(0), npassPho(0), npassLep(0), npassdR(0), npassZ(0), npassMETFilter(0);

  TFile *outputfile = TFile::Open(outputname,"RECREATE");
  outputfile->cd();

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
	float threeMass(0);
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

//*********** fake lepton *********************//
  TTree *fakeLeptree = new TTree("fakeLepTree","fakeLepTree");
  float fakeLepphoEt(0);
  float fakeLepphoEta(0);
  float fakeLepphoPhi(0);
  float fakeLepPt(0);
  float fakeLepEta(0);
  float fakeLepPhi(0);
	int   fakeLepIsStandardProxy(0);
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
  int   fakeLepnJet(0);
	std::vector<float> fakeLep_JetPt;
	std::vector<float> fakeLep_JetEta;
	std::vector<float> fakeLep_JetPhi;
  std::vector<int>   fakeLep_mcPID;
  std::vector<float> fakeLep_mcEta;
  std::vector<float> fakeLep_mcPhi;
  std::vector<float> fakeLep_mcPt;
  std::vector<int>   fakeLep_mcMomPID;
  std::vector<int>   fakeLep_mcGMomPID;
  
  
  fakeLeptree->Branch("phoEt",     &fakeLepphoEt);
  fakeLeptree->Branch("phoEta",    &fakeLepphoEta);
  fakeLeptree->Branch("phoPhi",    &fakeLepphoPhi);
  fakeLeptree->Branch("lepPt",     &fakeLepPt);
  fakeLeptree->Branch("lepEta",    &fakeLepEta);
  fakeLeptree->Branch("lepPhi",    &fakeLepPhi);
  fakeLeptree->Branch("fakeLepMiniIso",&fakeLepMiniIso);
  fakeLeptree->Branch("fakeLepIsStandardProxy", &fakeLepIsStandardProxy);
  fakeLeptree->Branch("sigMT",     &fakeLepsigMT);
  fakeLeptree->Branch("sigMET",    &fakeLepsigMET);
  fakeLeptree->Branch("sigMETPhi", &fakeLepsigMETPhi);
  fakeLeptree->Branch("dPhiLepMET",&fakeLepdPhiLepMET);
	fakeLeptree->Branch("threeMass", &fakethreeMass);
  fakeLeptree->Branch("nVertex",   &fakeLepnVertex);
  fakeLeptree->Branch("dRPhoLep",  &fakeLepdRPhoLep);
  fakeLeptree->Branch("HT",        &fakeLepHT);
  fakeLeptree->Branch("nJet",      &fakeLepnJet);
  fakeLeptree->Branch("nBJet",     &nBJet);
	fakeLeptree->Branch("JetPt",     &fakeLep_JetPt);
	fakeLeptree->Branch("JetEta",    &fakeLep_JetEta);
	fakeLeptree->Branch("JetPhi",    &fakeLep_JetPhi);
	if(isMC){
  	fakeLeptree->Branch("mcPID",     &fakeLep_mcPID);
  	fakeLeptree->Branch("mcEta",     &fakeLep_mcEta);
  	fakeLeptree->Branch("mcPhi",     &fakeLep_mcPhi);
  	fakeLeptree->Branch("mcPt",      &fakeLep_mcPt);
  	fakeLeptree->Branch("mcMomPID",  &fakeLep_mcMomPID);
  	fakeLeptree->Branch("mcGMomPID", &fakeLep_mcGMomPID);
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

			bool hasPho(false);
			std::vector<recoJet>::iterator signalPho = JetCollection.begin();
			for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
				if(!itJet->passSignalSelection())continue;
				if(itJet->getPt() < 35 || fabs(itJet->getEta())>1.442)continue;
				if(!hasPho){
					signalPho = itJet;
					hasPho = true;
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
						if(fabs((signalPho->getP4()+signalLep->getP4()).M() - 91.188) > 10.0)npassZ+=1;

						float deltaPhi = DeltaPhi(signalLep->getPhi(), METPhi);
						float MT = sqrt(2*MET*signalLep->getPt()*(1-std::cos(deltaPhi)));
						float ThreeBodyMass = sqrt(2*MET*(signalPho->getP4()+ signalLep->getP4()).Pt()*(1-std::cos(DeltaR(0, (signalPho->getP4()+signalLep->getP4()).Phi(), 0, METPhi))));

						phoEt = signalPho->getEt();
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
		 


		if(hasPho && !hasLep){
			std::vector<recoJet>::iterator fakeLepPho = signalPho;
			for(unsigned ip(0); ip < fakeLepCollection.size(); ip++){
				std::vector<recoMuon>::iterator fakeMu = fakeLepCollection[ip];
				double dRlepphoton = DeltaR(fakeLepPho->getEta(), fakeLepPho->getPhi(), fakeMu->getEta(), fakeMu->getPhi());
				if(dRlepphoton>0.8){
					if(raw.passMETFilter(METFilter)){

						float fakeLep_deltaPhi = DeltaPhi(fakeMu->getPhi(), METPhi);
						float fakeLep_MT = sqrt(2*MET*fakeMu->getPt()*(1-std::cos(fakeLep_deltaPhi)));
						float fake_ThreeBodyMass = sqrt(2*MET*(fakeLepPho->getP4()+ fakeMu->getP4()).Pt()*(1-std::cos(DeltaR(0, (fakeLepPho->getP4()+fakeMu->getP4()).Phi(), 0, METPhi))));
						fakeLepphoEt = fakeLepPho->getEt();
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

						fakeLepMiniIso = fakeMu->getMiniIso();
						if(fakeMu->isMedium())fakeLepIsStandardProxy = 1;
						else fakeLepIsStandardProxy = 0;
						fakeLepnJet = 0;
						fakeLepHT = 0;
						fakeLep_JetPt.clear();	
						fakeLep_JetEta.clear();
						fakeLep_JetPhi.clear();
						for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
							fakeLep_JetPt.push_back(itJet->getPt());
							fakeLep_JetEta.push_back(itJet->getEta());
							fakeLep_JetPhi.push_back(itJet->getPhi());
							if(!itJet->passSignalSelection())continue;
							if(DeltaR(itJet->getEta(), itJet->getPhi(), fakeLepPho->getEta(), fakeLepPho->getPhi()) <= 0.4)continue;	
							if(DeltaR(itJet->getEta(), itJet->getPhi(), fakeMu->getEta(),  fakeMu->getPhi()) <= 0.4)continue;
							fakeLepnJet += 1;
							fakeLepHT += itJet->getPt();
						}	
          			fakeLep_mcPID.clear();
          		  fakeLep_mcEta.clear();
          		  fakeLep_mcPhi.clear();
          		  fakeLep_mcPt.clear();
          		  fakeLep_mcMomPID.clear();
          		  fakeLep_mcGMomPID.clear();
								if(isMC){
          		 	 for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
          		 	   if(itMC->getEt() < 1.0)continue;
          		 	   fakeLep_mcPID.push_back(itMC->getPID());
          		 	   fakeLep_mcMomPID.push_back(itMC->getMomPID());
          		 	   fakeLep_mcGMomPID.push_back(itMC->getGMomPID());
          		 	   fakeLep_mcEta.push_back(itMC->getEta());
          		 	   fakeLep_mcPhi.push_back(itMC->getPhi());
          		 	   fakeLep_mcPt.push_back(itMC->getEt());
          		 	 }
								}
						fakeLeptree->Fill();
					}//MET Filter
					}//dR filter
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


