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

#include "../include/analysis_rawData.h"
#include "../include/analysis_photon.h"
#include "../include/analysis_muon.h"
#include "../include/analysis_ele.h"
#include "../include/analysis_mcData.h"
#include "../include/analysis_tools.h"
#include "../include/analysis_jet.h"

bool passMETFilter(int filter){
  bool passfilter(true);
  for(int im(1); im <= 8; im++)
    if(((filter >> im)&1)!=0){passfilter = false; return passfilter;}

  if(((filter >> 9)&1)!=1){passfilter = false; return passfilter;}
  if(((filter >> 10)&1)!=1){passfilter = false; return passfilter;}

  return passfilter;
}


void analysis_mg(){//main  

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

  char outputname[100] = "/uscms_data/d3/mengleis/Sep1/resTree_mgsignal_MuonEG_MiniIso.root";
  ofstream logfile;
  logfile.open("/uscms_data/d3/mengleis/Sep1/resTree_mgsignal_MuonEG_MiniIso.log"); 

  logfile << "analysis_mg()" << std::endl;
  logfile << "miniIso; one lepton for fakephoton background" << std::endl;
	
  RunType datatype(MuonEG2016); 
  TChain* es = new TChain("ggNtuplizer/EventTree");
	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_MuonEG_FebReminiAOD_BCD.root");
	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_MuonEG_FebReminiAOD_EFG.root");
	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_MuonEG_FebReminiAOD_H.root");

  const unsigned nEvts = es->GetEntries(); 
  logfile << "Total event: " << nEvts << std::endl;
  std::cout << "Total event: " << nEvts << std::endl;
  logfile << "Output file: " << outputname << std::endl;

	int nTotal(0),npassHLT(0), npassPho(0), npassLep(0), npassdR(0), npassZ(0), npassMETFilter(0);

  TFile *outputfile = TFile::Open(outputname,"RECREATE");
  outputfile->cd();

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
	int   nJetClean(0);
	float HTClean(0);

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
	sigtree->Branch("trailPt",   &trailPt);
	sigtree->Branch("trailEta",  &trailEta);
	sigtree->Branch("trailPhi",  &trailPhi);
	sigtree->Branch("nJetClean", &nJetClean);
	sigtree->Branch("HTClean",   &HTClean);

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
  int   proxynJet(0);
	float proxytrailPt(0);
	float proxytrailEta(0);
	float proxytrailPhi(0);
	int   proxynJetClean(0);
	float proxyHTClean(0);
  
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
	proxytree->Branch("trailPt",   &proxytrailPt);
	proxytree->Branch("trailEta",  &proxytrailEta);
	proxytree->Branch("trailPhi",  &proxytrailPhi);
	proxytree->Branch("nJetClean", &proxynJetClean);
	proxytree->Branch("HTClean",   &proxyHTClean);

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
  int   jetnJet(0);
	float jettrailPt(0);
	float jettrailEta(0);
	float jettrailPhi(0);
	int   jetnJetClean(0);
	float jetHTClean(0);
  
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
	jettree->Branch("trailPt",   &jettrailPt);
	jettree->Branch("trailEta",  &jettrailEta);
	jettree->Branch("trailPhi",  &jettrailPhi);
	jettree->Branch("nJetClean", &jetnJetClean);
	jettree->Branch("HTClean",   &jetHTClean);
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
  int   fakeLepnJet(0);
	int   fakeLepnJetClean(0);
	float fakeLepHTClean(0);
  
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
	fakeLeptree->Branch("nJetClean", &fakeLepnJetClean);
	fakeLeptree->Branch("HTClean",   &fakeLepHTClean);

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
			if(datatype == MC)for(int iMC(0); iMC < raw.nMC; iMC++){MCData.push_back(mcData(raw, iMC));}
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

			bool hasPho(false);
			std::vector<recoPhoton>::iterator signalPho = Photon.begin();
			std::vector< std::vector<recoPhoton>::iterator >  proxyPhoCollection;
			proxyPhoCollection.clear();
			std::vector< std::vector<recoPhoton>::iterator >  jetPhoCollection;
			jetPhoCollection.clear();
			for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
				if(itpho->getR9() < 0.5)continue;
				if(!itpho->passHLTSelection())continue;
				if(!itpho->passBasicSelection())continue;
				bool passSigma = itpho->passSigma(1);
				bool passChIso = itpho->passChIso(1);
				bool PixelVeto = itpho->PixelSeed()==0? true: false;
				bool GSFveto(true);
				bool FSRVeto(true);
				for(std::vector<recoEle>::iterator ie = Ele.begin(); ie != Ele.end(); ie++){
					if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) <= 0.02)GSFveto = false;
					if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.3)FSRVeto=false;
				}
				for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++)
					if(DeltaR(itpho->getEta(), itpho->getPhi(), im->getEta(), im->getPhi()) < 0.3 && im->getEt()>2.0)FSRVeto=false;

				if(itpho->getChIso()<20 && itpho->getSigma()< 0.02 && itpho->isEB()){
					//if(itpho->isEB())  // Loose the proxy definition
					if(!passSigma || !passChIso){
						if(GSFveto && PixelVeto && FSRVeto)jetPhoCollection.push_back(itpho);
					}
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

			if(hasPho && hasLep && signalPho->isEB()){
				double dRlepphoton = DeltaR(signalPho->getEta(), signalPho->getPhi(), signalLep->getEta(), signalLep->getPhi());
				if(dRlepphoton > 0.8){
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
						nJetClean = 0;
						HTClean = 0;
						for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
							if(!itJet->passSignalSelection())continue;
							if(DeltaR(itJet->getEta(), itJet->getPhi(), signalPho->getEta(),signalPho->getPhi()) <= 0.4)continue;	
							if(DeltaR(itJet->getEta(), itJet->getPhi(), signalLep->getEta(),signalLep->getPhi()) <= 0.4)continue;
							nJet += 1;
							HT += itJet->getPt();
							bool isclean(true);
							for(unsigned iMu(0); iMu < miniisoLep.size(); iMu++){
								std::vector<recoMuon>::iterator itMu = miniisoLep[iMu];
								if(DeltaR(itJet->getEta(), itJet->getPhi(), itMu->getEta(), itMu->getPhi()) <= 0.4)isclean=false;
							}
							if(isclean){
								nJetClean += 1;
								HTClean += itJet->getPt();
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
				proxynJetClean = 0;
				proxyHTClean = 0;
				for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
					if(!itJet->passSignalSelection())continue;
					if(DeltaR(itJet->getEta(), itJet->getPhi(), proxyPho->getEta(),proxyPho->getPhi()) <= 0.4)continue;	
					if(DeltaR(itJet->getEta(), itJet->getPhi(), proxyMuon->getEta(),proxyMuon->getPhi()) <= 0.4)continue;
					proxynJet += 1;
					proxyHT += itJet->getPt();
					bool isclean(true);
					for(unsigned iMu(0); iMu < miniisoLep.size(); iMu++){
						std::vector<recoMuon>::iterator itMu = miniisoLep[iMu];
						if(DeltaR(itJet->getEta(), itJet->getPhi(), itMu->getEta(), itMu->getPhi()) <= 0.4)isclean=false;
					}
					if(isclean){
						proxynJetClean += 1;
						proxyHTClean += itJet->getPt();
					}
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
				jetnJetClean = 0;
				jetHTClean = 0;
				for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
					if(!itJet->passSignalSelection())continue;
					//if(DeltaR(itJet->getEta(), itJet->getPhi(), jetPho->getEta(),jetPho->getPhi()) <= 0.4)continue;	
					if(DeltaR(itJet->getEta(), itJet->getPhi(), jetMuon->getEta(),jetMuon->getPhi()) <= 0.4)continue;
					jetnJet += 1;
					jetHT += itJet->getPt();
					bool isclean(true);
					for(unsigned iMu(0); iMu < miniisoLep.size(); iMu++){
						std::vector<recoMuon>::iterator itMu = miniisoLep[iMu];
						if(DeltaR(itJet->getEta(), itJet->getPhi(), itMu->getEta(), itMu->getPhi()) <= 0.4)isclean=false;
					}
					if(isclean){
						jetnJetClean += 1;
						jetHTClean += itJet->getPt();
					}
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
				if(dRlepphoton>0.8){
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

						fakeLepMiniIso = fakeMu->getMiniIso();
						if(fakeMu->isMedium())fakeLepIsMedium = 1;
						else fakeLepIsMedium = 0; 
						//fakeLepMiniIso = fakeMu->getRelIso();
						fakeLepnJet = 0;
						fakeLepHT = 0;
						for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
							if(!itJet->passSignalSelection())continue;
							if(DeltaR(itJet->getEta(), itJet->getPhi(), fakeLepPho->getEta(), fakeLepPho->getPhi()) <= 0.4)continue;	
							//if(DeltaR(itJet->getEta(), itJet->getPhi(), fakeMu->getEta(),  fakeMu->getPhi()) <= 0.4)continue;
							fakeLepnJet += 1;
							fakeLepHT += itJet->getPt();
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


