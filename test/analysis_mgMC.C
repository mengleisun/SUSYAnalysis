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

//  if(((filter >> 9)&1)!=1){passfilter = false; return passfilter;}
//  if(((filter >> 10)&1)!=1){passfilter = false; return passfilter;}

  return passfilter;
}


void analysis_mgMC(){//main  

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

  char outputname[100] = "/uscms_data/d3/mengleis/test/resTree_mgsignal_DY.root";
  ofstream logfile;
  logfile.open("/uscms_data/d3/mengleis/test/resTree_mgsignal_DY.log"); 

  logfile << "analysis_mg()" << std::endl;

  RunType datatype(MC); 
  TChain* es = new TChain("ggNtuplizer/EventTree");
	es->Add("root://cmseos.fnal.gov//store/user/msun/MCSummer16/DYJetsToLL_M-50_nlo.root");

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
  float nJet(0);
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
  sigtree->Branch("mcPID",     &mcPID);
  sigtree->Branch("mcEta",     &mcEta);
  sigtree->Branch("mcPhi",     &mcPhi);
  sigtree->Branch("mcPt",      &mcPt);
  sigtree->Branch("mcMomPID",  &mcMomPID);
  sigtree->Branch("mcGMomPID", &mcGMomPID);

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

  TTree *fakeratetree = new TTree("fakerateTree","fakerateTree");
  float fakeratephoEt(0);
  float fakeratephoEta(0);
  float fakeratephoPhi(0);
	float fakeratephoSigma(0);
	float fakeratephoChIso(0);
  float fakeratelepPt(0);
  float fakeratelepEta(0);
  float fakeratelepPhi(0);
  float fakeratesigMT(0);
  float fakeratesigMET(0);
  float fakeratesigMETPhi(0);
  float fakeratedPhiLepMET(0);
	float fakeratethreeMass(0);
  int   fakeratenVertex(0);
  float fakeratedRPhoLep(0);

  fakeratetree->Branch("phoEt",     &fakeratephoEt);
  fakeratetree->Branch("phoEta",    &fakeratephoEta);
  fakeratetree->Branch("phoPhi",    &fakeratephoPhi);
	fakeratetree->Branch("phoSigma",  &fakeratephoSigma);
	fakeratetree->Branch("phoChIso",  &fakeratephoChIso);
  fakeratetree->Branch("lepPt",     &fakeratelepPt);
  fakeratetree->Branch("lepEta",    &fakeratelepEta);
  fakeratetree->Branch("lepPhi",    &fakeratelepPhi);
  fakeratetree->Branch("sigMT",     &fakeratesigMT);
  fakeratetree->Branch("sigMET",    &fakeratesigMET);
  fakeratetree->Branch("sigMETPhi", &fakeratesigMETPhi);
  fakeratetree->Branch("dPhiLepMET",&fakeratedPhiLepMET);
	fakeratetree->Branch("threeMass", &fakeratethreeMass);
  fakeratetree->Branch("nVertex",   &fakeratenVertex);
  fakeratetree->Branch("dRPhoLep",  &fakeratedRPhoLep);
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
			if(raw.nGoodVtx < 1)continue;
			npassHLT+=1;

			if(raw.nMu < 1 || raw.nPho <1)continue;

			bool hasPho(false);
			std::vector<recoPhoton>::iterator signalPho = Photon.begin();
			std::vector< std::vector<recoPhoton>::iterator >  jetPhoCollection;
			jetPhoCollection.clear();
			bool hasFakeRatePho(false);
			std::vector<recoPhoton>::iterator fakeratePho = Photon.begin();
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
						if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) <= 0.02)GSFveto = false;
						if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.3)FSRVeto=false;
					}
					for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++)
						if(DeltaR(itpho->getEta(), itpho->getPhi(), im->getEta(), im->getPhi()) < 0.3 && im->getEt()>2.0)FSRVeto=false;

					if(itpho->getChIso()<20 && itpho->getSigma()< 0.02 && itpho->isEB()){
						if(!hasFakeRatePho && GSFveto && PixelVeto && FSRVeto){
							hasFakeRatePho = true;
							fakeratePho = itpho;
						}
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

				}
			}

			bool hasLep(false);
			std::vector<recoMuon>::iterator signalLep = Muon.begin();
			std::vector< std::vector<recoMuon>::iterator > proxyLepCollection;
			proxyLepCollection.clear();
			for(std::vector<recoMuon>::iterator itMu = Muon.begin(); itMu != Muon.end(); itMu++){
				if(itMu->getPt() < 25)continue;
				if(itMu->fireSingleTrg(2) || itMu->fireSingleTrg(21) || itMu->fireSingleTrg(22)){  //Mu+gamma HLT
					if(itMu->passSignalSelection()){
						if(proxyLepCollection.size() == 0)proxyLepCollection.push_back(itMu);
						if(!hasLep){
							hasLep=true; 
							npassLep +=1;
							signalLep = itMu;
						}
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

						sigtree->Fill();

					}//MET Filter
				}//dR Filter
			}//Candidate Filter
		 

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

		if(hasFakeRatePho && hasLep){
			double dRlepphoton = DeltaR(fakeratePho->getEta(), fakeratePho->getPhi(), signalLep->getEta(), signalLep->getPhi());
			if(dRlepphoton > 0.8){
				if(passMETFilter(METFilter)){ 

					float deltaPhi = DeltaPhi(signalLep->getPhi(), METPhi);
					float MT = sqrt(2*MET*signalLep->getPt()*(1-std::cos(deltaPhi)));
					float ThreeBodyMass = sqrt(2*MET*(fakeratePho->getP4()+ signalLep->getP4()).Pt()*(1-std::cos(DeltaR(0, (fakeratePho->getP4()+signalLep->getP4()).Phi(), 0, METPhi))));

					fakeratephoEt = fakeratePho->getCalibEt();
					fakeratephoEta= fakeratePho->getEta();
					fakeratephoPhi= fakeratePho->getPhi();
					fakeratephoSigma = fakeratePho->getSigma();
					fakeratephoChIso = fakeratePho->getChIso();
					fakeratelepPt = signalLep->getPt();
					fakeratelepEta= signalLep->getEta();
					fakeratelepPhi= signalLep->getPhi();
					fakeratesigMT = MT;
					fakeratesigMET= MET;
					fakeratesigMETPhi = METPhi;
					fakeratedPhiLepMET = deltaPhi;
					fakeratethreeMass = ThreeBodyMass;
					fakeratenVertex = nVtx;
					fakeratedRPhoLep= dRlepphoton;

					fakeratetree->Fill();

				}//MET Filter
			}//dR Filter
		}//Candidate Filter

 
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


