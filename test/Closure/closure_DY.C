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

  return passfilter;
}


void closure_DY(){//main 

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

  char outputname[100] = "/uscms_data/d3/mengleis/closureTree_elefake_DY.root";
  ofstream logfile;
  logfile.open("/uscms_data/d3/mengleis/closureTree_elefake_DY.log"); 

  logfile << "analysis_eg()" << std::endl;
  logfile << "medium eleID+miniIso" << std::endl;
  //logfile << "Loose the proxy definition: no upper bounds for photon; LooseFakeProxy for electron" << std::endl;

  RunType datatype(MC); 
  TChain* es = new TChain("ggNtuplizer/EventTree");
  TFileCollection fc("dum","","DY.txt");
  es->AddFileInfoList((TCollection*)fc.GetList());

  const unsigned nEvts = es->GetEntries(); 
  logfile << "Total event: " << nEvts << std::endl;
  std::cout << "Total event: " << nEvts << std::endl;
  logfile << "Output file: " << outputname << std::endl;

  int nTotal(0),npassHLT(0), npassPho(0), npassLep(0), npassdR(0), npassZ(0), npassMETFilter(0);

  TFile *outputfile = TFile::Open(outputname,"RECREATE");
  outputfile->cd();

	float crosssection = 5670;
	float ntotalevent  = nEvts;
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
  float nJet(0);
  std::vector<int> mcPID;
  std::vector<float> mcEta;
  std::vector<float> mcPhi;
  std::vector<float> mcPt;
  std::vector<int> mcMomPID;
  std::vector<int> mcGMomPID;
 
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
  sigtree->Branch("mcPID",     &mcPID);
  sigtree->Branch("mcEta",     &mcEta);
  sigtree->Branch("mcPhi",     &mcPhi);
  sigtree->Branch("mcPt",      &mcPt);
  sigtree->Branch("mcMomPID",  &mcMomPID);
  sigtree->Branch("mcGMomPID", &mcGMomPID);
	sigtree->Branch("crosssection",&crosssection);
	sigtree->Branch("ntotalevent", &ntotalevent);

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
  float proxynJet(0);
  std::vector<int>   proxymcPID;
  std::vector<float> proxymcEta;
  std::vector<float> proxymcPhi;
  std::vector<float> proxymcPt;
  std::vector<int>   proxymcMomPID;
  std::vector<int>   proxymcGMomPID;
  
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
  proxytree->Branch("mcPID",     &proxymcPID);
  proxytree->Branch("mcEta",     &proxymcEta);
  proxytree->Branch("mcPhi",     &proxymcPhi);
  proxytree->Branch("mcPt",      &proxymcPt);
  proxytree->Branch("mcMomPID",  &proxymcMomPID);
  proxytree->Branch("mcGMomPID", &proxymcGMomPID);
	proxytree->Branch("crosssection",&crosssection);
	proxytree->Branch("ntotalevent", &ntotalevent);

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
  jettree->Branch("nVertex",   &jetnVertex);
  jettree->Branch("dRPhoLep",  &jetdRPhoLep);
  jettree->Branch("HT",        &jetHT);
  jettree->Branch("nJet",      &jetnJet);
//*********** fake lepton *********************//
  TTree *fakeEtree = new TTree("fakeETree","fakeETree");
  float fakeEphoEt(0);
  float fakeEphoEta(0);
  float fakeEphoPhi(0);
  float fakeElepPt(0);
  float fakeElepEta(0);
  float fakeElepPhi(0);
  float fakeEsigMT(0);
  float fakeEsigMET(0);
  float fakeEsigMETPhi(0);
  float fakeEdPhiLepMET(0);
  int   fakeEnVertex(0);
  float fakeEdRPhoLep(0);
  float fakeEHT(0);
  float fakeEnJet(0);
  
  fakeEtree->Branch("phoEt",     &fakeEphoEt);
  fakeEtree->Branch("phoEta",    &fakeEphoEta);
  fakeEtree->Branch("phoPhi",    &fakeEphoPhi);
  fakeEtree->Branch("lepPt",     &fakeElepPt);
  fakeEtree->Branch("lepEta",    &fakeElepEta);
  fakeEtree->Branch("lepPhi",    &fakeElepPhi);
  fakeEtree->Branch("sigMT",     &fakeEsigMT);
  fakeEtree->Branch("sigMET",    &fakeEsigMET);
  fakeEtree->Branch("sigMETPhi", &fakeEsigMETPhi);
  fakeEtree->Branch("dPhiLepMET",&fakeEdPhiLepMET);
  fakeEtree->Branch("nVertex",   &fakeEnVertex);
  fakeEtree->Branch("dRPhoLep",  &fakeEdRPhoLep);
  fakeEtree->Branch("HT",        &fakeEHT);
  fakeEtree->Branch("nJet",      &fakeEnJet);

//*********** histo list **********************//
  TH1F *p_photonEt = new TH1F("photonET","#gamma E_{T}; E_{T} (GeV)",500,0,1500);
  TH1F *p_photonEta = new TH1F("photonEta","#gamma #eta; #eta;",60,-3,3);
  TH1F *p_lepPt = new TH1F("eleET","electron P_{T}; P_{T} (GeV)",500,0,1000);
  TH1F *p_lepEta = new TH1F("lepEta","electron #eta; #eta;",60,-3,3);
  TH1F *p_Mt = new TH1F("Mt","M_{T}; M_{T} (GeV);",200,0,400); 
  TH1F *p_MET = new TH1F("MET","MET",100,0,100);
 
  TH1F *p_PhoELeDeltaR = new TH1F("p_PhoELeDeltaR","#DeltaR(e, #gamma); #DeltaR(e, #gamma);",100,0,10);
  TH1F *p_PhoEleMass = new TH1F("p_PhoEleMass","M_{e#gamma}; M_{e#gamma}(GeV);", 200,0,400);
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
			if(((raw.HLTPho >> 14) &1) ==0)continue;
			if(raw.nGoodVtx < 1)continue;
			npassHLT+=1;

			if(raw.nEle < 1 || raw.nPho <1)continue;

			bool hasPho(false);
			std::vector<recoPhoton>::iterator signalPho = Photon.begin();
			std::vector< std::vector<recoPhoton>::iterator >  proxyPhoCollection;
			proxyPhoCollection.clear();
			std::vector< std::vector<recoPhoton>::iterator >  jetPhoCollection;
			jetPhoCollection.clear();
			for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
				if(itpho->getR9() < 0.5)continue;
				if(!itpho->passBasicSelection())continue;
				bool passSigma = itpho->passSigma(1);
				bool passChIso = itpho->passChIso(1);
				bool PixelVeto = itpho->PixelSeed()==0? true: false;
				bool GSFveto(true);
				bool photonFSRVeto(true);
				bool eleFSRVeto(true);
				for(std::vector<recoEle>::iterator ie = Ele.begin(); ie != Ele.end(); ie++){
					if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) <= 0.03)GSFveto = false;
					if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.3 && ie->getEt()>2.0)photonFSRVeto=false;
					if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.3 && DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) > 0.05 && ie->getEt()>2.0)eleFSRVeto=false;
				}
				for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++){
					if(DeltaR(itpho->getEta(), itpho->getPhi(), im->getEta(), im->getPhi()) < 0.3 && im->getEt()>2.0){
						photonFSRVeto=false;
						eleFSRVeto=false;
					}
				}

				if(itpho->getChIso()<20 && itpho->getSigma()< 0.02 && itpho->isEB()){
				//if(itpho->isEB()){  // Loose the proxy definition
					if(!passSigma || !passChIso){
						if(GSFveto && PixelVeto && photonFSRVeto)jetPhoCollection.push_back(itpho);
					}
				}

				if(!itpho->passSignalSelection())continue;
				if(!itpho->fireDoubleTrg(5) && !itpho->fireDoubleTrg(6))continue;
				if(GSFveto && PixelVeto && photonFSRVeto){
					if(!hasPho){
						hasPho=true;
						npassPho +=1;
						signalPho = itpho;
					}
				}

				//if(eleFSRVeto && (!PixelVeto || !GSFveto)){
				if((!PixelVeto || !GSFveto)){
					if(itpho->isEB())proxyPhoCollection.push_back(itpho);
				}
			}

			bool hasEle(false);
			std::vector<recoEle>::iterator signalEle = Ele.begin();
			std::vector< std::vector<recoEle>::iterator > proxyEleCollection;
			proxyEleCollection.clear();
			std::vector< std::vector<recoEle>::iterator > fakeEleCollection;
			fakeEleCollection.clear();
			for(std::vector<recoEle>::iterator itEle = Ele.begin(); itEle != Ele.end(); itEle++){
				if((itEle->isEB() && itEle->getR9() < 0.5) || (itEle->isEE() && itEle->getR9() < 0.8))continue;
				if(!itEle->fireTrgs(21) && !itEle->fireTrgs(22))continue;
				if(itEle->getCalibPt() < 25)continue;
				if(itEle->isFakeProxy())fakeEleCollection.push_back(itEle);	
				//if(itEle->isLooseFakeProxy())fakeEleCollection.push_back(itEle);//Loose the proxy definition	
				if(itEle->passSignalSelection()){
					proxyEleCollection.push_back(itEle);
					if(!hasEle && hasPho){
						hasEle=true; 
						npassLep +=1;
						signalEle = itEle;
					}
				}
			}

			bool saveMC(false);
			 
			if(hasPho && hasEle && signalPho->isEB()){
				double dReg = DeltaR(signalPho->getEta(), signalPho->getPhi(), signalEle->getEta(), signalEle->getPhi()); 
				p_PhoELeDeltaR->Fill(dReg);
				if(dReg > 0.8){
					npassdR+=1;
					p_PhoEleMass->Fill((signalPho->getCalibP4()+signalEle->getCalibP4()).M());
					if(fabs((signalPho->getCalibP4()+signalEle->getCalibP4()).M() - 91.188) > 10.0){

						npassZ+=1;
						//if(passMETFilter(METFilter)){
							npassMETFilter +=1;

							float deltaPhi = DeltaPhi(signalEle->getPhi(), METPhi);
							float MT = sqrt(2*MET*signalEle->getCalibPt()*(1-std::cos(deltaPhi)));

							phoEt = signalPho->getCalibEt();
							phoEta= signalPho->getEta();
							phoPhi= signalPho->getPhi();
							lepPt = signalEle->getCalibPt();
							lepEta= signalEle->getEta();
							lepPhi= signalEle->getPhi();
							sigMT = MT;
							sigMET= MET;
							sigMETPhi = METPhi;
							dPhiLepMET = deltaPhi; 
							nVertex = nVtx; 
							dRPhoLep= dReg;

							nJet = 0;
							HT = 0;
							for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
								if(!itJet->passSignalSelection())continue;
								if(DeltaR(itJet->getEta(), itJet->getPhi(), signalPho->getEta(),signalPho->getPhi()) <= 0.4)continue;	
								if(DeltaR(itJet->getEta(), itJet->getPhi(), signalEle->getEta(),signalEle->getPhi()) <= 0.4)continue;
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
								float mcdRe = DeltaR(signalEle->getEta(),signalEle->getPhi(), itMC->getEta(), itMC->getPhi());
								if(mcdR < 0.3 || mcdRe < 0.3){
									mcPID.push_back(itMC->getPID());
									mcMomPID.push_back(itMC->getMomPID());
									mcEta.push_back(itMC->getEta());      
									mcPhi.push_back(itMC->getPhi());
									mcPt.push_back(itMC->getEt());
									mcGMomPID.push_back(itMC->getGMomPID());
								}
							}

							sigtree->Fill();

							if(datatype == MC)saveMC=true;
//						}//MET Filter
					}// Z mass Filter
				}//dR filter
			}// ele + pho candidate
	 

			for(unsigned ip(0); ip < proxyPhoCollection.size(); ip++){
				for(unsigned ie(0); ie < proxyEleCollection.size(); ie++){
					std::vector<recoPhoton>::iterator proxyPho = proxyPhoCollection[ip];
					std::vector<recoEle>::iterator proxyEle = proxyEleCollection[ie];
					double dReg = DeltaR(proxyPho->getEta(), proxyPho->getPhi(), proxyEle->getEta(), proxyEle->getPhi());
					if(dReg>0.8){
						if(fabs((proxyPho->getCalibP4()+proxyEle->getCalibP4()).M() - 91.188) > 10.0){
							if(passMETFilter(METFilter)){

								float proxy_deltaPhi = DeltaPhi(proxyEle->getPhi(), METPhi);
								float proxy_MT = sqrt(2*MET*proxyEle->getCalibPt()*(1-std::cos(proxy_deltaPhi)));
								proxyphoEt = proxyPho->getCalibEt();
								proxyphoEta= proxyPho->getEta();
								proxyphoPhi= proxyPho->getPhi();
								proxylepPt = proxyEle->getCalibPt();
								proxylepEta= proxyEle->getEta();
								proxylepPhi= proxyEle->getPhi();
								proxysigMT = proxy_MT;
								proxysigMET= MET;
								proxysigMETPhi = METPhi;
								proxydPhiLepMET = proxy_deltaPhi; 
								proxynVertex = nVtx; 
								proxydRPhoLep= dReg;
								proxytree->Fill();

								proxynJet = 0;
								proxyHT = 0;
								for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
									if(!itJet->passSignalSelection())continue;
									if(DeltaR(itJet->getEta(), itJet->getPhi(), proxyPho->getEta(),proxyPho->getPhi()) <= 0.4)continue;	
									if(DeltaR(itJet->getEta(), itJet->getPhi(), proxyEle->getEta(),proxyEle->getPhi()) <= 0.4)continue;
									proxynJet += 1;
									proxyHT += itJet->getPt();
								}

								proxymcPID.clear();
								proxymcEta.clear();
								proxymcPhi.clear();
								proxymcPt.clear();
								proxymcMomPID.clear();
								proxymcGMomPID.clear();
								for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
									if(itMC->getEt() < 1.0)continue;
									float mcdR = DeltaR(proxyPho->getEta(), proxyPho->getPhi(), itMC->getEta(), itMC->getPhi());
									float mcdRe = DeltaR(proxyEle->getEta(),proxyEle->getPhi(), itMC->getEta(), itMC->getPhi());
									if(mcdR < 0.3 || mcdRe < 0.3){
										proxymcPID.push_back(itMC->getPID());
										proxymcMomPID.push_back(itMC->getMomPID());
										proxymcEta.push_back(itMC->getEta());      
										proxymcPhi.push_back(itMC->getPhi());
										proxymcPt.push_back(itMC->getEt());
										proxymcGMomPID.push_back(itMC->getGMomPID());
									}
								}
	
								if(datatype == MC)saveMC=true;
							}//MET Filter
						}// Z mass Filter
					}//dR filter
				}// loop on ele collection
			} // loop on pho collection

			
			if(proxyEleCollection.size() > 0){
				std::vector<recoEle>::iterator jetEle = proxyEleCollection[0];
				for(unsigned ip(0); ip < jetPhoCollection.size(); ip++){
					std::vector<recoPhoton>::iterator jetPho = jetPhoCollection[ip];
					double dReg = DeltaR(jetPho->getEta(), jetPho->getPhi(), jetEle->getEta(), jetEle->getPhi());
					if(dReg>0.8){
						if(fabs((jetPho->getCalibP4()+jetEle->getCalibP4()).M() - 91.188) > 10.0){
							if(passMETFilter(METFilter)){


								float jet_deltaPhi = DeltaPhi(jetEle->getPhi(), METPhi);
								float jet_MT = sqrt(2*MET*jetEle->getCalibPt()*(1-std::cos(jet_deltaPhi)));
								jetphoEt = jetPho->getCalibEt();
								jetphoEta= jetPho->getEta();
								jetphoPhi= jetPho->getPhi();
								jetlepPt = jetEle->getCalibPt();
								jetlepEta= jetEle->getEta();
								jetlepPhi= jetEle->getPhi();
								jetsigMT = jet_MT;
								jetsigMET= MET;
								jetsigMETPhi = METPhi;
								jetdPhiLepMET = jet_deltaPhi; 
								jetnVertex = nVtx; 
								jetdRPhoLep= dReg;
								jettree->Fill();

								jetnJet = 0;
								jetHT = 0;
								for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
									if(!itJet->passSignalSelection())continue;
									if(DeltaR(itJet->getEta(), itJet->getPhi(), jetPho->getEta(),jetPho->getPhi()) <= 0.4)continue;	
									if(DeltaR(itJet->getEta(), itJet->getPhi(), jetEle->getEta(),jetEle->getPhi()) <= 0.4)continue;
									jetnJet += 1;
									jetHT += itJet->getPt();
								}

								if(datatype == MC)saveMC=true;
							}//MET Filter
						}// Z mass Filter
					}//dR filter
				} // loop on pho collection
			}

			if(hasPho && !hasEle && signalPho->isEB()){
				std::vector<recoPhoton>::iterator fakeEPho = signalPho;
				for(unsigned ip(0); ip < fakeEleCollection.size(); ip++){
					std::vector<recoEle>::iterator fakeEle = fakeEleCollection[ip];
					double dReg = DeltaR(fakeEPho->getEta(), fakeEPho->getPhi(), fakeEle->getEta(), fakeEle->getPhi());
					if(dReg>0.8){
						if(fabs((fakeEPho->getCalibP4()+fakeEle->getCalibP4()).M() - 91.188) > 10.0){
							if(passMETFilter(METFilter)){

								float fakeE_deltaPhi = DeltaPhi(fakeEle->getPhi(), METPhi);
								float fakeE_MT = sqrt(2*MET*fakeEle->getCalibPt()*(1-std::cos(fakeE_deltaPhi)));
								fakeEphoEt = fakeEPho->getCalibEt();
								fakeEphoEta= fakeEPho->getEta();
								fakeEphoPhi= fakeEPho->getPhi();
								fakeElepPt = fakeEle->getCalibPt();
								fakeElepEta= fakeEle->getEta();
								fakeElepPhi= fakeEle->getPhi();
								fakeEsigMT = fakeE_MT;
								fakeEsigMET= MET;
								fakeEsigMETPhi = METPhi;
								fakeEdPhiLepMET = fakeE_deltaPhi; 
								fakeEnVertex = nVtx; 
								fakeEdRPhoLep= dReg;
								fakeEtree->Fill();
									
								fakeEnJet = 0;
								fakeEHT = 0;
								for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
									if(!itJet->passSignalSelection())continue;
									if(DeltaR(itJet->getEta(), itJet->getPhi(), fakeEPho->getEta(), fakeEPho->getPhi()) <= 0.4)continue;	
									if(DeltaR(itJet->getEta(), itJet->getPhi(), fakeEle->getEta(),  fakeEle->getPhi()) <= 0.4)continue;
									fakeEnJet += 1;
									fakeEHT += itJet->getPt();
								}	

								if(datatype == MC)saveMC=true;
							}//MET Filter
						}// Z mass Filter
					}//dR filter
				} // loop on pho collection
			}
	}//loop on  events

  p_eventcount->GetXaxis()->SetBinLabel(1,"nTotal");
  p_eventcount->GetXaxis()->SetBinLabel(2,"npassHLT");
  p_eventcount->GetXaxis()->SetBinLabel(3,"npassPho");
  p_eventcount->GetXaxis()->SetBinLabel(4,"npassLep");
  p_eventcount->GetXaxis()->SetBinLabel(5,"npassdR");
  p_eventcount->GetXaxis()->SetBinLabel(6,"npassZ");
  p_eventcount->GetXaxis()->SetBinLabel(7,"npassMETFilter");
  p_eventcount->Fill(0.5, nTotal);
  p_eventcount->Fill(1.5, npassHLT);
  p_eventcount->Fill(2.5, npassPho);
  p_eventcount->Fill(3.5, npassLep);
  p_eventcount->Fill(4.5, npassdR);
  p_eventcount->Fill(5.5, npassZ);
  p_eventcount->Fill(6.5, npassMETFilter);

	outputfile->Write();
	outputfile->Close();
	logfile.close();
}


