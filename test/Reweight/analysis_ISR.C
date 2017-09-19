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

  if(((filter >> 9)&1)!=1){passfilter = false; return passfilter;}
  if(((filter >> 10)&1)!=1){passfilter = false; return passfilter;}

  return passfilter;
}


void analysis_ISR(){//main  

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

  char outputname[100] = "/uscms_data/d3/mengleis/test/resTree_ISR_data.root";
  ofstream logfile;
  logfile.open("/uscms_data/d3/mengleis/test/resTree_ISR_data.log"); 

  logfile << "analysis_mg()" << std::endl;

  RunType datatype(MuonEG2016); 
  TChain* es = new TChain("ggNtuplizer/EventTree");
	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_MuonEG_FebReminiAOD_BCD.root");
	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_MuonEG_FebReminiAOD_EFG.root");
	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_MuonEG_FebReminiAOD_H.root");


//  RunType datatype(SingleMuon2016); 
//  TChain* es = new TChain("ggNtuplizer/EventTree");
//	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_SingleMuon_FebReminiAOD_B.root");
//	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_SingleMuon_FebReminiAOD_C.root");
//	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_SingleMuon_FebReminiAOD_D.root");
//	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_SingleMuon_FebReminiAOD_E.root");
//	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_SingleMuon_FebReminiAOD_F.root");
//	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_SingleMuon_FebReminiAOD_G.root");
//	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_SingleMuon_FebReminiAOD_H.root");

  const unsigned nEvts = es->GetEntries(); 
	float MCweight = 1;	
  logfile << "Total event: " << nEvts << std::endl;
  std::cout << "Total event: " << nEvts << std::endl;
  logfile << "Output file: " << outputname << std::endl;

	int nTotal(0),npassHLT(0), npassPho(0), npassLep(0), npassdR(0), npassZ(0), npassMETFilter(0);

  TFile *outputfile = TFile::Open(outputname,"RECREATE");
  outputfile->cd();

//************ ZGamma Tree **********************//
  TTree *Ztree = new TTree("ZTree","ZTree");
  float ZphoEt(0);
  float ZphoEta(0);
  float ZphoPhi(0);
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
	float ZdilepMass(0);
  int   ZnVertex(0);
  float ZdRPhoLep(0);
  float ZHT(0);
  float ZnJet(0);
	float ZJetPt(0);
 
	Ztree->Branch("MCweight",  &MCweight); 
  Ztree->Branch("phoEt",     &ZphoEt);
  Ztree->Branch("phoEta",    &ZphoEta);
  Ztree->Branch("phoPhi",    &ZphoPhi);
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
	Ztree->Branch("dilepMass", &ZdilepMass);
  Ztree->Branch("nVertex",   &ZnVertex);
  Ztree->Branch("dRPhoLep",  &ZdRPhoLep);
  Ztree->Branch("HT",        &ZHT);
  Ztree->Branch("nJet",      &ZnJet);
	Ztree->Branch("JetPt",     &ZJetPt);
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
  int   ZJetnVertex(0);
  float ZJetdRPhoLep(0);
  float ZJetHT(0);
  float ZJetnJet(0);
	float ZJetJetPt(0);
  
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
  ZJettree->Branch("nVertex",   &ZJetnVertex);
  ZJettree->Branch("dRPhoLep",  &ZJetdRPhoLep);
  ZJettree->Branch("HT",        &ZJetHT);
  ZJettree->Branch("nJet",      &ZJetnJet);
	ZJettree->Branch("JetPt",     &ZJetJetPt);

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

			nTotal+=1;
			if(!raw.passHLT())continue;
			if(raw.nGoodVtx < 1)continue;
			npassHLT+=1;

			int nGoodMu(0);
			for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++){
				if(im->isMedium() && im->getMiniIso() < 0.2 &&  im->getPt() > 15){nGoodMu+=1;}
			}

			if(nGoodMu != 2 || raw.nPho <1)continue;

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
					if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) <= 0.03)GSFveto = false;
					if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.3)FSRVeto=false;
				}
				for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++)
					if(DeltaR(itpho->getEta(), itpho->getPhi(), im->getEta(), im->getPhi()) < 0.3 && im->getEt()>2.0)FSRVeto=false;

				if(itpho->getChIso()<20 && itpho->getSigma()< 0.02 && itpho->isEB()){
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

			bool hasLep(false);
			std::vector<recoMuon>::iterator signalLep = Muon.begin();
			for(std::vector<recoMuon>::iterator itMu = Muon.begin(); itMu != Muon.end(); itMu++){
				if(itMu->getPt() < 25)continue;
				if(!itMu->passHLTSelection())continue;
				if(itMu->passSignalSelection()){
					if(!hasLep){
						hasLep=true; 
						npassLep +=1;
						signalLep = itMu;
					}
				}
			}

			/*********  ZG tree************/ 
			if(hasPho && hasLep){
				double dRlepphoton = DeltaR(signalPho->getEta(), signalPho->getPhi(), signalLep->getEta(), signalLep->getPhi());
				if(dRlepphoton > 0.8){
					if(passMETFilter(METFilter)){ 
						
						bool foundZG(false); 
						std::vector<recoMuon>::iterator trailLep=Muon.begin();	
						for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++){
							if(im == signalLep)continue;
							if(!im->isMedium() || im->getMiniIso() > 0.2)continue;
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
							float ThreeBodyMass = (trailLep->getP4()+signalLep->getP4()+signalPho->getCalibP4()).M(); 

							ZphoEt = signalPho->getCalibEt();
							ZphoEta= signalPho->getEta();
							ZphoPhi= signalPho->getPhi();
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

							ZnJet = 0;
							ZHT = 0;
							TLorentzVector JetVec(0,0,0,0);	
							ZJetPt = 0;
							for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
								if(!itJet->passSignalSelection())continue;
								if(DeltaR(itJet->getEta(), itJet->getPhi(), trailLep->getEta(),trailLep->getPhi()) <= 0.4)continue;	
								if(DeltaR(itJet->getEta(), itJet->getPhi(), signalLep->getEta(),signalLep->getPhi()) <= 0.4)continue;
								ZnJet += 1;
								ZHT += itJet->getPt();
								JetVec = JetVec + itJet->getP4();
							}
							ZJetPt = JetVec.Pt();

							
							Ztree->Fill();
						}	

					}//MET Filter
				}//dR Filter
			}//Candidate Filter

		if(hasLep){ 
			std::vector<recoMuon>::iterator jetMuon = signalLep; 
			for(unsigned ip(0); ip < jetPhoCollection.size(); ip++){
				std::vector<recoPhoton>::iterator jetPho = jetPhoCollection[ip];
				if(!jetPho->isEB())continue;
				double dRlepphoton = DeltaR(jetPho->getEta(), jetPho->getPhi(), jetMuon->getEta(), jetMuon->getPhi());
				if(dRlepphoton>0.8){
				if(passMETFilter(METFilter)){

					bool foundZG(false); 
					std::vector<recoMuon>::iterator trailLep=Muon.begin();	
					for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++){
						if(im == jetMuon)continue;
						if(!im->isMedium() || im->getMiniIso() > 0.2)continue;
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
						TLorentzVector JetVec(0,0,0,0);	
						ZJetJetPt = 0;
						for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
							if(!itJet->passSignalSelection())continue;
							if(DeltaR(itJet->getEta(), itJet->getPhi(), trailLep->getEta(),trailLep->getPhi()) <= 0.4)continue;	
							if(DeltaR(itJet->getEta(), itJet->getPhi(), jetMuon->getEta(),jetMuon->getPhi()) <= 0.4)continue;
							ZJetnJet += 1;
							ZJetHT += itJet->getPt();
							JetVec = JetVec + itJet->getP4();
						}
						ZJetJetPt = JetVec.Pt();

						ZJettree->Fill();
					}

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


