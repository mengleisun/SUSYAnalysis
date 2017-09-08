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


void analysis_ZISR(){//main  

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

  char outputname[100] = "/uscms_data/d3/mengleis/test/resTree_ZISR_data.root";
  ofstream logfile;
  logfile.open("/uscms_data/d3/mengleis/test/resTree_ZISR_data.log"); 

  logfile << "analysis_mg()" << std::endl;

//  RunType datatype(MuonEG2016); 
//  TChain* es = new TChain("ggNtuplizer/EventTree");
//	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_MuonEG_FebReminiAOD_BCD.root");
//	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_MuonEG_FebReminiAOD_EFG.root");
//	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_MuonEG_FebReminiAOD_H.root");

  RunType datatype(SingleMuon2016); 
  TChain* es = new TChain("ggNtuplizer/EventTree");
	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_SingleMuon_FebReminiAOD_B.root");
	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_SingleMuon_FebReminiAOD_C.root");
	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_SingleMuon_FebReminiAOD_D.root");
	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_SingleMuon_FebReminiAOD_E.root");
	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_SingleMuon_FebReminiAOD_F.root");
	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_SingleMuon_FebReminiAOD_G.root");
	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_SingleMuon_FebReminiAOD_H.root");

  const unsigned nEvts = es->GetEntries(); 
	float MCweight = 1;	
  logfile << "Total event: " << nEvts << std::endl;
  std::cout << "Total event: " << nEvts << std::endl;
  logfile << "Output file: " << outputname << std::endl;

	int nTotal(0),npassHLT(0), npassPho(0), npassLep(0), npassdR(0), npassZ(0), npassMETFilter(0);

  TFile *outputfile = TFile::Open(outputname,"RECREATE");
  outputfile->cd();

	TH1F *p_nMu = new TH1F("p_nMu","p_nMu",10,0,10);

//************ ZGamma Tree **********************//
  TTree *Ztree = new TTree("ZTree","ZTree");
  float ZlepPt(0);
  float ZlepEta(0);
  float ZlepPhi(0);
	float ZtrailPt(0);
	float ZtrailEta(0);
  float ZsigMT(0);
  float ZsigMET(0);
  float ZsigMETPhi(0);
  float ZdPhiLepMET(0);
	float ZdilepMass(0);
  int   ZnVertex(0);
  float ZHT(0);
  float ZnJet(0);
	int   ZnMu(0);
	int   ZnGoodMu(0);
	std::vector<float> ZJetPt;
	std::vector<int>   ZJetMatchPho;
 
	Ztree->Branch("MCweight",  &MCweight); 
  Ztree->Branch("lepPt",     &ZlepPt);
  Ztree->Branch("lepEta",    &ZlepEta);
  Ztree->Branch("lepPhi",    &ZlepPhi);
	Ztree->Branch("trailPt",   &ZtrailPt);
	Ztree->Branch("trailEta",  &ZtrailEta);
  Ztree->Branch("sigMT",     &ZsigMT);
  Ztree->Branch("sigMET",    &ZsigMET);
  Ztree->Branch("sigMETPhi", &ZsigMETPhi);
  Ztree->Branch("dPhiLepMET",&ZdPhiLepMET);
	Ztree->Branch("dilepMass", &ZdilepMass);
  Ztree->Branch("nVertex",   &ZnVertex);
  Ztree->Branch("HT",        &ZHT);
  Ztree->Branch("nJet",      &ZnJet);
	Ztree->Branch("nMu",       &ZnMu);
	Ztree->Branch("nGoodMu",   &ZnGoodMu);
	Ztree->Branch("JetPt",     &ZJetPt);
	Ztree->Branch("JetMatchPho",&ZJetMatchPho);
//*********** histo list **********************//
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
			  if(im->isMedium()){nGoodMu+=1;}
			}	
			p_nMu->Fill(nGoodMu);

			if(raw.nMu < 2)continue;

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
			if(hasLep && raw.nMu >= 2 ){
					if(passMETFilter(METFilter)){ 
						
						bool foundZ(false); std::vector<recoMuon>::iterator trailLep=Muon.begin();	
						for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++){
							if(im == signalLep)continue;
							if(!im->isMedium())continue;
							double llmass = (im->getP4()+signalLep->getP4()).M();
							if( (80 < llmass && llmass < 100)){
								foundZ = true;
								trailLep = im;
							}
						}	
					
						if(foundZ){	
							float deltaPhi = DeltaPhi(signalLep->getPhi(), METPhi);
							float MT = sqrt(2*MET*signalLep->getPt()*(1-std::cos(deltaPhi)));

							ZlepPt = signalLep->getPt();
							ZlepEta= signalLep->getEta();
							ZlepPhi= signalLep->getPhi();
							ZtrailPt = trailLep->getPt();
							ZtrailEta = trailLep->getEta();
							ZsigMT = MT;
							ZsigMET= MET;
							ZsigMETPhi = METPhi;
							ZdPhiLepMET = deltaPhi;
							ZdilepMass = (trailLep->getP4()+signalLep->getP4()).M();
							ZnVertex = nVtx;
							ZnMu = raw.nMu;
							ZnGoodMu = nGoodMu;

							ZnJet = 0;
							ZHT = 0;
							ZJetPt.clear();
							ZJetMatchPho.clear();
							for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
								if(!itJet->passSignalSelection())continue;
								if(DeltaR(itJet->getEta(), itJet->getPhi(), trailLep->getEta(),trailLep->getPhi()) <= 0.4)continue;	
								if(DeltaR(itJet->getEta(), itJet->getPhi(), signalLep->getEta(),signalLep->getPhi()) <= 0.4)continue;
								ZJetPt.push_back(itJet->getPt());
								ZnJet += 1;
								ZHT += itJet->getPt();
							}
							
							Ztree->Fill();
						}	

					}//MET Filter
			}//Candidate Filter


	}//loop on  events

outputfile->Write();
logfile.close();
}


