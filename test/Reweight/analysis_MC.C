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


void analysis_MC(){//main  

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

  char outputname[100] = "/uscms_data/d3/mengleis/test/resTree_MC_NLO130.root";
  ofstream logfile;
  logfile.open("/uscms_data/d3/mengleis/test/resTree_MC_NLO130.log"); 

  logfile << "analysis_mg()" << std::endl;

  RunType datatype(MC); 
  TChain* es = new TChain("ggNtuplizer/EventTree");
	//es->Add("root://cmseos.fnal.gov//store/user/msun/MCSummer16/ZGTo2LG_RunIISummer16MiniAODv2-TrancheIV_v6-v1.root");
	//es->Add("root://cmseos.fnal.gov//store/user/msun/MCSummer16/ZLLGJets_MonoPhoton_PtG-130.root");
	es->Add("root://cmseos.fnal.gov//store/user/msun/MCSummer16/ZGTo2LG_PtG-130_RunIISummer16MiniAODv2-TrancheIV_v6-v1.root");

  const unsigned nEvts = es->GetEntries(); 
	float MCweight = 0.14*35.8*1000.0/nEvts;	
  logfile << "Total event: " << nEvts << std::endl;
  std::cout << "Total event: " << nEvts << std::endl;
  logfile << "Output file: " << outputname << std::endl;

  TFile *outputfile = TFile::Open(outputname,"RECREATE");
  outputfile->cd();

//************ MC Tree **********************//
  TTree *MCtree = new TTree("MCTree","MCTree");
  std::vector<int> MC_mcPID;
  std::vector<float> MC_mcEta;
  std::vector<float> MC_mcPhi;
  std::vector<float> MC_mcPt;
  std::vector<int> MC_mcMomPID;
  std::vector<int> MC_mcGMomPID;
 
	MCtree->Branch("MCweight",  &MCweight); 
  MCtree->Branch("mcPID",    &MC_mcPID);
  MCtree->Branch("mcEta",    &MC_mcEta);
  MCtree->Branch("mcPhi",    &MC_mcPhi);
  MCtree->Branch("mcPt",     &MC_mcPt);
  MCtree->Branch("mcMomPID", &MC_mcMomPID);
  MCtree->Branch("mcGMomPID",&MC_mcGMomPID);
//************ reco Tree **********************//
  TTree *recoTree = new TTree("recoTree","recoTree");
  std::vector<float> recophoEt;
  std::vector<float> recophoEta;
  std::vector<float> recophoPhi;
  
  recoTree->Branch("phoEt",     &recophoEt);
  recoTree->Branch("phoEta",    &recophoEta);
  recoTree->Branch("phoPhi",    &recophoPhi);
//*********** histo list **********************//
  rawData raw(es, datatype);
  std::vector<mcData>  MCData;
  std::vector<recoPhoton> Photon;
  std::vector<recoMuon>   Muon;
  std::vector<recoEle>   Ele;
  std::vector<recoJet>   JetCollection;

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

			recophoEt.clear();
			recophoEta.clear();
			recophoPhi.clear();
			for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
					bool PixelVeto = itpho->PixelSeed()==0? true: false;
					if(!itpho->passSignalSelection())continue;
					if(PixelVeto){
						recophoEt.push_back(itpho->getEt());
						recophoEta.push_back(itpho->getEta());
						recophoPhi.push_back(itpho->getPhi());
					}
			}
			recoTree->Fill();


			MC_mcPID.clear();
			MC_mcEta.clear();
			MC_mcPhi.clear();
			MC_mcPt.clear();
			MC_mcMomPID.clear();
			MC_mcGMomPID.clear();
			for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
				if(itMC->getEt() < 1.0)continue;
				if(itMC->getPID() == 22){
					MC_mcPID.push_back(itMC->getPID());
					MC_mcMomPID.push_back(itMC->getMomPID());
					MC_mcGMomPID.push_back(itMC->getGMomPID());
					MC_mcEta.push_back(itMC->getEta());
					MC_mcPhi.push_back(itMC->getPhi());
					MC_mcPt.push_back(itMC->getEt());
				}
			}
			MCtree->Fill();


	}//loop on  events

outputfile->Write();
logfile.close();
}


