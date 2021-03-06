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
#include "TProfile2D.h"
#include "TLorentzVector.h"
#include "TFileCollection.h"

#include "../../include/analysis_rawData.h"
#include "../../include/analysis_photon.h"
#include "../../include/analysis_muon.h"
#include "../../include/analysis_ele.h"
#include "../../include/analysis_tools.h"
#include "../../include/analysis_mcData.h"

void analysis_mgTriggerMC(){//main  

	gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

	char outputname[100] = "/uscms_data/d3/mengleis/plot_MuonTrigger_DY.root";
	ofstream logfile;
	logfile.open("/uscms_data/d3/mengleis/plot_MuonTrigger_DY.log"); 

	logfile << "analysis_mgTrigger()" << std::endl;

	RunType datatype(MC); 

	TChain* es = new TChain("ggNtuplizer/EventTree");
	es->Add("root://cmseos.fnal.gov//store/user/msun/MCSummer16/DYJetsToLL_M-50_NLO.root");
	
	TFile *outputfile = TFile::Open(outputname,"RECREATE");
	outputfile->cd();

	TTree *mgtree = new TTree("mgTree","mgTree");
	float mg_phoEt(0);
	float mg_phoEta(0);
	float mg_phoPhi(0);
	float mg_muPt(0);
	float mg_muMiniIso(0);
	float mg_dR(0);
	float mg_R9(0); 
	int   mg_phofireL1;
	int   mg_mufireL1;
	int   mg_phofireL1_2;
	int   mg_mufireL1_2;
	int   mg_phofireHLT;
	int   mg_mufireHLT;
	int   mg_phofireHLT2;
	int   mg_mufireHLT2;

	mgtree->Branch("phoEt",     &mg_phoEt);
	mgtree->Branch("phoEta",    &mg_phoEta);
	mgtree->Branch("phoPhi",    &mg_phoPhi);
	mgtree->Branch("muPt",      &mg_muPt);
	mgtree->Branch("muMiniIso", &mg_muMiniIso);
	mgtree->Branch("dR",        &mg_dR);
	mgtree->Branch("R9",        &mg_R9);
	mgtree->Branch("phofireL1", &mg_phofireL1);
	mgtree->Branch("mufireL1",  &mg_mufireL1); 
	mgtree->Branch("phofireL1_2",&mg_phofireL1_2);
	mgtree->Branch("mufireL1_2", &mg_mufireL1_2); 
	mgtree->Branch("phofireHLT", &mg_phofireHLT);
	mgtree->Branch("mufireHLT",  &mg_mufireHLT);
	mgtree->Branch("phofireHLT2",&mg_phofireHLT2);
	mgtree->Branch("mufireHLT2", &mg_mufireHLT2);

	TH1F *p_dimuon = new TH1F("p_dimuon","di-muon invmass; #mu#mu mass(GeV);",200,0,200);
	TH1F *p_invmass = new TH1F("p_invmass","mu#mu#gamma invmass; #mu#mu#gamma mass(GeV);",60,60,120);

	const unsigned nEvts = es->GetEntries(); 
	std::cout << "total=" << es->GetEntries() << std::endl;

	rawData raw(es, datatype);
	std::vector<recoPhoton> Photon;
	std::vector<recoMuon>   Muon;
	std::vector<recoEle>   Ele;
  std::vector<mcData>  MCData;
	float MET(0);

	for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
  
		if(ievt%10000==0) std::cout << " -- Processing event " << ievt << std::endl;

			raw.GetData(es, ievt);
			Photon.clear();
			Muon.clear();
			Ele.clear();
			MCData.clear();
			for(int iPho(0); iPho < raw.nPho; iPho++){Photon.push_back(recoPhoton(raw, iPho));}
			for(int iMu(0); iMu < raw.nMu; iMu++){Muon.push_back(recoMuon(raw, iMu));}
			for(int iEle(0); iEle < raw.nEle; iEle++){Ele.push_back(recoEle(raw, iEle));}
			if(datatype == MC)for(int iMC(0); iMC < raw.nMC; iMC++){MCData.push_back(mcData(raw, iMC));}
			MET = raw.pfMET;

			std::vector<std::vector<mcData>::iterator> mcMuonCollection;
			std::vector<std::vector<mcData>::iterator> mcPhotonCollection;
			mcMuonCollection.clear();
			mcPhotonCollection.clear();
			for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
				if(itMC->getEt() < 10)continue;		
				if(abs(itMC->getPID()) == 13 && abs(itMC->getMomPID()) == 23)mcMuonCollection.push_back(itMC);
				if(abs(itMC->getPID()) == 22 && abs(itMC->getMomPID()) == 13)mcPhotonCollection.push_back(itMC);
			}

			std::vector<std::vector<recoMuon>::iterator> probeMuVec;
			std::vector<std::vector<recoPhoton>::iterator> probePhoVec; 
			probeMuVec.clear();
			probePhoVec.clear();
			for(std::vector<recoMuon>::iterator itMu = Muon.begin(); itMu!= Muon.end(); ++itMu){
				if(itMu->isMedium() && itMu->getMiniIso() < 0.2 && itMu->getD0() < 0.05 && itMu->getDz() < 0.1){
					for(unsigned it(0); it < mcMuonCollection.size(); it++){
						std::vector<mcData>::iterator itMcMuon = mcMuonCollection[it];
						if(DeltaR(itMu->getEta(), itMu->getPhi(), itMcMuon->getEta(), itMcMuon->getPhi()) < 0.1)probeMuVec.push_back(itMu);
					}
				}
			}

			for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
				if(!itpho->isLoose())continue;
				if(!itpho->isEB())continue;
				bool PixelVeto = itpho->PixelSeed()==0? true: false;
				bool GSFveto(true);
				bool FSRVeto(true);
				for(std::vector<recoEle>::iterator ie = Ele.begin(); ie != Ele.end(); ie++){
					if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) <= 0.05)GSFveto = false;
					if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.3 && DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) > 0.05 && ie->getCalibPt()>2.0)FSRVeto=false;
				}
				for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++)
					if(DeltaR(itpho->getEta(), itpho->getPhi(), im->getEta(), im->getPhi()) < 0.3 && im->getPt()>2.0)FSRVeto=false;
				if(GSFveto && PixelVeto && FSRVeto){
					for(unsigned it(0); it < mcPhotonCollection.size(); it++){
						std::vector<mcData>::iterator itMcPho = mcPhotonCollection[it];
						if(DeltaR(itpho->getEta(), itpho->getPhi(), itMcPho->getEta(), itMcPho->getPhi()) < 0.1)probePhoVec.push_back(itpho);
					}
				}
			}

			for(int iMu(0); iMu < probeMuVec.size(); iMu++){
				for(int ipho(0); ipho < probePhoVec.size(); ipho++){

					std::vector<recoMuon>::iterator probeMu = probeMuVec[iMu]; 
					std::vector<recoPhoton>::iterator probePho = probePhoVec[ipho];
      
					mg_R9 = probePho->getR9();
					mg_phoEt = probePho->getEt();
					mg_phoEta= probePho->getEta();
					mg_phoPhi= probePho->getPhi();
					mg_muPt=probeMu->getPt();
					mg_muMiniIso=probeMu->getMiniIso();

					if(probePho->fireL1Trg(12))mg_phofireL1=12;
					else if(probePho->fireL1Trg(17))mg_phofireL1=17;
					else mg_phofireL1=0;

					if(probeMu->fireL1Trg(12))mg_mufireL1=12;
					else if(probeMu->fireL1Trg(17))mg_mufireL1=17;
					else mg_mufireL1=0;

					if(probePho->fireL1Trg(29))mg_phofireL1_2=29;
					else mg_phofireL1_2=0;

					if(probeMu->fireL1Trg(29))mg_mufireL1_2=29;
					else mg_mufireL1_2=0;

					if(probePho->fireDoubleTrg(28) || probePho->fireDoubleTrg(29))mg_phofireHLT=1;
					else mg_phofireHLT=0;

					if(probeMu->fireSingleTrg(2) || probeMu->fireSingleTrg(21))mg_mufireHLT=1;
					else mg_mufireHLT=0;

					if(probePho->fireDoubleTrg(30))mg_phofireHLT2=1;
					else mg_phofireHLT2=0;

					if(probeMu->fireSingleTrg(22))mg_mufireHLT2=1;
					else mg_mufireHLT2=0;

					mgtree->Fill();
				}
			}

    }

outputfile->Write();
}


