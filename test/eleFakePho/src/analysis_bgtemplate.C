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
#include "TLorentzVector.h"
#include "TFileCollection.h"

#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_tools.h"

void analysis_bgtemplate(){//main  

  gSystem->Load("../../../lib/libAnaClasses.so");

  char outputname[100] = "/uscms_data/d3/mengleis/Sep1/plot_bgtemplate_FullEcal.root";
  ofstream logfile;
  logfile.open("/uscms_data/d3/mengleis/Sep1/plot_bgtemplate_FullEcal.log"); 

  logfile << "analysis_bgtemplate()" << std::endl;
	logfile << "no FSR before push photon into collection" << std::endl;

  RunType datatype(SingleMuon2016); 

  TChain* es = new TChain("ggNtuplizer/EventTree");
	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_SingleMuon_FebReminiAOD_B.root");
	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_SingleMuon_FebReminiAOD_C.root");
	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_SingleMuon_FebReminiAOD_D.root");
	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_SingleMuon_FebReminiAOD_E.root");
	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_SingleMuon_FebReminiAOD_F.root");
	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_SingleMuon_FebReminiAOD_G.root");
	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/private_SingleMuon_FebReminiAOD_H.root");
  TFile *outputfile = TFile::Open(outputname,"NEW");
  outputfile->cd();

  int   tracks(0);
  int   nVertex(0); 
  int   mcType = MCType::NOMC;
  if(datatype == MC && mcType == MCType::NOMC){std::cout << "wrong MC type" << std::endl; throw;} 
  logfile << "mcType" << mcType << std::endl;
 
  TTree *mtree = new TTree("BGTree","BGTree");
  float tagEta_mg;
  float tagPhi_mg;
  float tagPt_mg;
  float probeEta_mg;
  float probePhi_mg;
  float probePt_mg;
  float probeUncalibPt_mg;
  float invmass_mg;
  bool  vetovalue_mg;
	bool  FSRveto_mg;

  mtree->Branch("tracks",              &tracks);
  mtree->Branch("nVertex",             &nVertex);
  mtree->Branch("tagEta",              &tagEta_mg);
  mtree->Branch("tagPhi",              &tagPhi_mg);
  mtree->Branch("tagPt",               &tagPt_mg);
  mtree->Branch("probeEta",            &probeEta_mg);
  mtree->Branch("probePhi",            &probePhi_mg);
  mtree->Branch("probePt",             &probePt_mg);
  mtree->Branch("probeUncalibPt",      &probeUncalibPt_mg);
  mtree->Branch("invmass",             &invmass_mg);
  mtree->Branch("vetovalue",           &vetovalue_mg);
	mtree->Branch("FSRveto",             &FSRveto_mg);
  
  const unsigned nEvts = es->GetEntries(); 
  logfile << "Total event: " << nEvts << std::endl;
  logfile << "Output file: " << outputname << std::endl;
  
  rawData raw(es, datatype);
  std::vector<recoPhoton> Photon;
  std::vector<recoMuon>   Muon;
  std::vector<recoEle>   Ele;
  float MET(0);
  float METPhi(0);
  int   ntrks(0);
  int   nvtx(0);

    for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
  
      if (ievt%1000000==0) std::cout << " -- Processing event " << ievt << std::endl;

        raw.GetData(es, ievt);
        Photon.clear();
        Muon.clear();
        Ele.clear();
        for(int iPho(0); iPho < raw.nPho; iPho++){Photon.push_back(recoPhoton(raw, iPho));}
        for(int iMu(0); iMu < raw.nMu; iMu++){Muon.push_back(recoMuon(raw, iMu));}
        for(int iEle(0); iEle < raw.nEle; iEle++){Ele.push_back(recoEle(raw, iEle));}
        MET = raw.pfMET;
        METPhi = raw.pfMETPhi;
        nvtx = raw.nVtx;

        tracks = ntrks;
        nVertex = nvtx;
        if(MET > 70.0)continue;
        if(!raw.passHLT())continue;

        std::vector< std::vector<recoMuon>::iterator > bgMuCollection;
        bgMuCollection.clear();
				for(std::vector<recoMuon>::iterator itMu = Muon.begin(); itMu != Muon.end(); itMu++){
					if(itMu->getPt() < 30)continue;
					if(fabs(itMu->getEta() > 2.1))continue;
					if(itMu->passSignalSelection()){
						bgMuCollection.push_back(itMu);
					}
				}

				std::vector< std::vector<recoPhoton>::iterator > bgPhoCollection;
				bgPhoCollection.clear();
				for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
					if(itpho->getCalibEt() < 30)continue;
					bool muFSRveto(true);
					for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++)
				 		if(DeltaR(itpho->getEta(), itpho->getPhi(), im->getEta(), im->getPhi()) < 0.3 && im->getEt()>2.0)muFSRveto=false;
					if(muFSRveto && itpho->isLoose()){
				 		bgPhoCollection.push_back(itpho);
				 	}
			 	}

     		for(unsigned iPho(0); iPho < bgPhoCollection.size(); iPho++){
          for(unsigned iMu(0); iMu < bgMuCollection.size(); iMu++){
            std::vector<recoMuon>::iterator bgMu = bgMuCollection[iMu];
            std::vector<recoPhoton>::iterator bgPho = bgPhoCollection[iPho];
       			bool PixelVeto = bgPho->PixelSeed()==0? true: false;
       			bool GSFveto(true);
						bool FSRVeto(true);
       			for(std::vector<recoEle>::iterator ire = Ele.begin(); ire != Ele.end(); ire++){
         			if(ire->getCalibEt() > 2.0 && DeltaR(bgPho->getEta(), bgPho->getPhi(), ire->getEta(), ire->getPhi()) < 0.02)GSFveto=false;
							if(DeltaR(bgPho->getEta(), bgPho->getPhi(), ire->getEta(), ire->getPhi()) < 0.3 && ire->getCalibEt()>2.0)FSRVeto=false;
       			}

						 tagEta_mg = bgMu->getEta();
						 tagPhi_mg = bgMu->getPhi();
						 tagPt_mg = bgMu->getPt();
						 probeEta_mg = bgPho->getEta();
						 probePhi_mg = bgPho->getPhi();
						 probePt_mg = bgPho->getCalibEt();
						 probeUncalibPt_mg = bgPho->getEt();
						 invmass_mg = (bgPho->getCalibP4()+bgMu->getP4()).M();
						 vetovalue_mg = (PixelVeto && GSFveto);
						 FSRveto_mg = FSRVeto;

						 mtree->Fill();
           }
     		}

  }//loop on  events

outputfile->Write();
}


