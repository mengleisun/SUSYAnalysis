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

#include "../../include/analysis_rawData.h"
#include "../../include/analysis_photon.h"
#include "../../include/analysis_muon.h"
#include "../../include/analysis_ele.h"
#include "../../include/analysis_tools.h"

void analysis_singleMu(){//main  

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

  char outputname[100] = "/uscms_data/d3/mengleis/plot_SingleMuonTrigger.root";
  ofstream logfile;
  logfile.open("/uscms_data/d3/mengleis/plot_SingleMuonTrigger.log"); 

  logfile << "analysis_mgTrigger()" << std::endl;

  RunType datatype(SingleMuon2016); 

  TChain* es = new TChain("ggNtuplizer/EventTree");
  es->Add("/uscmst1b_scratch/lpc1/3DayLifetime/mengleis/SingleMu_FebReminiAOD.root");
//  logfile << "root://cmseos.fnal.gov//store/user/msun/2016ggNtuple/skim-SingleMu_ReReco2016B_private.root"<< std::endl;

  //std::vector<float> *L1Index = 0;
  //std::vector<float> *L1Pt = 0;
  //std::vector<float> *L1Eta = 0;
  //std::vector<float> *L1Phi = 0;
  //es->SetBranchAddress("L1Index", &L1Index);
  //es->SetBranchAddress("L1Pt",    &L1Pt);
  //es->SetBranchAddress("L1Eta",   &L1Eta);
  //es->SetBranchAddress("L1Phi",   &L1Phi);

  TFile *outputfile = TFile::Open(outputname,"RECREATE");
  outputfile->cd();
 
  TTree *mgtree = new TTree("mgTree","mgTree");
  int   mg_run(0);
  float mg_muPt(0);
  float mg_muEta(0);
  bool  isTight;
  float mg_muMiniIso(0);
  float mg_dR(0); 
  int   mg_phofireL1;
  int   mg_mufireL1;
  bool  mg_mgfireHLT;
  bool  passL1;
  bool  passIsoMuHLT;
  bool  passTrkMuHLT;
 
  mgtree->Branch("runN",      &mg_run);
  mgtree->Branch("muPt",      &mg_muPt);
  mgtree->Branch("muEta",     &mg_muEta);
  mgtree->Branch("isTight",      &isTight);
  mgtree->Branch("muMiniIso", &mg_muMiniIso);
  mgtree->Branch("dR",        &mg_dR);
  mgtree->Branch("phofireL1", &mg_phofireL1);
  mgtree->Branch("mufireL1",  &mg_mufireL1); 
  mgtree->Branch("mgfireHLT", &mg_mgfireHLT);
  mgtree->Branch("passL1",    &passL1);
  mgtree->Branch("passIsoMuHLT",   &passIsoMuHLT);
  mgtree->Branch("passTrkMuHLT",   &passTrkMuHLT);

  TH1F *p_dimuon = new TH1F("p_dimuon","di-muon invmass; #mu#mu mass(GeV);",200,0,200);
  TProfile *p_MuHLTeff = new TProfile("p_MuHLTeff",";#mu P_{T} (GeV);",100,0,100);

  const unsigned nEvts = es->GetEntries(); 
  std::cout << "total=" << es->GetEntries() << std::endl;
  
  rawData raw(es, datatype);
  std::vector<recoPhoton> Photon;
  std::vector<recoMuon>   Muon;
  std::vector<recoEle>   Ele;
  float MET(0);

    for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
  
      if (ievt%10000==0) std::cout << " -- Processing event " << ievt << std::endl;

        raw.GetData(es, ievt);
        Photon.clear();
        Muon.clear();
        Ele.clear();
        for(int iPho(0); iPho < raw.nPho; iPho++){Photon.push_back(recoPhoton(raw, iPho));}
        for(int iMu(0); iMu < raw.nMu; iMu++){Muon.push_back(recoMuon(raw, iMu));}
        for(int iEle(0); iEle < raw.nEle; iEle++){Ele.push_back(recoEle(raw, iEle));}
        MET = raw.pfMET;

        if(!raw.passHLT())continue;
        if(raw.nMu<2)continue;
        mg_run = raw.run; 

        std::vector<std::vector<recoMuon>::iterator> tagMuVec;
        std::vector<recoMuon>::iterator probeMu;
        tagMuVec.clear();

        for(std::vector<recoMuon>::iterator itLeadMu = Muon.begin(); itLeadMu!= Muon.end(); ++itLeadMu){
          if(!itLeadMu->isLoose() || itLeadMu->getPt() < 25.0 )continue;
          if(!itLeadMu->passHLTSelection())continue;
          tagMuVec.push_back(itLeadMu);
        }

        for(unsigned itag(0); itag < tagMuVec.size(); itag++){
          std::vector<recoMuon>::iterator tagMu = tagMuVec[itag];
		  for(std::vector<recoMuon>::iterator itMu = Muon.begin(); itMu!= Muon.end(); itMu++){
            if(itMu == tagMu)continue;
			bool sameCharge = ((tagMu->isPosi() && itMu->isPosi()) || (!tagMu->isPosi() && !itMu->isPosi()));        
			if( DeltaR(tagMu->getEta(), tagMu->getPhi(), itMu->getEta(), itMu->getPhi()) < 0.8)continue;
            if(sameCharge|| !itMu->isMedium() || itMu->getMiniIso()>0.2)continue; 
			float dimuonmass = (tagMu->getP4()+itMu->getP4()).M(); 
			p_dimuon->Fill(dimuonmass);
			if(dimuonmass > 80 && dimuonmass < 100){
              probeMu = itMu;
              if(probeMu->fireSingleTrg(1))p_MuHLTeff->Fill(itMu->getPt(), 1);
              else p_MuHLTeff->Fill(itMu->getPt(), 0);
			  mg_muPt = probeMu->getPt();
			  mg_muEta= probeMu->getEta();
			  isTight = probeMu->isTight();
			  mg_muMiniIso=probeMu->getMiniIso();
			  mg_dR = DeltaR(tagMu->getEta(), tagMu->getPhi(), probeMu->getEta(), probeMu->getPhi());
              passL1= probeMu->fireL1Trg(10);
              passIsoMuHLT = probeMu->fireSingleTrg(1); 
			  passTrkMuHLT = probeMu->fireSingleTrg(19);
              mgtree->Fill();     
              break;          
            }
		  }
        }
    }

outputfile->Write();

}


