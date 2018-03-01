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
#include "TRandom3.h"
#include "TFileCollection.h"

#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_tools.h"
#include "../../../include/analysis_mcData.h"

void analysis_elefakepho(){//main  

  gSystem->Load("../../../lib/libAnaClasses.so");

  char outputname[100] = "/uscms_data/d3/mengleis/FullStatusOct/plot_elefakepho_DYTnP_dR05.root";
  ofstream logfile;
  logfile.open("/uscms_data/d3/mengleis/FullStatusOct/plot_elefakepho_DYTnP_dR05.log"); 

  logfile << "analysis_elefakepho()" << std::endl;

  RunType datatype(MC); 

  TChain* es = new TChain("ggNtuplizer/EventTree");
  es->Add("root://cmseos.fnal.gov//store/user/msun/MCSummer16/skim-DYJetsToLL_M-50_nlo.root");
  
  TFile *outputfile = TFile::Open(outputname,"NEW");
  outputfile->cd();

  int   tracks(0);
  int   nVertex(0); 
  int   mcType = MCType::DYLL50;
  if(datatype == MC && mcType == MCType::NOMC){std::cout << "wrong MC type" << std::endl; throw;} 
  logfile << "mcType" << mcType << std::endl;

  TTree *etree = new TTree("FakeRateTree","FakeRateTree");
  float tagEta_bothcount;
  float tagPhi_bothcount;
  float tagPt_bothcount;
  float tagUncalibPt_bothcount;
  float probeEta_bothcount;
  float probePhi_bothcount;
  float probePt_bothcount;
  float probeUncalibPt_bothcount;
  float invmass_bothcount;
  float invmassUncalib_bothcount;
  bool  vetovalue_bothcount;
	bool  FSRveto_bothcount;
  std::vector<int> mcPID_bothcount;
  std::vector<float> mcEta_bothcount;
  std::vector<float> mcPhi_bothcount;
  std::vector<float> mcPt_bothcount;
  std::vector<int> mcMomPID_bothcount;
  std::vector<int> mcGMomPID_bothcount;

  etree->Branch("tracks",              &tracks);
  etree->Branch("nVertex",             &nVertex);
  etree->Branch("tagEta",              &tagEta_bothcount);
  etree->Branch("tagPhi",              &tagPhi_bothcount);
  etree->Branch("tagPt",               &tagPt_bothcount);
  etree->Branch("tagUncalibPt",        &tagUncalibPt_bothcount);
  etree->Branch("probeEta",            &probeEta_bothcount);
  etree->Branch("probePhi",            &probePhi_bothcount);
  etree->Branch("probePt",             &probePt_bothcount);
  etree->Branch("probeUncalibPt",      &probeUncalibPt_bothcount);
  etree->Branch("invmass",             &invmass_bothcount);
  etree->Branch("invmassUncalib",      &invmassUncalib_bothcount);
  etree->Branch("vetovalue",           &vetovalue_bothcount);
	etree->Branch("FSRveto",       &FSRveto_bothcount);
  etree->Branch("mcPID",			   &mcPID_bothcount);
  etree->Branch("mcEta",			   &mcEta_bothcount);
  etree->Branch("mcPhi",			   &mcPhi_bothcount);
  etree->Branch("mcPt",				   &mcPt_bothcount);
  etree->Branch("mcMomPID",			   &mcMomPID_bothcount);
  etree->Branch("mcGMomPID",		   &mcGMomPID_bothcount);
  etree->Branch("mcType",			   &mcType);

  TTree *rantree = new TTree("FakeRateRandomTree","FakeRateRandomTree");
  float tagEta_random;
  float tagPhi_random;
  float tagPt_random;
  float tagUncalibPt_random;
  float probeEta_random;
  float probePhi_random;
  float probePt_random;
  float probeUncalibPt_random;
  float invmass_random;
  float invmassUncalib_random;
  bool  vetovalue_random;
	bool  FSRveto_random;
  std::vector<int> mcPID_random;
  std::vector<float> mcEta_random;
  std::vector<float> mcPhi_random;
  std::vector<float> mcPt_random;
  std::vector<int> mcMomPID_random;
  std::vector<int> mcGMomPID_random;

  rantree->Branch("tracks",              &tracks);
  rantree->Branch("nVertex",             &nVertex);
  rantree->Branch("tagEta",              &tagEta_random);
  rantree->Branch("tagPhi",              &tagPhi_random);
  rantree->Branch("tagPt",               &tagPt_random);
  rantree->Branch("tagUncalibPt",        &tagUncalibPt_random);
  rantree->Branch("probeEta",            &probeEta_random);
  rantree->Branch("probePhi",            &probePhi_random);
  rantree->Branch("probePt",             &probePt_random);
  rantree->Branch("probeUncalibPt",      &probeUncalibPt_random);
  rantree->Branch("invmass",             &invmass_random);
  rantree->Branch("invmassUncalib",      &invmassUncalib_random);
  rantree->Branch("vetovalue",           &vetovalue_random);
	rantree->Branch("FSRveto",         &FSRveto_random);
  rantree->Branch("mcPID",			     &mcPID_random);
  rantree->Branch("mcEta",			     &mcEta_random);
  rantree->Branch("mcPhi",			     &mcPhi_random);
  rantree->Branch("mcPt",				 &mcPt_random);
  rantree->Branch("mcMomPID",			 &mcMomPID_random);
  rantree->Branch("mcGMomPID",		     &mcGMomPID_random);
  rantree->Branch("mcType",			     &mcType);

  const unsigned nEvts = es->GetEntries(); 
  std::cout << "Total event: " << nEvts << std::endl;
  logfile << "Total event: " << nEvts << std::endl;
  logfile << "Output file: " << outputname << std::endl;
  
  rawData raw(es, datatype);
  std::vector<mcData>  MCData;
  std::vector<recoPhoton> Photon;
  std::vector<recoMuon>   Muon;
  std::vector<recoEle>   Ele;
  float MET(0);
  float METPhi(0);
  int   ntrks(0);
  int   nvtx(0);

  TRandom3 ran(0);

    std::cout << "total: " << nEvts << std::endl;
    for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
  
      if (ievt%10000==0) std::cout << " -- Processing event " << ievt << std::endl;

        raw.GetData(es, ievt);
        MCData.clear();
        Photon.clear();
        Muon.clear();
        Ele.clear();
        if(datatype == MC)for(int iMC(0); iMC < raw.nMC; iMC++){MCData.push_back(mcData(raw, iMC));}
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
				if(((raw.HLTEleMuX >> 1) &1) ==0)continue;

        std::vector<std::vector<recoEle>::iterator> ElectronCollection;
        ElectronCollection.clear();
        for(std::vector<recoEle>::iterator itEle = Ele.begin(); itEle != Ele.end(); itEle++){
           if(itEle->getCalibEt() < 30 || fabs(itEle->getEta())>2.1)continue;
           if(!itEle->passHLTSelection())continue;
					 if(!itEle->fireTrgs(11))continue;
           if(itEle->passSignalSelection())ElectronCollection.push_back(itEle);
        }
        
        std::vector<std::vector<recoPhoton>::iterator> signalPho;
        std::vector<bool> PhoPixelVeto;
        std::vector<bool> PhoEleVeto;
				std::vector<bool> PhoFSRVeto;
        signalPho.clear(); 
        PhoPixelVeto.clear();
        PhoEleVeto.clear();
				PhoFSRVeto.clear();
        if(ElectronCollection.size() > 0){
	      for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
	        if(itpho->getCalibEt() < 30)continue;
	        if(itpho->isLoose()){
              bool PixelVeto = itpho->PixelSeed()==0? true: false;
              bool GSFveto(true);
		      bool FSRVeto(true);
			  for(std::vector<recoEle>::iterator ie = Ele.begin(); ie != Ele.end(); ie++){
				 //if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.02)GSFveto = false;
				 if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.05)GSFveto = false;
				 if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.3 && ie->getCalibEt()>2.0)FSRVeto=false;
			  }
			  for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++)
				 if(DeltaR(itpho->getEta(), itpho->getPhi(), im->getEta(), im->getPhi()) < 0.3 && im->getEt()>2.0)FSRVeto=false;
					signalPho.push_back(itpho);
					PhoPixelVeto.push_back(PixelVeto);
					PhoEleVeto.push_back(GSFveto);
					PhoFSRVeto.push_back(FSRVeto);
            }
	      }
        }
        
        int nTagEle = ElectronCollection.size();
	    int randomTag(0);
	    if(nTagEle>1)randomTag = ran.Integer(nTagEle);
	    for(unsigned iTag(0); iTag < nTagEle; iTag++){ 
	      for(unsigned iProbe(0); iProbe < signalPho.size(); iProbe++){
	        double dRTagProbe = DeltaR(ElectronCollection[iTag]->getEta(), ElectronCollection[iTag]->getPhi(), signalPho[iProbe]->getEta(), signalPho[iProbe]->getPhi());
	        double dETagProbe = fabs(ElectronCollection[iTag]->getCalibEt() - signalPho[iProbe]->getCalibEt())/signalPho[iProbe]->getCalibEt();
	        double InvMass = (ElectronCollection[iTag]->getCalibP4() + signalPho[iProbe]->getCalibP4()).M();
	        if(dRTagProbe < 0.05 && dETagProbe < 0.1)continue;
	        if(InvMass < 20 || InvMass > 160)continue;
              
	        if(iTag == randomTag){
			  tagEta_random = ElectronCollection[iTag]->getEta();
			  tagPhi_random = ElectronCollection[iTag]->getPhi();
			  tagPt_random = ElectronCollection[iTag]->getCalibPt();
			  tagUncalibPt_random = ElectronCollection[iTag]->getPt();
			  probeEta_random = signalPho[iProbe]->getEta();
			  probePhi_random = signalPho[iProbe]->getPhi();
			  probePt_random = signalPho[iProbe]->getCalibEt();
			  probeUncalibPt_random = signalPho[iProbe]->getEt();
			  invmass_random = InvMass;
			  invmassUncalib_random = (ElectronCollection[iTag]->getP4() + signalPho[iProbe]->getP4()).M();
			  vetovalue_random = (PhoPixelVeto[iProbe] && PhoEleVeto[iProbe]);
				FSRveto_random = PhoFSRVeto[iProbe];
              if(datatype == MC){
				mcPID_random.clear();
				mcEta_random.clear();
				mcPhi_random.clear();
				mcPt_random.clear();
				mcMomPID_random.clear();
				mcGMomPID_random.clear();

				for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
				  if(itMC->getEt() < 5)continue;
				  float mcdR = DeltaR(signalPho[iProbe]->getEta(),signalPho[iProbe]->getPhi(), itMC->getEta(), itMC->getPhi());
				  float mcdRe= DeltaR(ElectronCollection[iTag]->getEta(),ElectronCollection[iTag]->getPhi(),itMC->getEta(), itMC->getPhi());
				  if(mcdR < 0.3 || mcdRe < 0.3){
					mcPID_random.push_back(itMC->getPID());
					mcMomPID_random.push_back(itMC->getMomPID());
					mcGMomPID_random.push_back(itMC->getGMomPID());
					mcEta_random.push_back(itMC->getEta());      
					mcPhi_random.push_back(itMC->getPhi());
					mcPt_random.push_back(itMC->getEt());
				  }
				}
              }
		
              rantree->Fill();
            }
			tagEta_bothcount = ElectronCollection[iTag]->getEta();
			tagPhi_bothcount = ElectronCollection[iTag]->getPhi();
			tagPt_bothcount = ElectronCollection[iTag]->getCalibPt();
			tagUncalibPt_bothcount = ElectronCollection[iTag]->getPt();
			probeEta_bothcount = signalPho[iProbe]->getEta();
			probePhi_bothcount = signalPho[iProbe]->getPhi();
			probePt_bothcount = signalPho[iProbe]->getCalibEt();
			probeUncalibPt_bothcount = signalPho[iProbe]->getEt();
			invmass_bothcount = InvMass;
			invmassUncalib_bothcount = (ElectronCollection[iTag]->getP4() + signalPho[iProbe]->getP4()).M();
			vetovalue_bothcount = (PhoPixelVeto[iProbe] && PhoEleVeto[iProbe]);
			FSRveto_bothcount = PhoFSRVeto[iProbe];

            if(vetovalue_bothcount == false){
              bool isTag(false);
              for(unsigned iTag(0); iTag < nTagEle; iTag++){
                double testdR = DeltaR(ElectronCollection[iTag]->getEta(), ElectronCollection[iTag]->getPhi(), signalPho[iProbe]->getEta(), signalPho[iProbe]->getPhi());
                if(testdR < 0.05)isTag = true;
              }
            }

			if(datatype == MC){
			  mcPID_bothcount.clear();
			  mcEta_bothcount.clear();
			  mcPhi_bothcount.clear();
			  mcPt_bothcount.clear();
			  mcMomPID_bothcount.clear();
			  mcGMomPID_bothcount.clear();

			  for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
				if(itMC->getEt() < 5)continue;
				float mcdR = DeltaR(signalPho[iProbe]->getEta(),signalPho[iProbe]->getPhi(), itMC->getEta(), itMC->getPhi());
				float mcdRe= DeltaR(ElectronCollection[iTag]->getEta(),ElectronCollection[iTag]->getPhi(),itMC->getEta(), itMC->getPhi());
				if(mcdR < 0.3 || mcdRe < 0.3){
				  mcPID_bothcount.push_back(itMC->getPID());
				  mcMomPID_bothcount.push_back(itMC->getMomPID());
				  mcGMomPID_bothcount.push_back(itMC->getGMomPID());
				  mcEta_bothcount.push_back(itMC->getEta());      
				  mcPhi_bothcount.push_back(itMC->getPhi());
				  mcPt_bothcount.push_back(itMC->getEt());
				}
			  }
			}
            etree->Fill();
	   }
	 }
           

   }//loop on  events

outputfile->Write();
outputfile->Close(); 
logfile.close();
}


