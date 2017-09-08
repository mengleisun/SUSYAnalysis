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

void test_elefakepho(){//main  

  gSystem->Load("../../../lib/libAnaClasses.so");

	TH1F *pixel_dR = new TH1F("pixel_dR","",100,0,1);

  RunType datatype(MC); 

  TChain* es = new TChain("ggNtuplizer/EventTree");
  TFileCollection fc("dum","","DY.txt");
  es->AddFileInfoList((TCollection*)fc.GetList());

  int   mcType = MCType::DYJet;
  if(datatype == MC && mcType == MCType::NOMC){std::cout << "wrong MC type" << std::endl; throw;} 

  const unsigned nEvts = es->GetEntries()/500; 
  std::cout << "Total event: " << nEvts << std::endl;
  
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
        
	      for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
	        if(itpho->getCalibEt() < 30 || !itpho->isEB())continue;
	        if(itpho->isLoose()){

						bool isEle(false);
						bool isPho(false);
						for(std::vector<mcData>::iterator it = MCData.begin(); it!= MCData.end(); it++){
							if(it->getEt() < 10)continue;
							if(DeltaR(it->getEta(), it->getPhi(), itpho->getEta(), itpho->getPhi()) < 0.1 && fabs(itpho->getEt() - it->getEt())/itpho->getEt() < 0.2){
								if(fabs(it->getPID()) == 11)isEle = true;
								if(fabs(it->getPID()) == 22)isPho = true;
							}
						}
	
              bool PixelVeto = itpho->PixelSeed()==0? true: false;
					
						if(!isEle && isPho && PixelVeto){
							double pixledR(3);
							for(std::vector<recoEle>::iterator ie = Ele.begin(); ie != Ele.end(); ie++){
									if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < pixledR)pixledR = DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi());

							}
							pixel_dR->Fill(pixledR);
						}
					}
        }
        
 
   }//loop on  events

	pixel_dR->Draw();
}


