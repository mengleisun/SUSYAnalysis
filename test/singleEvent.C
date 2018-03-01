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


void singleEvent(){//main 

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

  RunType datatype(DoubleEG2016); 
  TChain* es = new TChain("ggNtuplizer/EventTree");
	es->Add("root://cmseos.fnal.gov//store/group/lpcsusystealth/ggNtuple_leppho/FebReminiAOD/skim-DoubleEG_FebReminiAOD.root");	

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

	{//loop on entries

			raw.GetData(es, 70672403);
			MCData.clear();
			Photon.clear();
			Muon.clear();
			Ele.clear();
			JetCollection.clear();
			for(int iPho(0); iPho < raw.nPho; iPho++){Photon.push_back(recoPhoton(raw, iPho));}
			for(int iMu(0); iMu < raw.nMu; iMu++){Muon.push_back(recoMuon(raw, iMu));}
			for(int iEle(0); iEle < raw.nEle; iEle++){Ele.push_back(recoEle(raw, iEle));}
			for(int iJet(0); iJet < raw.nJet; iJet++){JetCollection.push_back(recoJet(raw, iJet));}
			MET = raw.pfMET;
			METPhi = raw.pfMETPhi;
			METFilter = raw.metFilters;
			nVtx = raw.nVtx;


			/******************************************************************************************************************************************************************************/
			/***********************************                                  Select Photon                                              **********************************************/
			bool hasPho(false);
			std::vector<recoPhoton>::iterator signalPho = Photon.begin();
			for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
				if(itpho->getR9() < 0.5)continue;
				if(!itpho->passHLTSelection())continue;
				if(!itpho->passBasicSelection())continue;
				bool passSigma = itpho->passSigma(1);
				bool passChIso = itpho->passChIso(1);
				bool PixelVeto = itpho->PixelSeed()==0? true: false;
				bool GSFveto(true);
				bool photonFSRVeto(true);
				bool eleFSRVeto(true);
				for(std::vector<recoEle>::iterator ie = Ele.begin(); ie != Ele.end(); ie++){
					if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) <= 0.02)GSFveto = false;
					if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.3 && ie->getEt()>2.0)photonFSRVeto=false;
					if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.3 && DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) > 0.02)eleFSRVeto=false;
				}
				for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++){
					if(DeltaR(itpho->getEta(), itpho->getPhi(), im->getEta(), im->getPhi()) < 0.3 && im->getEt()>2.0){
						photonFSRVeto=false;
						eleFSRVeto=false;
					}
				}
        // ****************  standard ID ************************************//
				if(!itpho->passSignalSelection())continue;
				std::cout << " pho " << itpho->getEt() << " pixel " << PixelVeto << " FSR " << photonFSRVeto << std::endl;
				if(GSFveto && PixelVeto && photonFSRVeto){
					if(!hasPho){
						hasPho=true;
						std::cout << "sigma = " << itpho->getSigma() << " R9 = " << itpho->getR9() << std::endl;	
						signalPho = itpho;
					}
				}
			}


			/******************************************************************************************************************************************************************************/
			/***********************************                                  Select Lepton                                              **********************************************/
			bool hasLep(false);
			std::vector<recoEle>::iterator signalLep = Ele.begin();
			for(std::vector<recoEle>::iterator itEle = Ele.begin(); itEle != Ele.end(); itEle++){
				if(itEle->getCalibPt() < 25)continue;
				if((itEle->isEB() && itEle->getR9() < 0.5) || (itEle->isEE() && itEle->getR9() < 0.8))continue;

				if(!itEle->passHLTSelection())continue;
				if(itEle->passSignalSelection()){
					std::cout << " ele " << itEle->getPt() << " eta " << itEle->getEta() << std::endl;
					if(!hasLep){
						hasLep=true; 
						signalLep = itEle;
						std::cout << "sigma = " << itEle->getSigma() << " R9 = " << itEle->getR9() << std::endl;	
					}
				}
			}

			if(hasPho && hasLep){
				double dRlepphoton = DeltaR(signalPho->getEta(), signalPho->getPhi(), signalLep->getEta(), signalLep->getPhi()); 
				if(dRlepphoton > 0.8){
					if(fabs((signalPho->getCalibP4()+signalLep->getCalibP4()).M() - 91.188) > 10.0){

						if(raw.passMETFilter(METFilter)){

							float deltaPhi = DeltaPhi(signalLep->getPhi(), METPhi);
							float MT = sqrt(2*MET*signalLep->getCalibPt()*(1-std::cos(deltaPhi)));

							std::cout << " un_phoEt = " << signalPho->getEt() << std::endl;
							std::cout << " phoEt = " << signalPho->getCalibEt() << std::endl;
							std::cout << " phoEta= " << signalPho->getEta() << std::endl;
							std::cout << " phoPhi= " << signalPho->getPhi() << std::endl;
							std::cout << " lepPt = " << signalLep->getCalibPt() << std::endl;
							std::cout << " lepEta= " << signalLep->getEta() << std::endl;
							std::cout << " lepPhi= " << signalLep->getPhi() << std::endl;
							std::cout << " sigMT = " << MT << std::endl;
							std::cout << " sigMET= " << MET << std::endl;
							std::cout << " sigMETPhi = " << METPhi << std::endl;
							std::cout << " dPhiLepMET = " << deltaPhi << std::endl; 
							std::cout << " nVertex = " << nVtx << std::endl; 
							std::cout << " dRPhoLep= " << dRlepphoton << std::endl;
						//	HT = 0;
						//	for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
						//		if(!itJet->passSignalSelection())continue;
						//		if(DeltaR(itJet->getEta(), itJet->getPhi(), signalPho->getEta(),signalPho->getPhi()) <= 0.4)continue;	
						//		if(DeltaR(itJet->getEta(), itJet->getPhi(), signalLep->getEta(),signalLep->getPhi()) <= 0.4)continue;
						//		nJet += 1;
						//		HT += itJet->getPt();
						//	}	
						int nBJet = 0;
						for(std::vector<recoJet>::iterator itJet = JetCollection.begin() ; itJet != JetCollection.end(); ++itJet){
							if(itJet->getPt() < 20)continue;
							if(itJet->isBJet())nBJet+=1;
							if(itJet->isBJet())std::cout << "b-jet " << itJet->getPt() << std::endl;
							std::cout << "jet " << itJet->getPt() << std::endl;
						}
						std::cout << "BJet " << nBJet << std::endl;

						}//MET Filter
					}// Z mass Filter
				}//dR filter
			}// ele + pho candidate
	 
	}//loop on  events

}


