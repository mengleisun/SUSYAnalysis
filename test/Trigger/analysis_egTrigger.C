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


void analysis_egTrigger(){//main  

	gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

	char outputname[100] = "/uscms_data/d3/mengleis/plot_egTrigger_DY.root";
	ofstream logfile;
	logfile.open("/uscms_data/d3/mengleis/plot_egTrigger_DY.log");

	logfile << "analysis_egTrigger()" << std::endl;

	RunType datatype(MC);

	TChain* es = new TChain("ggNtuplizer/EventTree");
	TFileCollection fc("dum","","/uscmst1b_scratch/lpc1/3DayLifetime/mengleis/DY_1.txt");
	es->AddFileInfoList((TCollection*)fc.GetList());

	TFile *outputfile = TFile::Open(outputname,"NEW");
	outputfile->cd();

	logfile << "output: " << outputname <<  std::endl;

	float tagPt(0);
	float tagEta(0);
	float tagPhi(0);
	float tagR9(0);
	std::vector<int> mcPID;
	std::vector<float> mcEta;
	std::vector<float> mcPhi;
	std::vector<float> mcPt;
	std::vector<int> mcMomPID;
	std::vector<int> mcGMomPID;

	TTree *egtree = new TTree("egTree","egTree");
	float probePhoEt(0);
	float probePhoEta(0);
	float probePhoPhi(0);
	float probePhoR9(0);
	bool  probePhoMatchLeading;
	bool  probePhoMatchTrailing;
	float invmass;

	egtree->Branch("tagPt",                 &tagPt);
	egtree->Branch("tagEta",                &tagEta);
	egtree->Branch("tagPhi",                &tagPhi);
	egtree->Branch("tagR9",                 &tagR9);
	egtree->Branch("probeEt",            		&probePhoEt);
	egtree->Branch("probeEta",           		&probePhoEta);
	egtree->Branch("probePhi",           		&probePhoPhi);
	egtree->Branch("probeR9",            		&probePhoR9);
	egtree->Branch("probeMatchLeading",  		&probePhoMatchLeading);
	egtree->Branch("probeMatchTrailing", 		&probePhoMatchTrailing);
	egtree->Branch("invmass",               &invmass);
	egtree->Branch("mcPID",			   					&mcPID);
	egtree->Branch("mcEta",			   					&mcEta);
	egtree->Branch("mcPhi",			   					&mcPhi);
	egtree->Branch("mcPt",				   				&mcPt);
	egtree->Branch("mcMomPID",			   			&mcMomPID);
	egtree->Branch("mcGMomPID",		   				&mcGMomPID);

	TTree *eetree = new TTree("eeTree","eeTree");
	float probeElePt(0);
	float probeEleEta(0);
	float probeElePhi(0);
	float probeEleR9(0);
	bool  probeEleMatchLeading(false);
	bool  probeEleMatchTrailing;
	float diElectronMass;

	eetree->Branch("tagPt",                 &tagPt);
	eetree->Branch("tagEta",                &tagEta);
	eetree->Branch("tagPhi",                &tagPhi);
	eetree->Branch("tagR9",                 &tagR9);
	eetree->Branch("probeEt",            		&probeElePt);
	eetree->Branch("probeEta",           		&probeEleEta);
	eetree->Branch("probePhi",           		&probeElePhi);
	eetree->Branch("probeR9",            		&probeEleR9);
	eetree->Branch("probeMatchLeading",     &probeEleMatchLeading);
	eetree->Branch("probeMatchTrailing", 		&probeEleMatchTrailing);
	eetree->Branch("invmass",        				&diElectronMass);
	eetree->Branch("mcPID",			   					&mcPID);
	eetree->Branch("mcEta",			   					&mcEta);
	eetree->Branch("mcPhi",			   					&mcPhi);
	eetree->Branch("mcPt",				   				&mcPt);
	eetree->Branch("mcMomPID",			   			&mcMomPID);
	eetree->Branch("mcGMomPID",		   				&mcGMomPID);

	const unsigned nEvts = es->GetEntries(); 
	logfile << "total event "<< nEvts << std::endl;
 
	rawData raw(es, datatype);
	std::vector<mcData>  MCData;
	std::vector<recoPhoton> Photon;
	std::vector<recoMuon>   Muon;
	std::vector<recoEle>   Ele;
	float MET(0);

	for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries

		if (ievt%10000==0){
			std::cout << " -- Processing event " << ievt << std::endl;
			logfile 	<< " -- Processing event " << ievt << std::endl;
		}

		raw.GetData(es, ievt);
		Photon.clear();
		Muon.clear();
		Ele.clear();
		MCData.clear();
		if(datatype == MC)for(int iMC(0); iMC < raw.nMC; iMC++){MCData.push_back(mcData(raw, iMC));}
		for(int iPho(0); iPho < raw.nPho; iPho++){Photon.push_back(recoPhoton(raw, iPho));}
		for(int iMu(0); iMu < raw.nMu; iMu++){Muon.push_back(recoMuon(raw, iMu));}
		for(int iEle(0); iEle < raw.nEle; iEle++){Ele.push_back(recoEle(raw, iEle));}
		MET = raw.pfMET;

		if(!raw.passHLT())continue;
		if(raw.nPho <1)continue;
		if(MET>70)continue;
		if(((raw.HLTEleMuX >> 1) &1) ==0)continue;

		std::vector<std::vector<recoEle>::iterator> tagEleVec;
		tagEleVec.clear();
		for(std::vector<recoEle>::iterator itEle = Ele.begin(); itEle != Ele.end(); itEle++){
			if(itEle->getEt() < 30 || fabs(itEle->getEta())>2.1)continue;
			if(!itEle->passHLTSelection())continue;
			if(!itEle->fireTrgs(11))continue;
			if(!itEle->passSignalSelection())continue;
			tagEleVec.push_back(itEle);
		}

		for(unsigned itag(0); itag < tagEleVec.size(); itag++){
			std::vector<recoEle>::iterator tagEle = tagEleVec[itag];
			bool foundZ(false);
			bool tagMatchLead(false);
			bool foundtrailEle(false);
			std::vector<recoEle>::iterator trailEle;
			std::vector<recoPhoton>::iterator leadpho;
			std::vector<recoPhoton>::iterator phoMatchele;

			for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
				if(DeltaR(tagEle->getEta(), tagEle->getPhi(), itpho->getEta(), itpho->getPhi())<0.05){
					phoMatchele = itpho; 
					tagMatchLead = (itpho->fireDoubleTrg(5) || itpho->fireDoubleTrg(6));
					if(tagMatchLead){
						for(std::vector<recoEle>::iterator itEle = Ele.begin(); itEle != Ele.end(); itEle++){
							if(DeltaR(phoMatchele->getEta(), phoMatchele->getPhi(), itEle->getEta(), itEle->getPhi()) < 0.1)continue;
							if(!itEle->isMiniMedium())continue;
							if((itEle->isEB() && itEle->getR9() < 0.5) || (itEle->isEE() && itEle->getR9() < 0.8))continue;
							diElectronMass = (phoMatchele->getP4()+itEle->getP4()).M();
							if(!foundtrailEle && diElectronMass  > 60 && diElectronMass < 120){
								foundtrailEle = true;
								trailEle = itEle;
							}
						}
					}
					continue;
				} // check if the tag fires the leading leg
				if(!itpho->isEB() && !itpho->isEE())continue;
				if(itpho->getR9() < 0.5)continue;
				if(itpho->isLoose()){
					invmass = (itpho->getP4()+tagEle->getP4()).M();
					if(!foundZ && invmass > 60 && invmass < 120){ 
						foundZ = true; 
						leadpho = itpho;
						tagEle = tagEle;
					}//found tag+probe 
				}//check loose photon id
			}//loop on reco photons

			mcPID.clear();
			mcEta.clear();
			mcPhi.clear();
			mcPt.clear();
			mcMomPID.clear();
			mcGMomPID.clear();
			if(datatype == MC){
				for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
					if(itMC->getEt() < 5)continue;
					float mcdR(3.0);
					if(foundZ)mcdR = DeltaR(leadpho->getEta(),leadpho->getPhi(), itMC->getEta(), itMC->getPhi());
					float mcdRe= DeltaR(tagEle->getEta(),tagEle->getPhi(),itMC->getEta(), itMC->getPhi());
					float mcdRe2(3.0);
					if(foundtrailEle)mcdRe2 = DeltaR(trailEle->getEta(),trailEle->getPhi(),itMC->getEta(), itMC->getPhi());
					if(mcdR < 0.3 || mcdRe < 0.3 || mcdRe2 < 0.3){
						mcPID.push_back(itMC->getPID());
						mcMomPID.push_back(itMC->getMomPID());
						mcGMomPID.push_back(itMC->getGMomPID());
						mcEta.push_back(itMC->getEta());      
						mcPhi.push_back(itMC->getPhi());
						mcPt.push_back(itMC->getEt());
					}
				}
			}

			tagPt = tagEle->getPt(); 
			tagEta= tagEle->getEta();
			tagPhi= tagEle->getPhi();
			tagR9 = tagEle->getR9();
			if(foundZ){
				invmass = (leadpho->getP4()+tagEle->getP4()).M();
				probePhoEt = leadpho->getEt();
				probePhoEta= leadpho->getEta();
				probePhoPhi= leadpho->getPhi();
				probePhoR9 = leadpho->getR9();
				probePhoMatchLeading = (leadpho->fireDoubleTrg(5) || leadpho->fireDoubleTrg(6));
				probePhoMatchTrailing= (leadpho->fireDoubleTrg(1) || leadpho->fireDoubleTrg(2));
				egtree->Fill();
			}
			if(foundtrailEle){
				diElectronMass = (trailEle->getP4()+tagEle->getP4()).M();
				probeElePt=trailEle->getPt();
				probeEleEta=trailEle->getEta();
				probeElePhi=trailEle->getPhi();
				probeEleR9=trailEle->getR9();
				probeEleMatchTrailing=(trailEle->fireTrgs(21) || trailEle->fireTrgs(22));
				eetree->Fill();
			}


		}// loop over tag electrons

	}//loop on  events

	outputfile->Write();
	logfile.close();
}


