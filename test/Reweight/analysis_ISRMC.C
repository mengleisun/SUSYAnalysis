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

RunYear =2016;

void analysis_ISRMC(){//main  

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

  char outputname[100] = "/uscms_data/d3/mengleis/FullStatusOct/resTree_ISR_WWG.root";
  ofstream logfile;
  logfile.open("/uscms_data/d3/mengleis/FullStatusOct/resTree_ISR_WWG.log"); 

  logfile << "analysis_mg()" << std::endl;

  RunType datatype(MCMuonEG2016);
	bool  isMC(false);
	if(datatype == MC || datatype == MCDoubleEG2016 || datatype == MCMuonEG2016||  datatype == MCSingleElectron2016 || datatype == MCSingleMuon2016||  datatype == MCDoubleMuon2016 || datatype == MCMET2016)isMC=true;
  TChain* es = new TChain("ggNtuplizer/EventTree");
	//es->Add("root://cmseos.fnal.gov//store/user/msun/MCSummer16/TTGJets_RunIISummer16MiniAODv2-TrancheIV_v6_v1ANDext1.root");
	//es->Add("root://cmseos.fnal.gov//store/user/msun/MCSummer16/TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8.root");
	//es->Add("root://cmseos.fnal.gov//store/user/msun/MCSummer16/ZGTo2LG_RunIISummer16MiniAODv2-TrancheIV_v6-v1.root");
	es->Add("root://cmseos.fnal.gov//store/user/msun/MCSummer16/WWG_RunIISummer16MiniAODv2-TrancheIV_v6_ext1.root");
	//es->Add("root://cmseos.fnal.gov//store/user/msun/MCSummer16/WZG_RunIISummer16MiniAODv2-TrancheIV_v6.root");
  //int mcType = MCType::TTG;
  //int mcType = MCType::TT;
	//int mcType = MCType::ZGInclusive; 
	int mcType = MCType::WWG;
	//int mcType = MCType::WZG;


  const unsigned nEvts = es->GetEntries(); 
	float PUweight(1);

  logfile << "mcType" << mcType << std::endl;
  float crosssection = MC_XS[mcType];
	float MCweight = crosssection*35.87*1000.0/nEvts;

	logfile << "crosssection = " << crosssection << std::endl;
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
	float ZbosonPt(0);
  std::vector<int> Z_mcPID;
  std::vector<float> Z_mcEta;
  std::vector<float> Z_mcPhi;
  std::vector<float> Z_mcPt;
  std::vector<int> Z_mcMomPID;
  std::vector<int> Z_mcGMomPID;
	std::vector<int> Z_mcStatus;
 
	Ztree->Branch("MCweight",  &MCweight); 
	Ztree->Branch("PUweight",  &PUweight);
  Ztree->Branch("mcType",    &mcType);
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
	Ztree->Branch("bosonPt",     &ZbosonPt);
  Ztree->Branch("mcPID",    &Z_mcPID);
  Ztree->Branch("mcEta",    &Z_mcEta);
  Ztree->Branch("mcPhi",    &Z_mcPhi);
  Ztree->Branch("mcPt",     &Z_mcPt);
  Ztree->Branch("mcMomPID", &Z_mcMomPID);
  Ztree->Branch("mcGMomPID",&Z_mcGMomPID);
	Ztree->Branch("mcStatus", &Z_mcStatus);

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
			if(RunYear==2016)PUweight = getPUESF16(nVtx);
			if(RunYear==2017)PUweight = getPUESF17(nVtx);
			if(RunYear==2018)PUweight = getPUESF18(nVtx);

			nTotal+=1;
			if(!raw.passHLT())continue;
			if(raw.nGoodVtx < 1)continue;
			npassHLT+=1;

			int nGoodMu(0);
			for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++){
				if(im->isMedium() && im->getMiniIso() < 0.2 &&  im->getPt() > 15){nGoodMu+=1;}
			}

			if(nGoodMu != 2 || raw.nPho <1)continue;

			ZbosonPt = 0;	
      for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
       switch(itMC->getPID()){
         case 23: ZbosonPt = itMC->getPt();  break;
         case 24: ZbosonPt = itMC->getPt();  break;
         case -24: ZbosonPt = itMC->getPt(); break;
         default: break;
       }
      }

			bool hasPho(false);
			std::vector<recoPhoton>::iterator signalPho = Photon.begin();
			for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
				if(itpho->getR9() < 0.5)continue;
				if(!itpho->passHLTSelection())continue;
				if(itpho->fireDoubleTrg(28) || itpho->fireDoubleTrg(29) || itpho->fireDoubleTrg(30)){
					if(!itpho->passBasicSelection())continue;
					bool passSigma = itpho->passSigma(1);
					bool passChIso = itpho->passChIso(1);
					bool PixelVeto = itpho->PixelSeed()==0? true: false;
					bool GSFveto(true);
					bool FSRVeto(true);
					for(std::vector<recoEle>::iterator ie = Ele.begin(); ie != Ele.end(); ie++){
						if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) <= 0.02)GSFveto = false;
						if(DeltaR(itpho->getEta(), itpho->getPhi(), ie->getEta(), ie->getPhi()) < 0.3)FSRVeto=false;
					}
					for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++)
						if(DeltaR(itpho->getEta(), itpho->getPhi(), im->getEta(), im->getPhi()) < 0.3 && im->getEt()>2.0)FSRVeto=false;

					if(!itpho->passSignalSelection())continue;
					if(GSFveto && PixelVeto && FSRVeto){
						if(!hasPho){
							hasPho=true;
							npassPho +=1;
							signalPho = itpho;
						}
					}

				}
			}

			bool hasLep(false);
			std::vector<recoMuon>::iterator signalLep = Muon.begin();
			std::vector< std::vector<recoMuon>::iterator > proxyLepCollection;
			proxyLepCollection.clear();
			for(std::vector<recoMuon>::iterator itMu = Muon.begin(); itMu != Muon.end(); itMu++){
				if(itMu->getPt() < 25)continue;
				if(!itMu->passHLTSelection())continue;
				if(itMu->passSignalSelection()){
					proxyLepCollection.push_back(itMu);
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
					if(raw.passMETFilter(METFilter)){ 
						
						bool foundZG(false); 
						std::vector<recoMuon>::iterator trailLep=Muon.begin();	
						for(std::vector<recoMuon>::iterator im = Muon.begin(); im != Muon.end(); im++){
							if(im == signalLep)continue;
							if(!im->isMedium() || im->getMiniIso() > 0.2)continue;
							double llmass = (im->getP4()+signalLep->getP4()).M();
							double llgmass = (im->getP4()+signalLep->getP4()+signalPho->getP4()).M();
							if( 80 < llmass && llmass < 100){ 
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

						  Z_mcPID.clear();
						  Z_mcEta.clear();
						  Z_mcPhi.clear();
						  Z_mcPt.clear();
						  Z_mcMomPID.clear();
						  Z_mcGMomPID.clear();
							Z_mcStatus.clear();
						  for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){
						 	  if(itMC->getEt() < 1.0)continue;
							  float mcdRmu = DeltaR(signalLep->getEta(), signalLep->getPhi(), itMC->getEta(), itMC->getPhi());
						 	  float mcdR = DeltaR(signalPho->getEta(), signalPho->getPhi(), itMC->getEta(), itMC->getPhi());
						 	  if(mcdR < 0.3 || mcdRmu < 0.3){
						 		  Z_mcPID.push_back(itMC->getPID());
						 		  Z_mcMomPID.push_back(itMC->getMomPID());
						 		  Z_mcGMomPID.push_back(itMC->getGMomPID());
						 		  Z_mcEta.push_back(itMC->getEta());
						 		  Z_mcPhi.push_back(itMC->getPhi());
						 		  Z_mcPt.push_back(itMC->getEt());
									Z_mcStatus.push_back(itMC->getStatus());
						 	  }
						  }
							
							Ztree->Fill();
						}	

					}//MET Filter
				}//dR Filter
			}//Candidate Filter


	}//loop on  events

outputfile->Write();
logfile.close();
}


