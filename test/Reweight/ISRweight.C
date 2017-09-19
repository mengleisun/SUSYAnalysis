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
#include "TProfile2D.h"

#include "../../include/analysis_rawData.h"
#include "../../include/analysis_photon.h"
#include "../../include/analysis_muon.h"
#include "../../include/analysis_ele.h"
#include "../../include/analysis_mcData.h"
#include "../../include/analysis_tools.h"
#include "../../include/analysis_fakes.h"
#include "../../include/analysis_scalefactor.h"

bool doEB=false;

void ISRweight(){//main 

	int channelType = 2;
  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");


	esfScaleFactor  objectESF;

	Double_t plotEtBins[]={35,50,80,100,125,150,200,250,300,400,450,800};
	Double_t plotPtBins[]={25, 35,50,100,150,200,250,300,400,600,800};
	TH1F *p_phoEt_data     = new TH1F("p_phoEt_data","",11,plotEtBins);
	TH1F *p_phoEt_up     	 = new TH1F("p_phoEt_up","",11,plotEtBins);
	TH1F *p_phoEt_down     = new TH1F("p_phoEt_down","",11,plotEtBins);
	TH1F *p_JetPt_data     = new TH1F("p_JetPt_data","",11,plotEtBins);
	TH1F *p_JetPt_up       = new TH1F("p_JetPt_up","",11,plotEtBins);
	TH1F *p_JetPt_down     = new TH1F("p_JetPt_down","",11,plotEtBins);
	TH1F *p_JetPt_total    = new TH1F("p_JetPt_total","",11,plotEtBins);
	TH1F *p_llmass_data    = new TH1F("p_llmass_data","",100,30,130);
	TH1F *p_llmass_up      = new TH1F("p_llmass_up","",100,30,130);
	TH1F *p_llmass_down    = new TH1F("p_llmass_down","",100,30,130);
	TH1F *p_phoEt_data_highEt   = new TH1F("p_phoEt_data_highEt","",11,plotEtBins);
	TH1F *p_JetPt_data_highEt   = new TH1F("p_JetPt_data_highEt","",11,plotEtBins);
	TH1F *p_llmass_data_highEt  = new TH1F("p_llmass_data_highEt","",100,30,130);

	TH1F *p_phoEt_ZG     = new TH1F("p_phoEt_ZG","",11,plotEtBins);
	TH1F *p_JetPt_ZG     = new TH1F("p_JetPt_ZG","",11,plotEtBins);
	TH1F *p_JetPt_ZG_up  = new TH1F("p_JetPt_ZG_up","",11,plotEtBins);
	TH1F *p_JetPt_ZG_down= new TH1F("p_JetPt_ZG_down","",11,plotEtBins);
	TH1F *p_llmass_ZG    = new TH1F("p_llmass_ZG","",100,30,130);
	TH1F *p_phoEt_ZG_highEt   = new TH1F("p_phoEt_ZG_highEt","",11,plotEtBins);
	TH1F *p_JetPt_ZG_highEt   = new TH1F("p_JetPt_ZG_highEt","",11,plotEtBins);
	TH1F *p_llmass_ZG_highEt  = new TH1F("p_llmass_ZG_highEt","",100,30,130);

	TH1F *p_phoEt_rare     = new TH1F("p_phoEt_rare","",11,plotEtBins);
	TH1F *p_JetPt_rare     = new TH1F("p_JetPt_rare","",11,plotEtBins);
	TH1F *p_llmass_rare    = new TH1F("p_llmass_rare","",100,30,130);
	TH1F *p_phoEt_rare_highEt   = new TH1F("p_phoEt_rare_highEt","",11,plotEtBins);
	TH1F *p_JetPt_rare_highEt   = new TH1F("p_JetPt_rare_highEt","",11,plotEtBins);
	TH1F *p_llmass_rare_highEt  = new TH1F("p_llmass_rare_highEt","",100,30,130);

	TProfile *p_scalefactor  = new TProfile("p_scalefactor","p_scalefactor",50,-2.5,2.5);
//************ Signal Tree **********************//
  TChain *tree = new TChain("ZTree");
  tree->Add("/uscms_data/d3/mengleis/test/resTree_ISR_data.root");
  float phoEt(0);
  float phoEta(0);
  float phoPhi(0);
  float lepPt(0);
  float lepEta(0);
  float lepPhi(0);
	float trailPt(0);
  float sigMT(0);
  float sigMET(0);
  float sigMETPhi(0);
  float dPhiLepMET(0);
  int   nVertex(0);
  float dRPhoLep(0);
  float HT(0);
  float nJet(0);
	float threeMass(0);
	float dilepMass(0);
	float JetPt=0;

  tree->SetBranchAddress("phoEt",     &phoEt);
  tree->SetBranchAddress("phoEta",    &phoEta);
  tree->SetBranchAddress("phoPhi",    &phoPhi);
  tree->SetBranchAddress("lepPt",     &lepPt);
  tree->SetBranchAddress("lepEta",    &lepEta);
  tree->SetBranchAddress("lepPhi",    &lepPhi);
  tree->SetBranchAddress("trailPt",   &trailPt);
  tree->SetBranchAddress("sigMT",     &sigMT);
  tree->SetBranchAddress("sigMET",    &sigMET);
  tree->SetBranchAddress("sigMETPhi", &sigMETPhi);
  tree->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
  tree->SetBranchAddress("nVertex",   &nVertex);
  tree->SetBranchAddress("dRPhoLep",  &dRPhoLep);
  tree->SetBranchAddress("HT",        &HT);
  tree->SetBranchAddress("nJet",      &nJet);
	tree->SetBranchAddress("threeMass", &threeMass);
	tree->SetBranchAddress("dilepMass", &dilepMass);
	tree->SetBranchAddress("JetPt",     &JetPt);

  for(unsigned ievt(0); ievt<tree->GetEntries(); ++ievt){//loop on entries
		tree->GetEntry(ievt);
		
		if(doEB && fabs(phoEta) > 1.4442)continue;
		//else if(!doEB && (fabs(phoEta) < 1.56 || fabs(phoEta) > 2.4))continue;

		if(phoEt > 799)phoEt = 799;
		if(JetPt > 799)JetPt = 799;
		if(dRPhoLep < 0.8)continue;
		if(dilepMass < 80 || dilepMass > 100)continue;
		double weight = 1;
		double weightup = 1;
		double weightdown = 1; 
		if(dilepMass > 80 && dilepMass < 100){

			if(doEB){
				if(phoEt > 35 && phoEt < 50)weight = 1-0.20;
				else if(phoEt > 50 && phoEt < 100)weight = 1-0.13;
				else if(phoEt > 100 && phoEt < 150)weight = 1-0.068;
				else if(phoEt > 150 && phoEt < 200)weight = 1-0.067;
				else if(phoEt > 200 && phoEt < 300)weight = 1-0.024;
				else weight=1;	
	
				if(phoEt > 35 && phoEt < 40)weightup = 1-0.22 + 0.07;
				else if(phoEt > 40 && phoEt < 50)weightup = 1- 0.18 + 0.04;
				else if(phoEt > 50 && phoEt < 60)weightup = 1- 0.17 + 0.13;
				else if(phoEt > 60 && phoEt < 70)weightup = 1- 0.13 + 0.05;
				else if(phoEt > 70 && phoEt < 80)weightup = 1- 0.11 + 0.11;
				else if(phoEt > 80 && phoEt < 100)weightup = 1- 0.056 + 0.034;
				else if(phoEt > 100 && phoEt < 150)weightup = 1- 0.068 + 0.037;
				else if(phoEt > 150 && phoEt < 200)weightup = 1- 0.067 + 0.07;
				else if(phoEt > 200 && phoEt < 250)weightup = 1- 0.033 + 0.075;
				else weightup=1;	
	
				if(phoEt > 35 && phoEt < 40)weightdown = 1-0.22 - 0.07;
				else if(phoEt > 40 && phoEt < 50)weightdown = 1- 0.18 - 0.04;
				else if(phoEt > 50 && phoEt < 60)weightdown = 1- 0.17 - 0.13;
				else if(phoEt > 60 && phoEt < 70)weightdown = 1- 0.13 - 0.05;
				else if(phoEt > 70 && phoEt < 80)weightdown = 1- 0.11 - 0.11;
				else if(phoEt > 80 && phoEt < 100)weightdown = 1- 0.056 - 0.034;
				else if(phoEt > 100 && phoEt < 150)weightdown = 1- 0.068 - 0.037;
				else if(phoEt > 150 && phoEt < 200)weightdown = 1- 0.067 - 0.07;
				else if(phoEt > 200 && phoEt < 250)weightdown = 1- 0.033 - 0.075;
				else weightdown=1;
			}	
			else{
				if(phoEt > 35 && phoEt < 40)weight = 1-0.145349 ;
				else if(phoEt > 40 && phoEt < 50)weight = 1- 0.172 ;
				else if(phoEt > 50 && phoEt < 60)weight = 1- 0.11;
				else if(phoEt > 60 && phoEt < 70)weight = 1- 0.26 ;
				else if(phoEt > 70 && phoEt < 80)weight = 1- 0.22 ;
				else if(phoEt > 80 && phoEt < 100)weight = 1- 0.20 ;
				else if(phoEt > 100 && phoEt < 150)weight = 1- 0.09 ;
				else if(phoEt > 150 && phoEt < 200)weight = 1-0.03;
				else if(phoEt > 200 && phoEt < 250)weight = 1-0.03;
				else weight=1;	
	
				if(phoEt > 35 && phoEt < 40)weightup = 1-0.145349 ;
				else if(phoEt > 40 && phoEt < 50)weightup = 1- 0.147485 ;
				else if(phoEt > 50 && phoEt < 60)weightup = 1- 0.0932913;
				else if(phoEt > 60 && phoEt < 70)weightup = 1- 0.191783 ;
				else if(phoEt > 70 && phoEt < 80)weightup = 1- 0.17711 ;
				else if(phoEt > 80 && phoEt < 100)weightup = 1- 0.17148 ;
				else if(phoEt > 100 && phoEt < 150)weightup = 1- 0.0335751 ;
				else if(phoEt > 150 && phoEt < 200)weightup = 1;
				else if(phoEt > 200 && phoEt < 250)weightup = 1;
				else weightup=1;	
				if(phoEt > 35 && phoEt < 40)weightdown = 1-0.145349 ;
				else if(phoEt > 40 && phoEt < 50)weightdown = 1- 0.147485 ;
				else if(phoEt > 50 && phoEt < 60)weightdown = 1- 0.0932913;
				else if(phoEt > 60 && phoEt < 70)weightdown = 1- 0.191783 ;
				else if(phoEt > 70 && phoEt < 80)weightdown = 1- 0.17711 ;
				else if(phoEt > 80 && phoEt < 100)weightdown = 1- 0.17148 ;
				else if(phoEt > 100 && phoEt < 150)weightdown = 1- 0.0335751 ;
				else if(phoEt > 150 && phoEt < 200)weightdown = 1;
				else if(phoEt > 200 && phoEt < 250)weightdown = 1;
				else weightdown=1;	
			}
		}

	 	p_JetPt_data->Fill(JetPt, weight);
		p_phoEt_data->Fill(phoEt, weight); 
		p_llmass_data->Fill(dilepMass, weight);

		p_JetPt_up->Fill(JetPt, weightup);	
		p_phoEt_up->Fill(phoEt, weightup); 
		p_llmass_up->Fill(dilepMass, weightup);

		p_JetPt_down->Fill(JetPt, weightdown);	
		p_phoEt_down->Fill(phoEt, weightdown); 
		p_llmass_down->Fill(dilepMass, weightdown);

		if(phoEt > 145){
			p_llmass_data_highEt->Fill(dilepMass, weight);
			p_JetPt_data_highEt->Fill(JetPt, weight);
			p_phoEt_data_highEt->Fill(phoEt, weight);
		}
	}//loop on  events




//************ Signal Tree **********************//
  TChain *ZGtree = new TChain("ZTree");
  ZGtree->Add("/uscms_data/d3/mengleis/test/resTree_ISR_ZG.root");
  float ZG_phoEt(0);
  float ZG_phoEta(0);
  float ZG_phoPhi(0);
  float ZG_lepPt(0);
  float ZG_lepEta(0);
  float ZG_lepPhi(0);
	float ZG_trailPt(0);
	float ZG_trailEta(0);
	float ZG_trailPhi(0);
	float ZG_bosonPt(0);
  float ZG_sigMT(0);
  float ZG_sigMET(0);
  float ZG_sigMETPhi(0);
  float ZG_dPhiLepMET(0);
  int   ZG_nVertex(0);
  float ZG_dRPhoLep(0);
  float ZG_dRPhoTrail(0);
  float ZG_HT(0);
  float ZG_nJet(0);
	float ZG_JetPt(0);
	float ZG_threeMass(0);
	float ZG_dilepMass(0);
  std::vector<int>   *ZG_mcPID=0;
  std::vector<float> *ZG_mcEta=0;
  std::vector<float> *ZG_mcPhi=0;
  std::vector<float> *ZG_mcPt=0;
  std::vector<int>   *ZG_mcMomPID=0;
  std::vector<int>   *ZG_mcGMomPID=0;

  ZGtree->SetBranchAddress("phoEt",     &ZG_phoEt);
  ZGtree->SetBranchAddress("phoEta",    &ZG_phoEta);
  ZGtree->SetBranchAddress("phoPhi",    &ZG_phoPhi);
  ZGtree->SetBranchAddress("lepPt",     &ZG_lepPt);
  ZGtree->SetBranchAddress("lepEta",    &ZG_lepEta);
  ZGtree->SetBranchAddress("lepPhi",    &ZG_lepPhi);
	ZGtree->SetBranchAddress("trailPt",   &ZG_trailPt);
	ZGtree->SetBranchAddress("trailEta",  &ZG_trailEta);
	ZGtree->SetBranchAddress("trailPhi",  &ZG_trailPhi);
	ZGtree->SetBranchAddress("bosonPt",   &ZG_bosonPt);
  ZGtree->SetBranchAddress("sigMT",     &ZG_sigMT);
  ZGtree->SetBranchAddress("sigMET",    &ZG_sigMET);
  ZGtree->SetBranchAddress("sigMETPhi", &ZG_sigMETPhi);
  ZGtree->SetBranchAddress("dPhiLepMET",&ZG_dPhiLepMET);
  ZGtree->SetBranchAddress("nVertex",   &ZG_nVertex);
  ZGtree->SetBranchAddress("dRPhoLep",  &ZG_dRPhoLep);
  ZGtree->SetBranchAddress("dRPhoTrail",&ZG_dRPhoTrail);
  ZGtree->SetBranchAddress("HT",        &ZG_HT);
  ZGtree->SetBranchAddress("nJet",      &ZG_nJet);
	ZGtree->SetBranchAddress("threeMass", &ZG_threeMass);
	ZGtree->SetBranchAddress("dilepMass", &ZG_dilepMass);
	ZGtree->SetBranchAddress("JetPt",     &ZG_JetPt);
  ZGtree->SetBranchAddress("mcPID",    &ZG_mcPID);
  ZGtree->SetBranchAddress("mcEta",    &ZG_mcEta);
  ZGtree->SetBranchAddress("mcPhi",    &ZG_mcPhi);
  ZGtree->SetBranchAddress("mcPt",     &ZG_mcPt);
  ZGtree->SetBranchAddress("mcMomPID", &ZG_mcMomPID);
  ZGtree->SetBranchAddress("mcGMomPID",&ZG_mcGMomPID);

	double totalevent(0);
	double reweightevent(0);

  for(unsigned ievt(0); ievt<ZGtree->GetEntries(); ++ievt){//loop on entries
		ZGtree->GetEntry(ievt);

		if(doEB && fabs(ZG_phoEta) > 1.4442)continue;
		//if(!doEB && (fabs(ZG_phoEta)< 1.56 || fabs(ZG_phoEta) > 2.4))continue;
		double weight = 35.8*1000*117.8/14372399;

		double scalefactor(0);
		if(channelType == 2){
			scalefactor = objectESF.getMuonESF(ZG_lepPt,ZG_lepEta)*objectESF.getPhotonESF(ZG_phoEt,ZG_phoEta)*objectESF.getMuonEGTRGESF(ZG_phoEt, ZG_lepPt);
		}

		if(ZG_phoEt > 799)ZG_phoEt = 799;
		if(ZG_JetPt > 799)ZG_JetPt = 799;
		if(ZG_dRPhoLep < 0.8)continue;
		if(ZG_dilepMass < 80 || ZG_dilepMass >100)continue;
		weight = weight*scalefactor;

		float PUweight = getPUESF(ZG_nVertex);
		weight = weight*PUweight;
 
		bool   isTruePho(false);
		double mindR(0.3);
		unsigned phoIndex(0);
		for(unsigned iMC(0); iMC< ZG_mcPID->size(); iMC++){
			double dR1 = DeltaR((*ZG_mcEta)[iMC], (*ZG_mcPhi)[iMC], ZG_phoEta, ZG_phoPhi);
			double dE1 = fabs((*ZG_mcPt)[iMC] - ZG_phoEt)/ZG_phoEt;
			if(dR1 < mindR && dE1 < 0.2){mindR=dR1; phoIndex=iMC;}
		}
		if(mindR < 0.1){
			if((*ZG_mcPID)[phoIndex] == 22 && (fabs((*ZG_mcMomPID)[phoIndex]) == 23 || fabs((*ZG_mcMomPID)[phoIndex]) == 13 || fabs((*ZG_mcMomPID)[phoIndex])==999)){
				isTruePho=true;
			}
		}

		p_JetPt_ZG->Fill(ZG_JetPt, weight);	
		p_JetPt_ZG_up->Fill(ZG_JetPt, weight);	
		p_JetPt_ZG_down->Fill(ZG_JetPt, weight);	
		p_phoEt_ZG->Fill(ZG_phoEt, weight); 
		p_llmass_ZG->Fill(ZG_dilepMass, weight);
	}//loop on  events


//************ Signal Tree **********************//
  TChain *raretree = new TChain("ZTree");
  raretree->Add("/uscms_data/d3/mengleis/test/resTree_ISR_TTG.root");
  raretree->Add("/uscms_data/d3/mengleis/test/resTree_ISR_WWG.root");
  raretree->Add("/uscms_data/d3/mengleis/test/resTree_ISR_WZG.root");
  raretree->Add("/uscms_data/d3/mengleis/test/resTree_ISR_TT.root");
	float MCweight(0);
  float rare_phoEt(0);
  float rare_phoEta(0);
  float rare_phoPhi(0);
  float rare_lepPt(0);
  float rare_lepEta(0);
  float rare_lepPhi(0);
  float rare_sigMT(0);
  float rare_sigMET(0);
  float rare_sigMETPhi(0);
  float rare_dPhiLepMET(0);
  int   rare_nVertex(0);
  float rare_dRPhoLep(0);
  float rare_HT(0);
  float rare_nJet(0);
	float rare_threeMass(0);
	float rare_dilepMass(0);
	float rare_JetPt(0);
  std::vector<int>   *rare_mcPID=0;
  std::vector<float> *rare_mcEta=0;
  std::vector<float> *rare_mcPhi=0;
  std::vector<float> *rare_mcPt=0;
  std::vector<int>   *rare_mcMomPID=0;
  std::vector<int>   *rare_mcGMomPID=0;

  raretree->SetBranchAddress("MCweight",  &MCweight);
  raretree->SetBranchAddress("phoEt",     &rare_phoEt);
  raretree->SetBranchAddress("phoEta",    &rare_phoEta);
  raretree->SetBranchAddress("phoPhi",    &rare_phoPhi);
  raretree->SetBranchAddress("lepPt",     &rare_lepPt);
  raretree->SetBranchAddress("lepEta",    &rare_lepEta);
  raretree->SetBranchAddress("lepPhi",    &rare_lepPhi);
  raretree->SetBranchAddress("sigMT",     &rare_sigMT);
  raretree->SetBranchAddress("sigMET",    &rare_sigMET);
  raretree->SetBranchAddress("sigMETPhi", &rare_sigMETPhi);
  raretree->SetBranchAddress("dPhiLepMET",&rare_dPhiLepMET);
  raretree->SetBranchAddress("nVertex",   &rare_nVertex);
  raretree->SetBranchAddress("dRPhoLep",  &rare_dRPhoLep);
  raretree->SetBranchAddress("HT",        &rare_HT);
  raretree->SetBranchAddress("nJet",      &rare_nJet);
	raretree->SetBranchAddress("JetPt",     &rare_JetPt);
	raretree->SetBranchAddress("threeMass", &rare_threeMass);
	raretree->SetBranchAddress("dilepMass", &rare_dilepMass);
  raretree->SetBranchAddress("mcPID",    &rare_mcPID);
  raretree->SetBranchAddress("mcEta",    &rare_mcEta);
  raretree->SetBranchAddress("mcPhi",    &rare_mcPhi);
  raretree->SetBranchAddress("mcPt",     &rare_mcPt);
  raretree->SetBranchAddress("mcMomPID", &rare_mcMomPID);
  raretree->SetBranchAddress("mcGMomPID",&rare_mcGMomPID);

  for(unsigned ievt(0); ievt<raretree->GetEntries(); ++ievt){//loop on entries
		raretree->GetEntry(ievt);
		
		if(doEB && fabs(rare_phoEta) > 1.4442)continue;
		//else if(!doEB && (fabs(rare_phoEta) < 1.56|| fabs(rare_phoEta) > 2.4))continue;

		double weight = MCweight; 
		double scalefactor(0);
		if(channelType == 1){
			scalefactor = objectESF.getElectronESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getegPhotonTRGESF(phoEt,phoEta)*objectESF.getElectronTRGESF(lepPt,lepEta);
		}
		if(channelType == 2){
			scalefactor = objectESF.getMuonESF(rare_lepPt,rare_lepEta)*objectESF.getPhotonESF(rare_phoEt,rare_phoEta)*objectESF.getMuonEGTRGESF(rare_phoEt, rare_lepPt);
		}

		if(rare_phoEt > 799)rare_phoEt = 799;
		if(rare_JetPt > 799)rare_JetPt = 799;
		if(rare_dRPhoLep < 0.8)continue;

		if(rare_dilepMass < 80 || rare_dilepMass > 100)continue;
		weight = weight*scalefactor;
		float PUweight = getPUESF(rare_nVertex);
		weight = weight*PUweight;

	//	bool   isTruePho(false);
	//	double mindR(0.3);
	//	unsigned phoIndex(0);
	//	for(unsigned iMC(0); iMC< rare_mcPID->size(); iMC++){
	//		double dR1 = DeltaR((*rare_mcEta)[iMC], (*rare_mcPhi)[iMC], rare_phoEta, rare_phoPhi);
	//		double dE1 = fabs((*rare_mcPt)[iMC] - rare_phoEt)/rare_phoEt;
	//		if(dR1 < mindR && dE1 < 0.2){mindR=dR1; phoIndex=iMC;}
	//	}
	//	if(mindR < 0.1){
	//		if((*rare_mcPID)[phoIndex] == 22 && (fabs((*rare_mcMomPID)[phoIndex])<25 || fabs((*rare_mcMomPID)[phoIndex])==999))isTruePho=true;
	//	}

	//	if(!isTruePho && mindR < 0.1)std::cout << (*rare_mcPID)[phoIndex]  << " " <<  (*rare_mcMomPID)[phoIndex] << " " << (*rare_mcPt)[phoIndex] << std::endl;
	//	else if(mindR >= 0.3)std::cout << "no match " << rare_phoEt << std::endl;


		p_phoEt_rare->Fill(rare_phoEt, weight); 
		p_JetPt_rare->Fill(rare_JetPt, weight);
		p_llmass_rare->Fill(rare_dilepMass, weight);

		if(rare_phoEt > 145){
			p_phoEt_rare_highEt->Fill(rare_phoEt, weight); 
			p_JetPt_rare_highEt->Fill(rare_JetPt, weight);
			p_llmass_rare_highEt->Fill(rare_dilepMass, weight);
		}
	}//loop on  events
	
  TChain *highEttree = new TChain("ZTree");
  highEttree->Add("/uscms_data/d3/mengleis/test/resTree_ISR_ZG130_NLO.root");
	float highEtweight(0);
  float highEt_phoEt(0);
  float highEt_phoEta(0);
  float highEt_phoPhi(0);
  float highEt_lepPt(0);
  float highEt_lepEta(0);
  float highEt_lepPhi(0);
  float highEt_sigMT(0);
  float highEt_sigMET(0);
  float highEt_sigMETPhi(0);
  float highEt_dPhiLepMET(0);
  int   highEt_nVertex(0);
  float highEt_dRPhoLep(0);
  float highEt_HT(0);
  float highEt_nJet(0);
	float highEt_threeMass(0);
	float highEt_dilepMass(0);
	float highEt_JetPt(0);
  std::vector<int>   *highEt_mcPID=0;
  std::vector<float> *highEt_mcEta=0;
  std::vector<float> *highEt_mcPhi=0;
  std::vector<float> *highEt_mcPt=0;
  std::vector<int>   *highEt_mcMomPID=0;
  std::vector<int>   *highEt_mcGMomPID=0;

  highEttree->SetBranchAddress("MCweight",  &highEtweight);
  highEttree->SetBranchAddress("phoEt",     &highEt_phoEt);
  highEttree->SetBranchAddress("phoEta",    &highEt_phoEta);
  highEttree->SetBranchAddress("phoPhi",    &highEt_phoPhi);
  highEttree->SetBranchAddress("lepPt",     &highEt_lepPt);
  highEttree->SetBranchAddress("lepEta",    &highEt_lepEta);
  highEttree->SetBranchAddress("lepPhi",    &highEt_lepPhi);
  highEttree->SetBranchAddress("sigMT",     &highEt_sigMT);
  highEttree->SetBranchAddress("sigMET",    &highEt_sigMET);
  highEttree->SetBranchAddress("sigMETPhi", &highEt_sigMETPhi);
  highEttree->SetBranchAddress("dPhiLepMET",&highEt_dPhiLepMET);
  highEttree->SetBranchAddress("nVertex",   &highEt_nVertex);
  highEttree->SetBranchAddress("dRPhoLep",  &highEt_dRPhoLep);
  highEttree->SetBranchAddress("HT",        &highEt_HT);
  highEttree->SetBranchAddress("nJet",      &highEt_nJet);
	highEttree->SetBranchAddress("JetPt",     &highEt_JetPt);
	highEttree->SetBranchAddress("threeMass", &highEt_threeMass);
	highEttree->SetBranchAddress("dilepMass", &highEt_dilepMass);
  highEttree->SetBranchAddress("mcPID",    &highEt_mcPID);
  highEttree->SetBranchAddress("mcEta",    &highEt_mcEta);
  highEttree->SetBranchAddress("mcPhi",    &highEt_mcPhi);
  highEttree->SetBranchAddress("mcPt",     &highEt_mcPt);
  highEttree->SetBranchAddress("mcMomPID", &highEt_mcMomPID);
  highEttree->SetBranchAddress("mcGMomPID",&highEt_mcGMomPID);

  for(unsigned ievt(0); ievt<highEttree->GetEntries(); ++ievt){//loop on entries
		highEttree->GetEntry(ievt);

		double weight = highEtweight;	
		double scalefactor(0);
		if(channelType == 2){
			scalefactor = objectESF.getMuonESF(highEt_lepPt,highEt_lepEta)*objectESF.getPhotonESF(highEt_phoEt,highEt_phoEta)*objectESF.getMuonEGTRGESF(highEt_phoEt, highEt_lepPt);
		}

		if(doEB && fabs(highEt_phoEta) > 1.4442)continue;
		//else if(!doEB && (fabs(highEt_phoEta) < 1.56 || fabs(highEt_phoEta) > 2.4))continue;

		if(highEt_phoEt > 799)highEt_phoEt = 799;
		if(highEt_JetPt > 799)highEt_JetPt = 799;
		if(highEt_dRPhoLep < 0.8)continue;
		if(highEt_phoEt < 145)continue;
		if(highEt_dilepMass < 80 || highEt_dilepMass > 100)continue;
		weight = weight*scalefactor;
		float PUweight = getPUESF(highEt_nVertex);
		weight = weight*PUweight;

		bool   isTruePho(false);
		double mindR(0.3);
		unsigned phoIndex(0);
		for(unsigned iMC(0); iMC< highEt_mcPID->size(); iMC++){
			double dR1 = DeltaR((*highEt_mcEta)[iMC], (*highEt_mcPhi)[iMC], highEt_phoEta, highEt_phoPhi);
			double dE1 = fabs((*highEt_mcPt)[iMC] - highEt_phoEt)/highEt_phoEt;
			if(dR1 < mindR && dE1 < 0.2){mindR=dR1; phoIndex=iMC;}
		}
		if(mindR < 0.1){
			if((*highEt_mcPID)[phoIndex] == 22 && (fabs((*highEt_mcMomPID)[phoIndex])<25 || fabs((*highEt_mcMomPID)[phoIndex])==999))isTruePho=true;
		}
		if(!isTruePho)continue;


		p_phoEt_ZG_highEt->Fill(highEt_phoEt, weight); 
		p_JetPt_ZG_highEt->Fill(highEt_JetPt, weight);
		p_llmass_ZG_highEt->Fill(highEt_dilepMass, weight);
	}//loop on  events

	TCanvas *can_phoEt     = new TCanvas("can_phoEt",       "can_phoEt", 600,600); 
	TCanvas *can_JetPt        = new TCanvas("can_JetPt",          "can_JetPt",600,600); 
	TCanvas *can_JetPt_up     = new TCanvas("can_JetPt_up",       "can_JetPt_up",600,600); 
	TCanvas *can_JetPt_down   = new TCanvas("can_JetPt_down",     "can_JetPt_down",600,600); 
	TCanvas *can_JetPt_alter  = new TCanvas("can_JetPt_alter",    "can_JetPt_alter",600,600); 
	TCanvas *can_llmass   = new TCanvas("can_llmass",     "can_llmass", 600,600);   


	p_llmass_ZG->Add(p_llmass_rare);
	p_llmass_ZG_highEt->Add(p_llmass_rare_highEt);
	float scalefactor = p_llmass_data->Integral(50,70)/p_llmass_ZG->Integral(50,70); 
	std::cout << "scale = " << scalefactor<<std::endl;
	float scalefactorup = p_llmass_up->Integral(50,70)/p_llmass_ZG->Integral(50,70); 
	std::cout << "scale = " << scalefactorup<<std::endl;
	float scalefactordown = p_llmass_down->Integral(50,70)/p_llmass_ZG->Integral(50,70); 
	std::cout << "scale = " << scalefactordown<<std::endl;
	float scalefactorhighEt = p_llmass_data_highEt->Integral(50,70)/p_llmass_ZG_highEt->Integral(50,70);
	std::cout << "highEt scale = " << scalefactorhighEt <<std::endl;

	can_phoEt->cd();
	TPad *phoEt_pad1 = new TPad("phoEt_pad1", "phoEt_pad1", 0, 0.3, 1, 1.0);
	phoEt_pad1->SetBottomMargin(0.1);
	phoEt_pad1->Draw();  
	phoEt_pad1->cd();  
	phoEt_pad1->SetLogy();
	p_phoEt_data->Draw();
	p_phoEt_ZG->Add(p_phoEt_rare);
	p_phoEt_ZG->SetLineColor(kRed);
	p_phoEt_ZG->Scale(scalefactor);
	p_phoEt_ZG->Draw("same");
	p_phoEt_rare->SetLineColor(kGreen);
	p_phoEt_rare->Scale(scalefactor);
	p_phoEt_rare->Draw("same");
	TLegend *leg=new TLegend(0.6,0.9,0.75,0.9);
	leg->AddEntry(p_phoEt_data,"observed");
	leg->AddEntry(p_phoEt_ZG,"ZG");
	leg->AddEntry(p_phoEt_rare,"other");
	leg->Draw("same");

	can_phoEt->cd();
	TPad *phoEt_pad2 = new TPad("phoEt_pad2", "phoEt_pad2", 0, 0.05, 1, 0.25);
	phoEt_pad2->Draw();
	phoEt_pad2->cd();
	TLine *flatratio_eventcount = new TLine(0,1,600,1);
	TH1F *ratio_phoEt = (TH1F*)p_phoEt_ZG->Clone("ratio_phoEt");
	ratio_phoEt->Divide(p_phoEt_data);
	ratio_phoEt->Draw();
	flatratio_eventcount->Draw("same");
	
	double ISRwgt_norm[11];
	double ISRwgt_up[11];
	double ISRwgt_down[11];
	double ISRwgt_alter[11];

	can_JetPt->cd();
	gPad->SetLogy();
	gStyle->SetOptStat(0);
	p_JetPt_data->GetXaxis()->SetTitle("ISR P_{T} (GeV)");
	p_JetPt_data->SetMinimum(0.1);
	p_JetPt_data->SetLineColor(kBlack);
	p_JetPt_data->SetMarkerStyle(20);
	p_JetPt_data->SetMarkerColor(kBlack);
	p_JetPt_data->Draw("P");
	p_JetPt_ZG->Add(p_JetPt_rare);
	p_JetPt_ZG->SetLineColor(6);
	p_JetPt_ZG->SetFillStyle(1001);
	p_JetPt_ZG->SetFillColor(6);
	p_JetPt_ZG->Scale(scalefactor);
	p_JetPt_ZG->Draw("hist same");
	p_JetPt_rare->SetLineColor(8);
	p_JetPt_rare->SetFillStyle(1001);
	p_JetPt_rare->SetFillColor(8);
	p_JetPt_rare->Scale(scalefactor);
	p_JetPt_rare->Draw("hist same");
	p_JetPt_data->Draw("P same");
	TLegend *leg_HT =  new TLegend(0.6,0.75,0.9,0.9);
	leg_HT->SetFillStyle(0);
	gStyle->SetLegendBorderSize(1);
	gStyle->SetLegendFillColor(0);
	leg_HT->AddEntry(p_JetPt_data, "Data");
	leg_HT->AddEntry(p_JetPt_ZG,  "ZG");
	leg_HT->AddEntry(p_JetPt_rare, "tt,WWG,WZG");   
	leg_HT->Draw("same");
	can_JetPt->SaveAs("ISRweight.pdf");
	for(int i(1); i<=11; i++){
		std::cout << "norm ratio " << i << " " << p_JetPt_data->GetBinContent(i)/p_JetPt_ZG->GetBinContent(i) 
							<< " stat " <<   p_JetPt_data->GetBinError(i)/p_JetPt_ZG->GetBinContent(i)
              << std::endl;
		ISRwgt_norm[i-1] = p_JetPt_data->GetBinContent(i)/p_JetPt_ZG->GetBinContent(i);
	}

	can_JetPt_up->cd();
	gPad->SetLogy();
	p_JetPt_up->Draw();
	p_JetPt_ZG_up->Add(p_JetPt_rare);
	p_JetPt_ZG_up->Scale(scalefactorup);
	p_JetPt_ZG_up->Draw("same");
	p_JetPt_rare->Draw("same");
	std::cout << std::endl;
	for(int i(1); i<=11; i++){
		std::cout << "up   ratio " << i << " " << p_JetPt_up->GetBinContent(i)/p_JetPt_ZG_up->GetBinContent(i) << std::endl;
		ISRwgt_up[i-1] = p_JetPt_up->GetBinContent(i)/p_JetPt_ZG_up->GetBinContent(i);	
	}

	can_JetPt_down->cd();
	gPad->SetLogy();
	p_JetPt_down->Draw();
	p_JetPt_ZG_down->Add(p_JetPt_rare);
	p_JetPt_ZG_down->Scale(scalefactordown);
	p_JetPt_ZG_down->Draw("same");
	p_JetPt_rare->Draw("same");
	std::cout << std::endl;
	for(int i(1); i<=11; i++){
		std::cout << "down  ratio " << i << " " << p_JetPt_down->GetBinContent(i)/p_JetPt_ZG_down->GetBinContent(i) << std::endl;
		ISRwgt_down[i-1] = p_JetPt_down->GetBinContent(i)/p_JetPt_ZG_down->GetBinContent(i);
	}
	
	p_JetPt_ZG_highEt->Scale(scalefactorhighEt);
	p_JetPt_ZG_highEt->Add(p_JetPt_rare_highEt, 1);
	p_JetPt_ZG_highEt->Sumw2();
	can_JetPt_alter->cd();
	gPad->SetLogy();
	p_JetPt_data_highEt->Draw();
	p_JetPt_ZG_highEt->SetLineColor(kRed);
	p_JetPt_ZG_highEt->Draw("same");
	p_JetPt_rare_highEt->SetLineColor(kGreen);
	p_JetPt_rare_highEt->Draw("same");
	std::cout << std::endl;
	for(int i(1); i<=11; i++){
		std::cout << "alter  ratio " << i << " " << p_JetPt_data_highEt->GetBinContent(i)/p_JetPt_ZG_highEt->GetBinContent(i) << std::endl;
		ISRwgt_alter[i-1] = p_JetPt_data_highEt->GetBinContent(i)/p_JetPt_ZG_highEt->GetBinContent(i);
	}
	
	can_llmass->cd();
	p_llmass_data->Draw();
	p_llmass_up->SetLineColor(kBlack);
	p_llmass_up->Draw("same");
	p_llmass_down->SetLineColor(kYellow);
	p_llmass_down->Draw("same");
	p_llmass_ZG->SetLineColor(kRed);
	p_llmass_ZG->Scale(scalefactor);
	p_llmass_ZG->Draw("same");
	p_llmass_rare->SetLineColor(kGreen);
	p_llmass_rare->Scale(scalefactor);
	p_llmass_rare->Draw("same");

}


