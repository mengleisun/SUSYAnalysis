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


void ISRweight(){//main 

	int channelType = 2;
  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");


	esfScaleFactor  objectESF;

	Double_t plotEtBins[]={35,50,100,150,200,250,300,400,600};
	Double_t plotPtBins[]={25, 35,50,100,150,200,250,300,400,600,800};
	TH1F *p_phoEt_data     = new TH1F("p_phoEt_data","",8,plotEtBins);
	TH1F *p_phoEt_up     	 = new TH1F("p_phoEt_up","",8,plotEtBins);
	TH1F *p_phoEt_down     = new TH1F("p_phoEt_down","",8,plotEtBins);
	TH1F *p_JetPt_data     = new TH1F("p_JetPt_data","",8,plotEtBins);
	TH1F *p_JetPt_up     = new TH1F("p_JetPt_up","",8,plotEtBins);
	TH1F *p_JetPt_down     = new TH1F("p_JetPt_down","",8,plotEtBins);
	TH1F *p_JetPt_fake     = new TH1F("p_JetPt_fake","",8,plotEtBins);
	TH1F *p_JetPt_total    = new TH1F("p_JetPt_total","",8,plotEtBins);
	TH1F *p_llmass_data   = new TH1F("p_llmass_data","",100,30,130);
	TH1F *p_llmass_up   = new TH1F("p_llmass_up","",100,30,130);
	TH1F *p_llmass_down   = new TH1F("p_llmass_down","",100,30,130);
	TH1F *p_llmass_highEt  = new TH1F("p_llmass_highEt","",100,30,130);

	TH1F *p_phoEta_data    = new TH1F("p_phoEta_data","",30,-3,3);
	TH1F *p_phoPhi_data    = new TH1F("p_phoPhi_data","",32,-3.2, 3.2);
	TH1F *p_lepPt_data     = new TH1F("p_lepPt_data","",10, plotPtBins);
	TH1F *p_trailPt_data     = new TH1F("p_trail_data","",8,plotEtBins);
	TH1F *p_lepEta_data    = new TH1F("p_lepEta_data","",30,-3,3);
	TH1F *p_lepPhi_data    = new TH1F("p_lepPhi_data","",32,-3.2, 3.2);
	TH1F *p_sigMET_data    = new TH1F("p_sigMET_data","",20,0,400);
	TH1F *p_sigMT_data     = new TH1F("p_sigMT_data","",20,0,800);
	TH1F *p_sigMETPhi_data = new TH1F("p_sigMETPhi_data","",64,-3.2, 3.2);
	TH1F *p_dPhiLepMET_data= new TH1F("p_dPhiLepMET_data","",64,-3.2, 3.2);
	TH1F *p_nVertex_data   = new TH1F("p_nVertex_data","",100,0,100);
	TH1F *p_dRPhoLep_data  = new TH1F("p_dRPhoLep_data","",60,0.3,3.3);
	TH1F *p_HT_data        = new TH1F("p_HT_data","",100,0,1000);
	TH1F *p_nJet_data      = new TH1F("p_nJet_data","",10,0,10);
	TH1F *p_invmass_data   = new TH1F("p_invmass_data","",120,80,200);

	TH1F *p_phoEt_ZG     = new TH1F("p_phoEt_ZG","",8,plotEtBins);
	TH1F *p_phoEta_ZG    = new TH1F("p_phoEta_ZG","",30,-3,3);
	TH1F *p_phoPhi_ZG    = new TH1F("p_phoPhi_ZG","",32,-3.2, 3.2);
	TH1F *p_lepPt_ZG     = new TH1F("p_lepPt_ZG","",10, plotPtBins);
	TH1F *p_trailPt_ZG     = new TH1F("p_trail_ZG","",8,plotEtBins);
	TH1F *p_lepEta_ZG    = new TH1F("p_lepEta_ZG","",30,-3,3);
	TH1F *p_lepPhi_ZG    = new TH1F("p_lepPhi_ZG","",32,-3.2, 3.2);
	TH1F *p_sigMET_ZG    = new TH1F("p_sigMET_ZG","",20,0,400);
	TH1F *p_sigMT_ZG     = new TH1F("p_sigMT_ZG","",20,0,800);
	TH1F *p_sigMETPhi_ZG = new TH1F("p_sigMETPhi_ZG","",64,-3.2, 3.2);
	TH1F *p_dPhiLepMET_ZG= new TH1F("p_dPhiLepMET_ZG","",64,-3.2, 3.2);
	TH1F *p_nVertex_ZG   = new TH1F("p_nVertex_ZG","",100,0,100);
	TH1F *p_dRPhoLep_ZG  = new TH1F("p_dRPhoLep_ZG","",60,0.3,3.3);
	TH1F *p_HT_ZG        = new TH1F("p_HT_ZG","",100,0,1000);
	TH1F *p_nJet_ZG      = new TH1F("p_nJet_ZG","",10,0,10);
	TH1F *p_invmass_ZG   = new TH1F("p_invmass_ZG","",120,80,200);
	TH1F *p_llmass_ZG   = new TH1F("p_llmass_ZG","",100,30,130);
	TH1F *p_JetPt_ZG     = new TH1F("p_JetPt_ZG","",8,plotEtBins);
	TH1F *p_JetPt_ZG_highEt = new TH1F("p_JetPt_ZG_highEt","",8,plotEtBins);

	TH1F *p_rewgt_phoEt_ZG     = new TH1F("p_rewgt_phoEt_ZG","",8,plotEtBins);
	TH1F *p_rewgt_lepPt_ZG     = new TH1F("p_rewgt_lepPt_ZG","",10, plotPtBins);
	TH1F *p_rewgt_sigMET_ZG    = new TH1F("p_rewgt_sigMET_ZG","",20,0,400);
	TH1F *p_rewgt_sigMT_ZG     = new TH1F("p_rewgt_sigMT_ZG","",20,0,800);
	TH1F *p_rewgt_dPhiLepMET_ZG= new TH1F("p_rewgt_dPhiLepMET_ZG","",64,-3.2, 3.2);
	TH1F *p_rewgt_dRPhoLep_ZG  = new TH1F("p_rewgt_dRPhoLep_ZG","",60,0.3,3.3);

	TH1F *p_phoEt_rare     = new TH1F("p_phoEt_rare","",8,plotEtBins);
	TH1F *p_phoEta_rare    = new TH1F("p_phoEta_rare","",30,-3,3);
	TH1F *p_phoPhi_rare    = new TH1F("p_phoPhi_rare","",32,-3.2, 3.2);
	TH1F *p_lepPt_rare     = new TH1F("p_lepPt_rare","",10, plotPtBins);
	TH1F *p_lepEta_rare    = new TH1F("p_lepEta_rare","",30,-3,3);
	TH1F *p_lepPhi_rare    = new TH1F("p_lepPhi_rare","",32,-3.2, 3.2);
	TH1F *p_sigMET_rare    = new TH1F("p_sigMET_rare","",20,0,400);
	TH1F *p_sigMT_rare     = new TH1F("p_sigMT_rare","",20,0,800);
	TH1F *p_sigMETPhi_rare = new TH1F("p_sigMETPhi_rare","",64,-3.2, 3.2);
	TH1F *p_dPhiLepMET_rare= new TH1F("p_dPhiLepMET_rare","",64,-3.2, 3.2);
	TH1F *p_nVertex_rare   = new TH1F("p_nVertex_rare","",100,0,100);
	TH1F *p_dRPhoLep_rare  = new TH1F("p_dRPhoLep_rare","",60,0.3,3.3);
	TH1F *p_HT_rare        = new TH1F("p_HT_rare","",100,0,1000);
	TH1F *p_nJet_rare      = new TH1F("p_nJet_rare","",10,0,10);
	TH1F *p_JetPt_rare     = new TH1F("p_JetPt_rare","",8,plotEtBins);
	TH1F *p_invmass_rare   = new TH1F("p_invmass_rare","",120,80,200);
	TH1F *p_llmass_rare   = new TH1F("p_llmass_rare","",100,30,130);

	TH1F *p_phoEt_DY     = new TH1F("p_phoEt_DY","",8,plotEtBins);
	TH1F *p_phoEta_DY    = new TH1F("p_phoEta_DY","",30,-3,3);
	TH1F *p_phoPhi_DY    = new TH1F("p_phoPhi_DY","",32,-3.2, 3.2);
	TH1F *p_lepPt_DY     = new TH1F("p_lepPt_DY","",10, plotPtBins);
	TH1F *p_lepEta_DY    = new TH1F("p_lepEta_DY","",30,-3,3);
	TH1F *p_lepPhi_DY    = new TH1F("p_lepPhi_DY","",32,-3.2, 3.2);
	TH1F *p_sigMET_DY    = new TH1F("p_sigMET_DY","",20,0,400);
	TH1F *p_sigMT_DY     = new TH1F("p_sigMT_DY","",20,0,800);
	TH1F *p_sigMETPhi_DY = new TH1F("p_sigMETPhi_DY","",64,-3.2, 3.2);
	TH1F *p_dPhiLepMET_DY= new TH1F("p_dPhiLepMET_DY","",64,-3.2, 3.2);
	TH1F *p_nVertex_DY   = new TH1F("p_nVertex_DY","",100,0,100);
	TH1F *p_dRPhoLep_DY  = new TH1F("p_dRPhoLep_DY","",60,0.3,3.3);
	TH1F *p_HT_DY        = new TH1F("p_HT_DY","",100,0,1000);
	TH1F *p_nJet_DY      = new TH1F("p_nJet_DY","",10,0,10);
	TH1F *p_JetPt_DY     = new TH1F("p_JetPt_DY","",8,plotEtBins);
	TH1F *p_invmass_DY   = new TH1F("p_invmass_DY","",120,80,200);
	TH1F *p_llmass_DY   = new TH1F("p_llmass_DY","",100,30,130);

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

		if(phoEt > 599)phoEt = 599;
		if(JetPt > 599)JetPt = 599;
		if(dRPhoLep < 0.8)continue;
		if(dilepMass < 80 || dilepMass > 100)continue;
		//if(dilepMass < 30 || dilepMass > 80 )continue;
		double weight = 1;
		double weightup = 1;
		double weightdown = 1; 
		if(dilepMass > 80 && dilepMass < 100){
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

		p_JetPt_data->Fill(JetPt, weight);
	 p_JetPt_fake->Fill(JetPt, 1-weight);
    	
		p_phoEt_data->Fill(phoEt, weight); 
		p_llmass_data->Fill(dilepMass, weight);
		p_JetPt_up->Fill(JetPt, weightup);	
		p_phoEt_up->Fill(phoEt, weightup); 
		p_llmass_up->Fill(dilepMass, weightup);
		p_JetPt_down->Fill(JetPt, weightdown);	
		p_phoEt_down->Fill(phoEt, weightdown); 
		p_llmass_down->Fill(dilepMass, weightdown);

		if(phoEt > 145)p_llmass_highEt->Fill(dilepMass, weight);

		p_phoEta_data->Fill(phoEta, weight);
		p_phoPhi_data->Fill(phoPhi, weight);
		p_lepPt_data->Fill(lepPt, weight);
		p_lepEta_data->Fill(lepEta, weight);
		p_lepPhi_data->Fill(lepPhi, weight);
		p_trailPt_data->Fill(trailPt, weight);
		p_sigMET_data->Fill(sigMET, weight);
		p_sigMT_data->Fill(sigMT, weight);
		p_sigMETPhi_data->Fill(sigMETPhi, weight);
		p_dPhiLepMET_data->Fill(dPhiLepMET, weight);
		p_nVertex_data->Fill(nVertex, weight);
		p_dRPhoLep_data->Fill(dRPhoLep, weight);
		p_HT_data->Fill(HT, weight);
		p_nJet_data->Fill(nJet, weight);
		p_invmass_data->Fill(threeMass, weight);
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

		double weight = 35.8*1000*117.8/14372399;

		double scalefactor(0);
		if(channelType == 1){
			scalefactor = objectESF.getElectronESF(ZG_lepPt,ZG_lepEta)*objectESF.getPhotonESF(ZG_phoEt,ZG_phoEta)*objectESF.getegPhotonTRGESF(ZG_phoEt,ZG_phoEta)*objectESF.getElectronTRGESF(ZG_lepPt,ZG_lepEta);
		}
		if(channelType == 2){
			scalefactor = objectESF.getMuonESF(ZG_lepPt,ZG_lepEta)*objectESF.getPhotonESF(ZG_phoEt,ZG_phoEta)*objectESF.getMuonEGTRGESF(ZG_phoEt, ZG_lepPt);
		}
		p_scalefactor->Fill(ZG_lepEta,objectESF.getMuonESF(ZG_lepPt,ZG_lepEta));

		if(ZG_phoEt > 599)ZG_phoEt = 599;
		if(ZG_JetPt > 599)ZG_JetPt = 599;
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

		totalevent += 1;
		double reweightF = 1;
		if(ZG_bosonPt < 50)reweightF = 1;
		else if(ZG_bosonPt >= 50 && ZG_bosonPt < 100)reweightF  = 1.09297;
		else if(ZG_bosonPt >= 100 && ZG_bosonPt < 150)reweightF = 0.833083;
		else if(ZG_bosonPt >= 150 && ZG_bosonPt < 200)reweightF = 0.713514;
		else if(ZG_bosonPt >= 200 && ZG_bosonPt < 250)reweightF =	0.744641;
		else if(ZG_bosonPt >= 250 && ZG_bosonPt < 300)reweightF =	0.737551;
		else if(ZG_bosonPt >= 300 && ZG_bosonPt < 400)reweightF = 0.697719;
		else if(ZG_bosonPt >= 400)reweightF =	0.36516;

		reweightevent += reweightF;
	
		p_rewgt_phoEt_ZG->Fill(ZG_phoEt, reweightF*weight);
		p_rewgt_lepPt_ZG->Fill(ZG_lepPt, reweightF*weight);
		p_rewgt_sigMET_ZG->Fill(ZG_sigMET, reweightF*weight);
		p_rewgt_sigMT_ZG->Fill(ZG_sigMT, reweightF*weight);
		p_rewgt_dPhiLepMET_ZG->Fill(ZG_dPhiLepMET, reweightF*weight);
		p_rewgt_dRPhoLep_ZG->Fill(ZG_dRPhoLep, reweightF*weight);

		if(ZG_phoEt > 145)p_JetPt_ZG_highEt->Fill(ZG_JetPt, weight);
		p_JetPt_ZG->Fill(ZG_JetPt, weight);	
		p_phoEt_ZG->Fill(ZG_phoEt, weight); 
		p_phoEta_ZG->Fill(ZG_phoEta, weight);
		p_phoPhi_ZG->Fill(ZG_phoPhi, weight);
		p_lepPt_ZG->Fill(ZG_lepPt, weight);
		p_lepEta_ZG->Fill(ZG_lepEta, weight);
		p_lepPhi_ZG->Fill(ZG_lepPhi, weight);
		p_trailPt_ZG->Fill(ZG_trailPt, weight);
		p_sigMET_ZG->Fill(ZG_sigMET, weight);
		p_sigMT_ZG->Fill(ZG_sigMT, weight);
		p_sigMETPhi_ZG->Fill(ZG_sigMETPhi, weight);
		p_dPhiLepMET_ZG->Fill(ZG_dPhiLepMET, weight);
		p_nVertex_ZG->Fill(ZG_nVertex, weight);
		p_dRPhoLep_ZG->Fill(ZG_dRPhoLep, weight);
		p_HT_ZG->Fill(ZG_HT, weight);
		p_nJet_ZG->Fill(ZG_nJet, weight);
		p_invmass_ZG->Fill(ZG_threeMass, weight);
		p_llmass_ZG->Fill(ZG_dilepMass, weight);
	}//loop on  events

	p_rewgt_phoEt_ZG->Scale(totalevent/reweightevent);
	p_rewgt_lepPt_ZG->Scale(totalevent/reweightevent);
	p_rewgt_sigMET_ZG->Scale(totalevent/reweightevent);
	p_rewgt_sigMT_ZG->Scale(totalevent/reweightevent);
	p_rewgt_dPhiLepMET_ZG->Scale(totalevent/reweightevent);
	p_rewgt_dRPhoLep_ZG->Scale(totalevent/reweightevent);
	

//************ Signal Tree **********************//
  TChain *raretree = new TChain("ZTree");
  raretree->Add("/uscms_data/d3/mengleis/test/resTree_ISR_TTG.root");
  raretree->Add("/uscms_data/d3/mengleis/test/resTree_ISR_WWG.root");
  raretree->Add("/uscms_data/d3/mengleis/test/resTree_ISR_WZG.root");
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

		double weight = MCweight; 
		double scalefactor(0);
		if(channelType == 1){
			scalefactor = objectESF.getElectronESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getegPhotonTRGESF(phoEt,phoEta)*objectESF.getElectronTRGESF(lepPt,lepEta);
		}
		if(channelType == 2){
			scalefactor = objectESF.getMuonESF(rare_lepPt,rare_lepEta)*objectESF.getPhotonESF(rare_phoEt,rare_phoEta)*objectESF.getMuonEGTRGESF(rare_phoEt, rare_lepPt);
		}

		if(rare_phoEt > 599)rare_phoEt = 599;
		if(rare_JetPt > 599)rare_JetPt = 599;
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
		p_phoEta_rare->Fill(rare_phoEta, weight);
		p_phoPhi_rare->Fill(rare_phoPhi, weight);
		p_lepPt_rare->Fill(rare_lepPt, weight);
		p_lepEta_rare->Fill(rare_lepEta, weight);
		p_lepPhi_rare->Fill(rare_lepPhi, weight);
		p_sigMET_rare->Fill(rare_sigMET, weight);
		p_sigMT_rare->Fill(rare_sigMT, weight);
		p_sigMETPhi_rare->Fill(rare_sigMETPhi, weight);
		p_dPhiLepMET_rare->Fill(rare_dPhiLepMET, weight);
		p_nVertex_rare->Fill(rare_nVertex, weight);
		p_dRPhoLep_rare->Fill(rare_dRPhoLep, weight);
		p_HT_rare->Fill(rare_HT, weight);
		p_nJet_rare->Fill(rare_nJet, weight);
		p_JetPt_rare->Fill(rare_JetPt, weight);
		p_invmass_rare->Fill(rare_threeMass, weight);
		p_llmass_rare->Fill(rare_dilepMass, weight);
	}//loop on  events
	
  TChain *DYtree = new TChain("ZTree");
  DYtree->Add("/uscms_data/d3/mengleis/test/resTree_ISR_ZGNLO130.root");
	float DYweight(0);
  float DY_phoEt(0);
  float DY_phoEta(0);
  float DY_phoPhi(0);
  float DY_lepPt(0);
  float DY_lepEta(0);
  float DY_lepPhi(0);
  float DY_sigMT(0);
  float DY_sigMET(0);
  float DY_sigMETPhi(0);
  float DY_dPhiLepMET(0);
  int   DY_nVertex(0);
  float DY_dRPhoLep(0);
  float DY_HT(0);
  float DY_nJet(0);
	float DY_threeMass(0);
	float DY_dilepMass(0);
	float DY_JetPt(0);
  std::vector<int>   *DY_mcPID=0;
  std::vector<float> *DY_mcEta=0;
  std::vector<float> *DY_mcPhi=0;
  std::vector<float> *DY_mcPt=0;
  std::vector<int>   *DY_mcMomPID=0;
  std::vector<int>   *DY_mcGMomPID=0;

  DYtree->SetBranchAddress("MCweight",  &DYweight);
  DYtree->SetBranchAddress("phoEt",     &DY_phoEt);
  DYtree->SetBranchAddress("phoEta",    &DY_phoEta);
  DYtree->SetBranchAddress("phoPhi",    &DY_phoPhi);
  DYtree->SetBranchAddress("lepPt",     &DY_lepPt);
  DYtree->SetBranchAddress("lepEta",    &DY_lepEta);
  DYtree->SetBranchAddress("lepPhi",    &DY_lepPhi);
  DYtree->SetBranchAddress("sigMT",     &DY_sigMT);
  DYtree->SetBranchAddress("sigMET",    &DY_sigMET);
  DYtree->SetBranchAddress("sigMETPhi", &DY_sigMETPhi);
  DYtree->SetBranchAddress("dPhiLepMET",&DY_dPhiLepMET);
  DYtree->SetBranchAddress("nVertex",   &DY_nVertex);
  DYtree->SetBranchAddress("dRPhoLep",  &DY_dRPhoLep);
  DYtree->SetBranchAddress("HT",        &DY_HT);
  DYtree->SetBranchAddress("nJet",      &DY_nJet);
	DYtree->SetBranchAddress("JetPt",     &DY_JetPt);
	DYtree->SetBranchAddress("threeMass", &DY_threeMass);
	DYtree->SetBranchAddress("dilepMass", &DY_dilepMass);
  DYtree->SetBranchAddress("mcPID",    &DY_mcPID);
  DYtree->SetBranchAddress("mcEta",    &DY_mcEta);
  DYtree->SetBranchAddress("mcPhi",    &DY_mcPhi);
  DYtree->SetBranchAddress("mcPt",     &DY_mcPt);
  DYtree->SetBranchAddress("mcMomPID", &DY_mcMomPID);
  DYtree->SetBranchAddress("mcGMomPID",&DY_mcGMomPID);

  for(unsigned ievt(0); ievt<DYtree->GetEntries(); ++ievt){//loop on entries
		DYtree->GetEntry(ievt);

		double weight = DYweight;	
		double scalefactor(0);
		if(channelType == 1){
			scalefactor = objectESF.getElectronESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getegPhotonTRGESF(phoEt,phoEta)*objectESF.getElectronTRGESF(lepPt,lepEta);
		}
		if(channelType == 2){
			scalefactor = objectESF.getMuonESF(DY_lepPt,DY_lepEta)*objectESF.getPhotonESF(DY_phoEt,DY_phoEta)*objectESF.getMuonEGTRGESF(DY_phoEt, DY_lepPt);
		}

		if(DY_phoEt > 599)DY_phoEt = 599;
		if(DY_JetPt > 599)DY_JetPt = 599;
		if(DY_dRPhoLep < 0.8)continue;
		if(DY_phoEt < 145)continue;
		if(DY_dilepMass < 80 || DY_dilepMass > 100)continue;
		weight = weight*scalefactor;
		float PUweight = getPUESF(DY_nVertex);
		weight = weight*PUweight;

	//	bool   isTruePho(false);
	//	double mindR(0.3);
	//	unsigned phoIndex(0);
	//	for(unsigned iMC(0); iMC< DY_mcPID->size(); iMC++){
	//		double dR1 = DeltaR((*DY_mcEta)[iMC], (*DY_mcPhi)[iMC], DY_phoEta, DY_phoPhi);
	//		double dE1 = fabs((*DY_mcPt)[iMC] - DY_phoEt)/DY_phoEt;
	//		if(dR1 < mindR && dE1 < 0.2){mindR=dR1; phoIndex=iMC;}
	//	}
	//	if(mindR < 0.1){
	//		if((*DY_mcPID)[phoIndex] == 22 && (fabs((*DY_mcMomPID)[phoIndex])<25 || fabs((*DY_mcMomPID)[phoIndex])==999))isTruePho=true;
	//	}
	//	if(isTruePho)continue;


		p_phoEt_DY->Fill(DY_phoEt, weight); 
		p_phoEta_DY->Fill(DY_phoEta, weight);
		p_phoPhi_DY->Fill(DY_phoPhi, weight);
		p_lepPt_DY->Fill(DY_lepPt, weight);
		p_lepEta_DY->Fill(DY_lepEta, weight);
		p_lepPhi_DY->Fill(DY_lepPhi, weight);
		p_sigMET_DY->Fill(DY_sigMET, weight);
		p_sigMT_DY->Fill(DY_sigMT, weight);
		p_sigMETPhi_DY->Fill(DY_sigMETPhi, weight);
		p_dPhiLepMET_DY->Fill(DY_dPhiLepMET, weight);
		p_nVertex_DY->Fill(DY_nVertex, weight);
		p_dRPhoLep_DY->Fill(DY_dRPhoLep, weight);
		p_HT_DY->Fill(DY_HT, weight);
		p_nJet_DY->Fill(DY_nJet, weight);
		p_JetPt_DY->Fill(DY_JetPt, weight);
		p_invmass_DY->Fill(DY_threeMass, weight);
		p_llmass_DY->Fill(DY_dilepMass, weight);
	}//loop on  events

	TCanvas *can_phoEt     = new TCanvas("can_phoEt",       "can_phoEt", 600,600); 
//	TCanvas *can_phoEta    = new TCanvas("can_phoEta",      "can_phoEta",600,600); 
//	TCanvas *can_phoPhi    = new TCanvas("can_phoPhi",      "can_phoPhi",600,600); 
//	TCanvas *can_lepPt     = new TCanvas("can_lepPt",       "can_lepPt",600,600); 
//	TCanvas *can_lepEta    = new TCanvas("can_lepEta",      "can_lepEta",600,600); 
//	TCanvas *can_lepPhi    = new TCanvas("can_lepPhi",      "can_lepPhi",600,600); 
//	TCanvas *can_sigMET    = new TCanvas("can_sigMET",      "can_sigMET",600,600); 
//	TCanvas *can_sigMT     = new TCanvas("can_sigMT",       "can_sigMT",600,600); 
//	TCanvas *can_sigMETPhi = new TCanvas("can_sigMETPhi",   "can_sigMETPhi",600,600); 
//	TCanvas *can_dPhiLepMET= new TCanvas("can_dPhiLepMET",  "can_dPhiLepMET",600,600); 
//	TCanvas *can_nVertex   = new TCanvas("can_nVertex",     "can_nVertex",600,600); 
//	TCanvas *can_dRPhoLep  = new TCanvas("can_dRPhoLep",    "can_dRPhoLep",600,600); 
	TCanvas *can_HT        = new TCanvas("can_HT",          "can_HT",600,600); 
	TCanvas *can_HT_up     = new TCanvas("can_HT_up",       "can_HT_up",600,600); 
	TCanvas *can_HT_down   = new TCanvas("can_HT_down",     "can_HT_down",600,600); 
	TCanvas *can_HT_alter  = new TCanvas("can_HT_alter",    "can_HT_alter",600,600); 
//	TCanvas *can_nJet      = new TCanvas("can_nJet",        "can_nJet",600,600); 
//	TCanvas *can_invmass   = new TCanvas("can_invmass",     "can_invmass", 600,600);    
	TCanvas *can_llmass   = new TCanvas("can_llmass",     "can_llmass", 600,600);   


	p_llmass_ZG->Add(p_llmass_rare);
	float scalefactor = p_llmass_data->Integral(50,70)/p_llmass_ZG->Integral(50,70); 
	std::cout << "scale = " << scalefactor<<std::endl;
	float scalefactorup = p_llmass_up->Integral(50,70)/p_llmass_ZG->Integral(50,70); 
	std::cout << "scale = " << scalefactorup<<std::endl;
	float scalefactordown = p_llmass_down->Integral(50,70)/p_llmass_ZG->Integral(50,70); 
	std::cout << "scale = " << scalefactordown<<std::endl;

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
	p_rewgt_phoEt_ZG->Add(p_phoEt_rare);
	p_rewgt_phoEt_ZG->SetLineColor(kMagenta);
	p_rewgt_phoEt_ZG->Scale(scalefactor);
	p_rewgt_phoEt_ZG->Draw("same");
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
	

//	can_phoEta->cd();
//	p_phoEta_data->Draw();
//	p_phoEta_ZG->Add(p_phoEta_rare);
//	p_phoEta_ZG->SetLineColor(kRed);
//	p_phoEta_ZG->Scale(scalefactor);
//	p_phoEta_ZG->Draw("same");
//	p_phoEta_rare->SetLineColor(kGreen);
//	p_phoEta_rare->Scale(scalefactor);
//	p_phoEta_rare->Draw("same");
//
//	can_phoPhi->cd();
//	p_phoPhi_data->Draw();
//	p_phoPhi_ZG->Add(p_phoPhi_rare);
//	p_phoPhi_ZG->SetLineColor(kRed);
//	p_phoPhi_ZG->Scale(scalefactor);
//	p_phoPhi_ZG->Draw("same");
//	p_phoPhi_rare->SetLineColor(kGreen);
//	p_phoPhi_rare->Scale(scalefactor);
//	p_phoPhi_rare->Draw("same");
//
//	can_lepPt->cd();
//	gPad->SetLogy();
//	p_lepPt_data->Draw();
//	p_lepPt_ZG->Add(p_lepPt_rare);
//	p_lepPt_ZG->SetLineColor(kRed);
//	p_lepPt_ZG->Scale(scalefactor);
//	p_lepPt_ZG->Draw("same");
//	p_lepPt_rare->SetLineColor(kGreen);
//	p_lepPt_rare->Scale(scalefactor);
//	p_lepPt_rare->Draw("same");
//	p_rewgt_lepPt_ZG->Add(p_lepPt_rare);
//	p_rewgt_lepPt_ZG->SetLineColor(kMagenta);
//	p_rewgt_lepPt_ZG->Scale(scalefactor);
//	p_rewgt_lepPt_ZG->Draw("same");
//
//	can_lepEta->cd();
//	p_lepEta_data->Draw();
//	p_lepEta_ZG->Add(p_lepEta_rare);
//	p_lepEta_ZG->SetLineColor(kRed);
//	p_lepEta_ZG->Scale(scalefactor);
//	p_lepEta_ZG->Draw("same");
//	p_lepEta_rare->SetLineColor(kGreen);
//	p_lepEta_rare->Scale(scalefactor);
//	p_lepEta_rare->Draw("same");
//
//	can_lepPhi->cd();
//	p_lepPhi_data->Draw();
//	p_lepPhi_ZG->Add(p_lepPhi_rare);
//	p_lepPhi_ZG->SetLineColor(kRed);
//	p_lepPhi_ZG->Scale(scalefactor);
//	p_lepPhi_ZG->Draw("same");
//	p_lepPhi_rare->SetLineColor(kGreen);
//	p_lepPhi_rare->Scale(scalefactor);
//	p_lepPhi_rare->Draw("same");
//
//	can_sigMET->cd();
//	p_sigMET_data->Draw();
//	p_sigMET_ZG->Add(p_sigMET_rare);
//	p_sigMET_ZG->SetLineColor(kRed);
//	p_sigMET_ZG->Scale(scalefactor);
//	p_sigMET_ZG->Draw("same");
//	p_sigMET_rare->SetLineColor(kGreen);
//	p_sigMET_rare->Scale(scalefactor);
//	p_sigMET_rare->Draw("same");
//	p_rewgt_sigMET_ZG->Add(p_sigMET_rare);
//	p_rewgt_sigMET_ZG->SetLineColor(kMagenta);
//	p_rewgt_sigMET_ZG->Scale(scalefactor);
//	p_rewgt_sigMET_ZG->Draw("same");
//
//	can_sigMT->cd();
//	p_sigMT_data->Draw();
//	p_sigMT_ZG->Add(p_sigMT_rare);
//	p_sigMT_ZG->SetLineColor(kRed);
//	p_sigMT_ZG->Scale(scalefactor);
//	p_sigMT_ZG->Draw("same");
//	p_sigMT_rare->SetLineColor(kGreen);
//	p_sigMT_rare->Scale(scalefactor);
//	p_sigMT_rare->Draw("same");
//	p_rewgt_sigMT_ZG->Add(p_sigMT_rare);
//	p_rewgt_sigMT_ZG->SetLineColor(kMagenta);
//	p_rewgt_sigMT_ZG->Scale(scalefactor);
//	p_rewgt_sigMT_ZG->Draw("same");
//
//	can_sigMETPhi->cd();
//	p_sigMETPhi_data->Draw();
//	p_sigMETPhi_ZG->Add(p_sigMETPhi_rare);
//	p_sigMETPhi_ZG->SetLineColor(kRed);
//	p_sigMETPhi_ZG->Scale(scalefactor);
//	p_sigMETPhi_ZG->Draw("same");
//	p_sigMETPhi_rare->SetLineColor(kGreen);
//	p_sigMETPhi_rare->Scale(scalefactor);
//	p_sigMETPhi_rare->Draw("same");
//
//	can_dPhiLepMET->cd();
//	p_dPhiLepMET_data->Draw();
//	p_dPhiLepMET_ZG->Add(p_dPhiLepMET_rare);
//	p_dPhiLepMET_ZG->SetLineColor(kRed);
//	p_dPhiLepMET_ZG->Scale(scalefactor);
//	p_dPhiLepMET_ZG->Draw("same");
//	p_dPhiLepMET_rare->SetLineColor(kGreen);
//	p_dPhiLepMET_rare->Scale(scalefactor);
//	p_dPhiLepMET_rare->Draw("same");
//	p_rewgt_dPhiLepMET_ZG->Add(p_dPhiLepMET_rare);
//	p_rewgt_dPhiLepMET_ZG->SetLineColor(kMagenta);
//	p_rewgt_dPhiLepMET_ZG->Scale(scalefactor);
//	p_rewgt_dPhiLepMET_ZG->Draw("same");
//
//	can_nVertex->cd();
//	p_nVertex_data->Draw();
//	p_nVertex_ZG->Add(p_nVertex_rare);
//	p_nVertex_ZG->SetLineColor(kRed);
//	p_nVertex_ZG->Scale(scalefactor);
//	p_nVertex_ZG->Draw("same");
//	p_nVertex_rare->SetLineColor(kGreen);
//	p_nVertex_rare->Scale(scalefactor);
//	p_nVertex_rare->Draw("same");
//
//	can_dRPhoLep->cd();
//	p_dRPhoLep_data->Draw();
//	p_dRPhoLep_ZG->Add(p_dRPhoLep_rare);
//	p_dRPhoLep_ZG->SetLineColor(kRed);
//	p_dRPhoLep_ZG->Scale(scalefactor);
//	p_dRPhoLep_ZG->Draw("same");
//	p_dRPhoLep_rare->SetLineColor(kGreen);
//	p_dRPhoLep_rare->Scale(scalefactor);
//	p_dRPhoLep_rare->Draw("same");
//	p_rewgt_dRPhoLep_ZG->Add(p_dRPhoLep_rare);
//	p_rewgt_dRPhoLep_ZG->SetLineColor(kMagenta);
//	p_rewgt_dRPhoLep_ZG->Scale(scalefactor);
//	p_rewgt_dRPhoLep_ZG->Draw("same");
//
	double ISRwgt_norm[8];
	double ISRwgt_up[8];
	double ISRwgt_down[8];
	double ISRwgt_alter[8];

	can_HT->cd();
	gPad->SetLogy();
	gStyle->SetOptStat(0);
	p_JetPt_data->GetXaxis()->SetTitle("ISR P_{T} (GeV)");
	p_JetPt_data->SetMinimum(1);
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
	can_HT->SaveAs("ISRweight.pdf");
	for(int i(1); i<9; i++){
		std::cout << "norm ratio " << i << " " << p_JetPt_data->GetBinContent(i)/p_JetPt_ZG->GetBinContent(i) << std::endl;
		ISRwgt_norm[i-1] = p_JetPt_data->GetBinContent(i)/p_JetPt_ZG->GetBinContent(i);
	}

	can_HT_up->cd();
	gPad->SetLogy();
	p_JetPt_up->Draw();
	p_JetPt_ZG->Scale(scalefactorup/scalefactor);
	p_JetPt_ZG->Draw("same");
	p_JetPt_rare->Draw("same");
	std::cout << std::endl;
	for(int i(1); i<9; i++){
		std::cout << "up   ratio " << i << " " << p_JetPt_up->GetBinContent(i)/p_JetPt_ZG->GetBinContent(i) << std::endl;
		ISRwgt_up[i-1] = p_JetPt_up->GetBinContent(i)/p_JetPt_ZG->GetBinContent(i);	
	}

	can_HT_down->cd();
	gPad->SetLogy();
	p_JetPt_down->Draw();
	p_JetPt_ZG->Scale(scalefactordown/scalefactorup);
	p_JetPt_ZG->Draw("same");
	p_JetPt_rare->Draw("same");
	std::cout << std::endl;
	for(int i(1); i<9; i++){
		std::cout << "down  ratio " << i << " " << p_JetPt_down->GetBinContent(i)/p_JetPt_ZG->GetBinContent(i) << std::endl;
		ISRwgt_down[i-1] = p_JetPt_down->GetBinContent(i)/p_JetPt_ZG->GetBinContent(i);
	}
	
	p_JetPt_ZG->Scale(scalefactor/scalefactordown);
	double scaleDY = p_llmass_highEt->Integral(50,70)/p_llmass_DY->Integral(50,70);
	p_JetPt_DY->Scale(scaleDY);
	p_JetPt_ZG->Add( p_JetPt_ZG_highEt, -1);
	p_JetPt_ZG->Add( p_JetPt_DY, 1);
	p_JetPt_ZG->Sumw2();
	can_HT_alter->cd();
	gPad->SetLogy();
	p_JetPt_data->Draw();
	p_JetPt_ZG->Draw("same");
	p_JetPt_rare->Draw("same");
	p_JetPt_DY->SetLineColor(kCyan);	
	p_JetPt_DY->Draw("same");
	std::cout << std::endl;
	for(int i(1); i<9; i++){
		std::cout << "alter  ratio " << i << " " << p_JetPt_data->GetBinContent(i)/p_JetPt_ZG->GetBinContent(i) << std::endl;
		ISRwgt_alter[i-1] = p_JetPt_data->GetBinContent(i)/p_JetPt_ZG->GetBinContent(i);
	}

	std::cout << std::endl;
	for(int i(0); i<8; i++){
		double err1 = max(fabs(ISRwgt_up[i] - ISRwgt_norm[i]), fabs(ISRwgt_down[i] - ISRwgt_norm[i]));
		//double err2 = fabs(ISRwgt_alter[i] - ISRwgt_norm[i]);
		double err2 = 0; 
		std::cout << ISRwgt_norm[i] << " " << (err1>err2? err1: err2) << std::endl;
	} 
	
//	can_nJet->cd();
//	p_nJet_data->Draw();
//	p_nJet_ZG->Add(p_nJet_rare);
//	p_nJet_ZG->SetLineColor(kRed);
//	p_nJet_ZG->Scale(scalefactor);
//	p_nJet_ZG->Draw("same");
//	p_nJet_rare->SetLineColor(kGreen);
//	p_nJet_rare->Scale(scalefactor);
//	p_nJet_rare->Draw("same");
//
//	can_invmass->cd();
//	p_invmass_data->Draw();
//	p_invmass_ZG->Add(p_invmass_rare);
//	p_invmass_ZG->SetLineColor(kRed);
//	p_invmass_ZG->Scale(scalefactor);
//	p_invmass_ZG->Draw("same");
//	p_invmass_rare->SetLineColor(kGreen);
//	p_invmass_rare->Scale(scalefactor);
//	p_invmass_rare->Draw("same");
//
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


