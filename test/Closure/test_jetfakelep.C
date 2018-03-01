#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TF3.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TChain.h"
#include "TSystem.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TGraphErrors.h"

#include "../../include/analysis_rawData.h"
#include "../../include/analysis_photon.h"
#include "../../include/analysis_ele.h"
#include "../../include/analysis_muon.h"
#include "../../include/analysis_mcData.h"
#include "../../include/analysis_tools.h"
#include "../../include/analysis_fakes.h"
#include "../../include/tdrstyle.C"

#define NTOY 1000


void test_jetfakelep(int ichannel){

	setTDRStyle();   
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFillColor(0);
  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");
  int channelType = ichannel; // eg = 1; mg =2;

	//MC Fake Tree //
	//*********** hist o list **********************//
	TH1D *mcpred_PhoEt = new TH1D("mcpred_PhoEt","#gamma E_{T}; E_{T} (GeV)",20,0,200);
	TH1D *mcpred_LepPt = new TH1D("mcpred_LepPt","mcpred_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *mcpred_MET = new TH1D("mcpred_MET","MET; MET (GeV);",20,0,100);
	TH1D *mcpred_Mt = new TH1D("mcpred_Mt","M_{T}; M_{T} (GeV);",40,0,200);
	TH1D *mcpred_HT = new TH1D("mcpred_HT","HT; HT (GeV);",40,0,400); 
	TH1D *mcpred_PhoEta = new TH1D("mcpred_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *mcpred_LepEta = new TH1D("mcpred_LepEta","mcpred_LepEta",60,-3,3);
	TH1D *mcpred_dPhiEleMET = new TH1D("mcpred_dPhiEleMET","dPhiEleMET",32,0,3.2); 
	TH1D *mcpred_nJet = new TH1D("mcpred_nJet","mcpred_nJet",10,0,10);

	//************ Proxy Tree **********************//
	TChain *mcproxytree = new TChain("signalTree");
	//if(channelType==1)mcproxytree->Add("/uscms_data/d3/mengleis/FullStatusOct/fakelep_egsignal_QCD.root");
	//if(channelType==1)mcproxytree->Add("/uscms_data/d3/mengleis/FullStatusOct/fakelep_egsignal_GJet.root");
	if(channelType==1)mcproxytree->Add("/uscms_data/d3/mengleis/FullStatusOct/fakelep_QCD.root");
	if(channelType==2)mcproxytree->Add("/uscms_data/d3/mengleis/FullStatusOct/fakelep_mgsignal_QCD.root");

	float mcproxyphoEt(0);
	float mcproxyphoEta(0);
	float mcproxyphoPhi(0);
	float mcproxylepPt(0);
	float mcproxylepEta(0);
	float mcproxylepPhi(0);
	float mcproxysigMT(0);
	float mcproxysigMET(0);
	float mcproxysigMETPhi(0);
	float mcproxydPhiLepMET(0);
	int   mcproxynVertex(0);
	float mcproxydRPhoLep(0);
	float mcproxyHT(0);
	float mcproxynJet(0);
	float mcfakeLepMiniIso(0);
  float fakeLepSigma(0);
  float fakeLepdEta(0);
  float fakeLepdPhi(0);
	std::vector<float> *jetPt=0;
	std::vector<float> *jetEta=0;
	std::vector<float> *jetPhi=0;
  std::vector<int>   *mcPID=0;
  std::vector<float> *mcEta=0;
  std::vector<float> *mcPhi=0;
  std::vector<float> *mcPt=0;
  std::vector<int>   *mcMomPID=0;
  std::vector<int>   *mcGMomPID=0;
	
	mcproxytree->SetBranchAddress("phoEt",     	 &mcproxyphoEt);
	mcproxytree->SetBranchAddress("phoEta",    	 &mcproxyphoEta);
	mcproxytree->SetBranchAddress("phoPhi",    	 &mcproxyphoPhi);
	mcproxytree->SetBranchAddress("lepPt",     	 &mcproxylepPt);
	mcproxytree->SetBranchAddress("lepEta",    	 &mcproxylepEta);
	mcproxytree->SetBranchAddress("lepPhi",    	 &mcproxylepPhi);
	mcproxytree->SetBranchAddress("sigMT",     	 &mcproxysigMT);
	mcproxytree->SetBranchAddress("sigMET",    	 &mcproxysigMET);
	mcproxytree->SetBranchAddress("sigMETPhi", 	 &mcproxysigMETPhi);
	mcproxytree->SetBranchAddress("dPhiLepMET",	 &mcproxydPhiLepMET);
//  mcproxytree->SetBranchAddress("fakeLepMiniIso", &mcfakeLepMiniIso);
//	mcproxytree->SetBranchAddress("nVertex",   	 &mcproxynVertex);
//	mcproxytree->SetBranchAddress("dRPhoLep",  	 &mcproxydRPhoLep);
//	mcproxytree->SetBranchAddress("HT",        	 &mcproxyHT);
//	mcproxytree->SetBranchAddress("nJet",      	 &mcproxynJet);
//	mcproxytree->SetBranchAddress("fakeLepSigma",&fakeLepSigma);
//	mcproxytree->SetBranchAddress("fakeLepdEta", &fakeLepdEta);
//	mcproxytree->SetBranchAddress("fakeLepdPhi", &fakeLepdPhi);
//	mcproxytree->SetBranchAddress("JetPt",     &jetPt);
//	mcproxytree->SetBranchAddress("JetEta",    &jetEta);
//	mcproxytree->SetBranchAddress("JetPhi",    &jetPhi);
  mcproxytree->SetBranchAddress("mcPID",     &mcPID);
  mcproxytree->SetBranchAddress("mcEta",     &mcEta);
  mcproxytree->SetBranchAddress("mcPhi",     &mcPhi);
  mcproxytree->SetBranchAddress("mcPt",      &mcPt);
  mcproxytree->SetBranchAddress("mcMomPID",  &mcMomPID);
  mcproxytree->SetBranchAddress("mcGMomPID", &mcGMomPID);

	for (unsigned ievt(0); ievt<mcproxytree->GetEntries(); ++ievt){//loop on entries
		mcproxytree->GetEntry(ievt);

		double weight = 1; 
		
		/** cut flow *****/
//		if(mcproxyphoEt < 40 || mcproxylepPt < 25)continue;
		if(fabs(mcproxyphoEta) > 1.4442 || fabs(mcproxylepEta) > 2.5)continue;


//		bool isProxy(false);
//		if(channelType==1){if(mcfakeLepMiniIso < 0.4)isProxy=true;}
//		else if(channelType==2){if((mcfakeLepMiniIso > 0.2 && mcfakeLepMiniIso < 0.4))isProxy=true;}
//		if(!isProxy)continue;

//		bool passFakeProxy(true);
//	  if( fabs(mcproxylepEta) < 1.4442){
//  		if( fakeLepSigma > 0.00998 || fabs(fakeLepdEta) > 0.00311 || fabs(fakeLepdPhi) > 0.103)passFakeProxy=false; 
//	  }
//	  else if( fabs(mcproxylepEta) > 1.56){
//  		if( fakeLepSigma > 0.0298 || fabs(fakeLepdEta) > 0.00609 || fabs(fakeLepdPhi) > 0.045)passFakeProxy=false;
//	  }
//		if(!passFakeProxy)continue;

		if(mcproxysigMET < 20)continue;	

		bool isB(false);
		for(unsigned iMC(0); iMC<mcPID->size(); iMC++){
			double dR = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], mcproxylepEta, mcproxylepPhi);
			if( dR < 0.3 && (fabs( (*mcPID)[iMC]) == 5 ||  (fabs( (*mcMomPID)[iMC]) > 500 && fabs( (*mcMomPID)[iMC]) < 600)))isB = true;
		}

		if(fabs(mcproxydPhiLepMET) < 0.3 || fabs(mcproxydPhiLepMET) > 2.9){
		std::cout << "event " << ievt << " dphi " << fabs(mcproxydPhiLepMET) <<  std::endl;
		std::cout << "lep " << mcproxylepPt << " " << mcproxylepEta << " " << mcproxylepPhi << std::endl;
	//	std::cout << "pho " << mcproxyphoEt << " " << mcproxyphoEta << " " << mcproxyphoPhi << std::endl;
	//	for(unsigned ij(0); ij < jetPt->size(); ij++)std::cout << " jet " << (*jetPt)[ij] << " eta " << (*jetEta)[ij] << " phi " << (*jetPhi)[ij] << std::endl;
		 std::cout << "  MET " << mcproxysigMET << " " << mcproxysigMETPhi <<  std::endl;
		for(unsigned iMC(0); iMC<mcPID->size(); iMC++){
			double dR = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], mcproxylepEta, mcproxylepPhi);
			if( dR < 0.3 || fabs(DeltaPhi(mcproxylepPt,(*mcPhi)[iMC]))>3.0  ) std::cout << (*mcPID)[iMC] << " " << (*mcMomPID)[iMC] << " " << (*mcPt)[iMC] << " " << (*mcEta)[iMC] << " " << (*mcPhi)[iMC] <<  std::endl;
		}
		if(isB)std::cout << "isB " << fabs(mcproxydPhiLepMET) << std::endl;
		 std::cout << std::endl; 
		}

		if(!isB)continue;

		mcpred_PhoEt->Fill(mcproxyphoEt,weight);
		mcpred_PhoEta->Fill(mcproxyphoEta, weight);
		mcpred_MET->Fill(mcproxysigMET, weight);
		mcpred_Mt->Fill(mcproxysigMT, weight);
		mcpred_HT->Fill(mcproxyHT, weight);
		mcpred_LepPt->Fill(mcproxylepPt, weight);
		mcpred_LepEta->Fill(mcproxylepEta, weight);
		mcpred_dPhiEleMET->Fill(fabs(mcproxydPhiLepMET), weight);
		mcpred_nJet->Fill(mcproxynJet, weight);
	}

	mcpred_dPhiEleMET->Draw();

}


