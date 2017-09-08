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

#include "../../include/analysis_rawData.h"
#include "../../include/analysis_photon.h"
#include "../../include/analysis_muon.h"
#include "../../include/analysis_ele.h"
#include "../../include/analysis_mcData.h"
#include "../../include/analysis_tools.h"
#include "../../include/analysis_fakes.h"
#include "../../include/analysis_scalefactor.h"


void analysis_DYFake(){
	
  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");
	
	TH1D *p_predEt = new TH1D("p_predEt","p_predEt",40,0,400);
	TH1D *p_predPt = new TH1D("p_predPt","p_predPt",40,0,400);
	TH1D *p_predMET= new TH1D("p_predMET","predMET",40,0,800);
	TH1D *p_predMT = new TH1D("p_predMT","p_predMT",40,0,800);
	TH1D *p_preddPhi = new TH1D("p_preddPhi"  ,"p_preddPhi",  16, 0,3.2);

	TH1D *p_scaledEt = new TH1D("p_scaledEt","p_scaledEt",40,0,400);
	TH1D *p_scaledPt = new TH1D("p_scaledPt","p_scaledPt",40,0,400);
	TH1D *p_scaledMET= new TH1D("p_scaledMET","scaledMET",40,0,800);
	TH1D *p_scaledMT = new TH1D("p_scaledMT","p_scaledMT",40,0,800);
	TH1D *p_scaleddPhi = new TH1D("p_scaleddPhi"  ,"p_scaleddPhi",  16, 0,3.2);

	TF1 *fitfunc_num = new TF1("fitfunc_num",jetfake_func,35,1000,4);
	TF1 *fitfunc_den = new TF1("fitfunc_den",jetfake_func,35,1000,4);
	std::ifstream jetfakefile("JetFakeRate-transferfactor-DY-Sep1.txt");
	std::string paratype;
	float paravalue;	
	for(int i(0); i < 4; i++){
		jetfakefile >> paratype >> paravalue;
		fitfunc_den->SetParameter(i, paravalue);
	}
	for(int i(0); i < 4; i++){
		jetfakefile >> paratype >> paravalue;
		fitfunc_num->SetParameter(i, paravalue);
	}

// ********  MC *************************//
  TChain *sig;
  sig = new TChain("signalTree");
	sig->Add("/uscms_data/d3/mengleis/test/resTree_mgsignal_DY.root");
  float phoEt(0);
	float phoEta(0);
	float phoPhi(0);
  float lepPt(0);
  float sigMT(0);
  float sigMET(0);
  float dPhiLepMET(0);
  std::vector<int> *mcPID=0;
  std::vector<float> *mcEta=0;
  std::vector<float> *mcPhi=0;
  std::vector<float> *mcPt=0;
  std::vector<int> *mcMomPID=0;

  sig->SetBranchAddress("phoEt",     &phoEt);
  sig->SetBranchAddress("phoEta",    &phoEta);
  sig->SetBranchAddress("phoPhi",    &phoPhi);
  sig->SetBranchAddress("lepPt",     &lepPt);
  sig->SetBranchAddress("sigMT",     &sigMT);
  sig->SetBranchAddress("sigMET",    &sigMET);
  sig->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
  sig->SetBranchAddress("mcPID",     &mcPID);
  sig->SetBranchAddress("mcEta",     &mcEta);
  sig->SetBranchAddress("mcPhi",     &mcPhi);
  sig->SetBranchAddress("mcPt",      &mcPt);
  sig->SetBranchAddress("mcMomPID",  &mcMomPID);

	for(unsigned ievt(0); ievt < sig->GetEntries(); ievt++){
		sig->GetEntry(ievt);

		bool istruepho(false);
		double  mindRpho(0.3);
		unsigned phoIndex(0);
		for(unsigned iMC(0); iMC<mcPID->size(); iMC++){
			double dR = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], phoEta,phoPhi);
			double dE = fabs((*mcPt)[iMC] - phoEt)/phoEt;
			if(dR < mindRpho && dE < 0.5){mindRpho=dR; phoIndex=iMC;}
		}
		if(mindRpho < 0.2){
			if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 23)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 24)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 1)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 2)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 3)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 4)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 5)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 6)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 11)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 13)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 15)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 21)istruepho=true;
			else if((*mcPID)[phoIndex] == 22 && fabs((*mcMomPID)[phoIndex])== 999)istruepho=true;
			else std::cout << (*mcPID)[phoIndex] << " " << fabs((*mcMomPID)[phoIndex]) << std::endl;
		}
		else std::cout << "no match" << std::endl;
		if(!istruepho){
			p_scaledEt->Fill(phoEt);
			p_scaledPt->Fill(lepPt);
			p_scaledMET->Fill(sigMET);
			p_scaledMT->Fill(sigMT);
			p_scaleddPhi->Fill( dPhiLepMET);
		}
	}

  TChain *jet;
  jet = new TChain("jetTree");
	jet->Add("/uscms_data/d3/mengleis/test/resTree_mgsignal_DY.root");
  float jet_phoEt(0);
  float jet_lepPt(0);
  float jet_sigMT(0);
  float jet_sigMET(0);
  float jet_dPhiLepMET(0);
  jet->SetBranchAddress("phoEt",     &jet_phoEt);
  jet->SetBranchAddress("lepPt",     &jet_lepPt);
  jet->SetBranchAddress("sigMT",     &jet_sigMT);
  jet->SetBranchAddress("sigMET",    &jet_sigMET);
  jet->SetBranchAddress("dPhiLepMET",&jet_dPhiLepMET);
	for(unsigned ievt(0); ievt < jet->GetEntries(); ievt++){
		jet->GetEntry(ievt);

		double weight(0);
		weight = fitfunc_num->Eval(phoEt)/fitfunc_den->Eval(phoEt);

		p_predEt->Fill(jet_phoEt, weight); 
		p_predPt->Fill(jet_lepPt, weight);
		p_predMET->Fill(jet_sigMET, weight);
		p_predMT->Fill(jet_sigMT, weight); 
		p_preddPhi->Fill( jet_dPhiLepMET, weight);
	}

	p_predEt->SetLineColor(kRed);
	p_predPt->SetLineColor(kRed);
	p_predMET->SetLineColor(kRed);
	p_predMT->SetLineColor(kRed);
	p_preddPhi->SetLineColor(kRed);

	p_scaledEt->SetLineColor(kBlue);
	p_scaledPt->SetLineColor(kBlue);
	p_scaledMET->SetLineColor(kBlue);
	p_scaledMT->SetLineColor(kBlue);
	p_scaleddPhi->SetLineColor(kBlue);


	TCanvas *can = new TCanvas("can","",1200,800);
	can->Divide(2,3);	
	can->cd(1);
	p_predEt->Draw();
	p_scaledEt->Draw("same");
	can->cd(2);
	p_predPt->Draw();
	p_scaledPt->Draw("same");
	can->cd(3);
	p_predMET->Draw();
	p_scaledMET->Draw("same");
	can->cd(4);
	p_predMT->Draw();
	p_scaledMT->Draw("same");
	can->cd(5);
	p_preddPhi->Draw();
	p_scaleddPhi->Draw("same");
	
	
}

