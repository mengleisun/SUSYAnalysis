#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<ctime>

#include "TFile.h"
#include "TTree.h"
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
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TProfile2D.h"

#include "../include/analysis_rawData.h"
#include "../include/analysis_photon.h"
#include "../include/analysis_muon.h"
#include "../include/analysis_ele.h"
#include "../include/analysis_mcData.h"
#include "../include/analysis_tools.h"
#include "../include/analysis_fakes.h"
#include "../include/analysis_scalefactor.h"

void compareTree(){//main 

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

	std::ostringstream treeName;
	treeName.str("");
	treeName << "jetTree";

	bool   EBOnly = true;
	double scalefactor1 = 1; 
	double scalefactor2 = 1; 

	Double_t plotEtBins[]={35,50,100,150,200,250,300,400,600,800};
	Double_t plotPtBins[]={25,30,35,40,45,50,55,60,65,70,75,80, 85,90,95,100,105,110,115,120,125,130, 135,140,146,152,158,164,170,177,184,192, 200,208,216,224,232,240,250,260,275,290, 305,325,345,370,400};
	TH1D *p_phoEt_data     = new TH1D("p_phoEt_data","",9,plotEtBins);
	TH1D *p_mcEt = new TH1D("p_mcEt","",100,0,400);
	TH1D *p_phoEta_data    = new TH1D("p_phoEta_data","",60,-3,3);
	TH1D *p_phoPhi_data    = new TH1D("p_phoPhi_data","",64,-3.2, 3.2);
	TH1D *p_lepPt_data     = new TH1D("p_lepPt_data","",9,plotEtBins);
	TH1D *p_lepEta_data    = new TH1D("p_lepEta_data","",60,-3,3);
	TH1D *p_lepPhi_data    = new TH1D("p_lepPhi_data","",64,-3.2, 3.2);
	TH1D *p_sigMET_data    = new TH1D("p_sigMET_data","",100,0,1000);
	TH1D *p_sigMT_data     = new TH1D("p_sigMT_data","",100,0,1000);
	TH1D *p_sigMETPhi_data = new TH1D("p_sigMETPhi_data","",64,-3.2, 3.2);
	TH1D *p_dPhiLepMET_data= new TH1D("p_dPhiLepMET_data","",64,-3.2, 3.2);
	TH1D *p_nVertex_data   = new TH1D("p_nVertex_data","",100,0,100);
	TH1D *p_dRPhoLep_data  = new TH1D("p_dRPhoLep_data","",60,0.3,3.3);
	TH1D *p_HT_data        = new TH1D("p_HT_data","",100,0,1000);
	TH1D *p_nJet_data      = new TH1D("p_nJet_data","",10,0,10);
	TH1D *p_llmass_data   = new TH1D("p_llmass_data","",130,0,130);

	TH1D *p_phoEt_2     = new TH1D("p_phoEt_2","",9,plotEtBins);
	TH1D *p_phoEta_2    = new TH1D("p_phoEta_2","",60,-3,3);
	TH1D *p_phoPhi_2    = new TH1D("p_phoPhi_2","",64,-3.2, 3.2);
	TH1D *p_lepPt_2     = new TH1D("p_lepPt_2","",9,plotEtBins);
	TH1D *p_lepEta_2    = new TH1D("p_lepEta_2","",60,-3,3);
	TH1D *p_lepPhi_2    = new TH1D("p_lepPhi_2","",64,-3.2, 3.2);
	TH1D *p_sigMET_2    = new TH1D("p_sigMET_2","",100,0,1000);
	TH1D *p_sigMT_2     = new TH1D("p_sigMT_2","",100,0,1000);
	TH1D *p_sigMETPhi_2 = new TH1D("p_sigMETPhi_2","",64,-3.2, 3.2);
	TH1D *p_dPhiLepMET_2= new TH1D("p_dPhiLepMET_2","",64,-3.2, 3.2);
	TH1D *p_nVertex_2   = new TH1D("p_nVertex_2","",100,0,100);
	TH1D *p_dRPhoLep_2  = new TH1D("p_dRPhoLep_2","",60,0.3,3.3);
	TH1D *p_HT_2        = new TH1D("p_HT_2","",100,0,1000);
	TH1D *p_nJet_2      = new TH1D("p_nJet_2","",10,0,10);
	TH1D *p_llmass_2   = new TH1D("p_llmass_2","",130,0,130);

	int  Nsig_1(0);
	int  Nsig_2(0);
//************ Signal Tree **********************//
  TChain *tree = new TChain(treeName.str().c_str());
  tree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_egsignal_DoubleEG_ReMiniAOD_FullEcal_newEta.root");
  float phoEt(0);
  float phoEta(0);
  float phoPhi(0);
  float lepPt(0);
  float lepEta(0);
  float lepPhi(0);
  float sigMT(0);
  float sigMET(0);
  float sigMETPhi(0);
  float dPhiLepMET(0);
  int   nVertex(0);
  float dRPhoLep(0);
  float HT(0);
  float nJet(0);
	float llmass(0);
	float threeMass(0);

  tree->SetBranchAddress("phoEt",     &phoEt);
  tree->SetBranchAddress("phoEta",    &phoEta);
  tree->SetBranchAddress("phoPhi",    &phoPhi);
  tree->SetBranchAddress("lepPt",     &lepPt);
  tree->SetBranchAddress("lepEta",    &lepEta);
  tree->SetBranchAddress("lepPhi",    &lepPhi);
  tree->SetBranchAddress("sigMT",     &sigMT);
  tree->SetBranchAddress("sigMET",    &sigMET);
  tree->SetBranchAddress("sigMETPhi", &sigMETPhi);
  tree->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
  tree->SetBranchAddress("nVertex",   &nVertex);
  tree->SetBranchAddress("dRPhoLep",  &dRPhoLep);
  tree->SetBranchAddress("HT",        &HT);
  tree->SetBranchAddress("nJet",      &nJet);
	tree->SetBranchAddress("llmass",    &llmass);
	tree->SetBranchAddress("threeMass", &threeMass);

  for(unsigned ievt(0); ievt<tree->GetEntries(); ++ievt){//loop on entries
		tree->GetEntry(ievt);
 
		double weight = scalefactor1; 
		if(EBOnly && fabs(phoEta) > 1.4442)continue;

		if(sigMT > 100 && sigMET > 120)Nsig_1 += 1;
		p_phoEt_data->Fill(phoEt, weight); 
		p_phoEta_data->Fill(phoEta, weight);
		p_phoPhi_data->Fill(phoPhi, weight);
		p_lepPt_data->Fill(lepPt, weight);
		p_lepEta_data->Fill(lepEta, weight);
		p_lepPhi_data->Fill(lepPhi, weight);
		p_sigMET_data->Fill(sigMET, weight);
		p_sigMT_data->Fill(sigMT, weight);
		p_sigMETPhi_data->Fill(sigMETPhi, weight);
		p_dPhiLepMET_data->Fill(dPhiLepMET, weight);
		p_nVertex_data->Fill(nVertex, weight);
		p_dRPhoLep_data->Fill(dRPhoLep, weight);
		p_HT_data->Fill(HT, weight);
		p_nJet_data->Fill(nJet, weight);
		p_llmass_data->Fill(llmass, weight);
	}//loop on  events


//************ Signal Tree **********************//
  TChain *tree2 = new TChain(treeName.str().c_str());
  tree2->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_egsignal_SingleEvent.root");
  float second_phoEt(0);
  float second_phoEta(0);
  float second_phoPhi(0);
  float second_lepPt(0);
  float second_lepEta(0);
  float second_lepPhi(0);
  float second_sigMT(0);
  float second_sigMET(0);
  float second_sigMETPhi(0);
  float second_dPhiLepMET(0);
  int   second_nVertex(0);
  float second_dRPhoLep(0);
  float second_dRPhoTrail(0);
  float second_HT(0);
  float second_nJet(0);
	float second_llmass(0);
	float second_threeMass(0);

  tree2->SetBranchAddress("phoEt",     &second_phoEt);
  tree2->SetBranchAddress("phoEta",    &second_phoEta);
  tree2->SetBranchAddress("phoPhi",    &second_phoPhi);
  tree2->SetBranchAddress("lepPt",     &second_lepPt);
  tree2->SetBranchAddress("lepEta",    &second_lepEta);
  tree2->SetBranchAddress("lepPhi",    &second_lepPhi);
  tree2->SetBranchAddress("sigMT",     &second_sigMT);
  tree2->SetBranchAddress("sigMET",    &second_sigMET);
  tree2->SetBranchAddress("sigMETPhi", &second_sigMETPhi);
  tree2->SetBranchAddress("dPhiLepMET",&second_dPhiLepMET);
  tree2->SetBranchAddress("nVertex",   &second_nVertex);
  tree2->SetBranchAddress("dRPhoLep",  &second_dRPhoLep);
  tree2->SetBranchAddress("dRPhoTrail",&second_dRPhoTrail);
  tree2->SetBranchAddress("HT",        &second_HT);
  tree2->SetBranchAddress("nJet",      &second_nJet);
	tree2->SetBranchAddress("llmass",    &second_llmass);
	tree2->SetBranchAddress("threeMass", &second_threeMass);

  for(unsigned ievt(0); ievt<tree2->GetEntries(); ++ievt){//loop on entries
		tree2->GetEntry(ievt);

		double weight = scalefactor2; 
		if(EBOnly && fabs(second_phoEta) > 1.4442)continue;

		if(second_sigMT > 100 && second_sigMET > 120)Nsig_2 += 1;
		p_phoEt_2->Fill(second_phoEt, weight); 
		p_phoEta_2->Fill(second_phoEta, weight);
		p_phoPhi_2->Fill(second_phoPhi, weight);
		p_lepPt_2->Fill(second_lepPt, weight);
		p_lepEta_2->Fill(second_lepEta, weight);
		p_lepPhi_2->Fill(second_lepPhi, weight);
		p_sigMET_2->Fill(second_sigMET, weight);
		p_sigMT_2->Fill(second_sigMT, weight);
		p_sigMETPhi_2->Fill(second_sigMETPhi, weight);
		p_dPhiLepMET_2->Fill(second_dPhiLepMET, weight);
		p_nVertex_2->Fill(second_nVertex, weight);
		p_dRPhoLep_2->Fill(second_dRPhoLep, weight);
		p_HT_2->Fill(second_HT, weight);
		p_nJet_2->Fill(second_nJet, weight);
		p_llmass_2->Fill(second_llmass, weight);
	}//loop on  events
	
	std::cout << "sig1 = " << Nsig_1 << "  sig2 = " << Nsig_2 << std::endl;

	TCanvas *can_phoEt     = new TCanvas("can_phoEt",       "can_phoEt", 600,600); 
	TCanvas *can_phoEta    = new TCanvas("can_phoEta",      "can_phoEta",600,600); 
	TCanvas *can_phoPhi    = new TCanvas("can_phoPhi",      "can_phoPhi",600,600); 
	TCanvas *can_lepPt     = new TCanvas("can_lepPt",       "can_lepPt",600,600); 
	TCanvas *can_lepEta    = new TCanvas("can_lepEta",      "can_lepEta",600,600); 
	TCanvas *can_lepPhi    = new TCanvas("can_lepPhi",      "can_lepPhi",600,600); 
	TCanvas *can_sigMET    = new TCanvas("can_sigMET",      "can_sigMET",600,600); 
	TCanvas *can_sigMT     = new TCanvas("can_sigMT",       "can_sigMT",600,600); 
	TCanvas *can_sigMETPhi = new TCanvas("can_sigMETPhi",   "can_sigMETPhi",600,600); 
	TCanvas *can_dPhiLepMET= new TCanvas("can_dPhiLepMET",  "can_dPhiLepMET",600,600); 
	TCanvas *can_nVertex   = new TCanvas("can_nVertex",     "can_nVertex",600,600); 
	TCanvas *can_dRPhoLep  = new TCanvas("can_dRPhoLep",    "can_dRPhoLep",600,600); 
	TCanvas *can_HT        = new TCanvas("can_HT",          "can_HT",600,600); 
	TCanvas *can_nJet      = new TCanvas("can_nJet",        "can_nJet",600,600); 
	TCanvas *can_llmass   = new TCanvas("can_llmass",     "can_llmass", 600,600);   

	float scalefactor = p_phoEt_data->Integral(4,5)/p_phoEt_2->Integral(4,5); 
	std::cout << "scale=" << scalefactor << std::endl;
	//float scalefactor = 1;	

	can_phoEt->cd();
	TPad *phoEt_pad1 = new TPad("phoEt_pad1", "phoEt_pad1", 0, 0.3, 1, 1.0);
	phoEt_pad1->SetBottomMargin(0.1);
	phoEt_pad1->Draw();  
	phoEt_pad1->cd();  
	phoEt_pad1->SetLogy();
	p_phoEt_data->Draw();
	p_phoEt_2->SetLineColor(kRed);
	p_phoEt_2->Scale(scalefactor);
	p_phoEt_2->Draw("same");
	TLegend *leg=new TLegend(0.6,0.9,0.75,0.9);
	leg->AddEntry(p_phoEt_data,"observed");
	leg->AddEntry(p_phoEt_2,"ZG");
	leg->Draw("same");

	can_phoEt->cd();
	TPad *phoEt_pad2 = new TPad("phoEt_pad2", "phoEt_pad2", 0, 0.05, 1, 0.25);
	phoEt_pad2->Draw();
	phoEt_pad2->cd();
	TLine *flatratio_eventcount = new TLine(0,1,600,1);
	TH1D *ratio_phoEt = (TH1D*)p_phoEt_2->Clone("ratio_phoEt");
	ratio_phoEt->Divide(p_phoEt_data);
	ratio_phoEt->Draw();
	flatratio_eventcount->Draw("same");
	for(int ibin(1); ibin < ratio_phoEt->GetSize(); ibin++){
		std::cout << ibin << " " << ratio_phoEt->GetBinContent(ibin) << std::endl;
	}
	
	std::cout << " tree1-tree2 = " << p_phoEta_data->GetEntries() - p_phoEta_2->GetEntries() << std::endl;
	std::cout << " reduce = " << (p_phoEta_data->GetEntries() - p_phoEta_2->GetEntries())/p_phoEta_data->GetEntries() << std::endl;

	can_phoEta->cd();
	p_phoEta_data->Draw();
	p_phoEta_2->SetLineColor(kRed);
	p_phoEta_2->Scale(scalefactor);
	p_phoEta_2->Draw("same");

	can_phoPhi->cd();
	p_phoPhi_data->Draw();
	p_phoPhi_2->SetLineColor(kRed);
	p_phoPhi_2->Scale(scalefactor);
	p_phoPhi_2->Draw("same");

	can_lepPt->cd();
	gPad->SetLogy();
	p_lepPt_data->Draw();
	p_lepPt_2->SetLineColor(kRed);
	p_lepPt_2->Scale(scalefactor);
	p_lepPt_2->Draw("same");

	can_lepEta->cd();
	p_lepEta_data->Draw();
	p_lepEta_2->SetLineColor(kRed);
	p_lepEta_2->Scale(scalefactor);
	p_lepEta_2->Draw("same");
	for(unsigned ibin(1); ibin < 60; ibin++){
		if(p_lepEta_data->GetBinContent(ibin) > 0)std::cout << p_lepEta_2->GetBinCenter(ibin) << " " << p_lepEta_2->GetBinContent(ibin)/p_lepEta_data->GetBinContent(ibin) << std::endl;
	} 

	can_lepPhi->cd();
	p_lepPhi_data->Draw();
	p_lepPhi_2->SetLineColor(kRed);
	p_lepPhi_2->Scale(scalefactor);
	p_lepPhi_2->Draw("same");

	can_sigMET->cd();
	p_sigMET_data->Draw();
	p_sigMET_2->SetLineColor(kRed);
	p_sigMET_2->Scale(scalefactor);
	p_sigMET_2->Draw("same");

	can_sigMT->cd();
	p_sigMT_data->Rebin(10);
	p_sigMT_2->Rebin(10);
	p_sigMT_data->Draw();
	p_sigMT_2->SetLineColor(kRed);
	p_sigMT_2->Scale(scalefactor);
	p_sigMT_2->Draw("same");

	can_sigMETPhi->cd();
	p_sigMETPhi_data->Draw();
	p_sigMETPhi_2->SetLineColor(kRed);
	p_sigMETPhi_2->Scale(scalefactor);
	p_sigMETPhi_2->Draw("same");

	can_dPhiLepMET->cd();
	p_dPhiLepMET_data->Draw();
	p_dPhiLepMET_2->SetLineColor(kRed);
	p_dPhiLepMET_2->Scale(scalefactor);
	p_dPhiLepMET_2->Draw("same");

	can_nVertex->cd();
	p_nVertex_data->Draw();
	p_nVertex_2->SetLineColor(kRed);
	p_nVertex_2->Scale(scalefactor);
	p_nVertex_2->Draw("same");

	can_dRPhoLep->cd();
	p_dRPhoLep_data->Draw();
	p_dRPhoLep_2->SetLineColor(kRed);
	p_dRPhoLep_2->Scale(scalefactor);
	p_dRPhoLep_2->Draw("same");

	can_HT->cd();
	gPad->SetLogy();
	p_HT_data->Draw();
	p_HT_2->SetLineColor(kRed);
	p_HT_2->Scale(scalefactor);
	p_HT_2->Draw("same");

	can_nJet->cd();
	p_nJet_data->Draw();
	p_nJet_2->SetLineColor(kRed);
	p_nJet_2->Scale(scalefactor);
	p_nJet_2->Draw("same");

	can_llmass->cd();
	p_llmass_data->Draw();
	p_llmass_2->SetLineColor(kRed);
	p_llmass_2->Scale(scalefactor);
	p_llmass_2->Draw("same");

}


