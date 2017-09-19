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

#include "../../include/analysis_rawData.h"
#include "../../include/analysis_photon.h"
#include "../../include/analysis_muon.h"
#include "../../include/analysis_ele.h"
#include "../../include/analysis_mcData.h"
#include "../../include/analysis_tools.h"
#include "../../include/analysis_fakes.h"
#include "../../include/analysis_scalefactor.h"

void test_WG(){//main 

	int channelType = 2;
  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

	esfScaleFactor  objectESF;

	Double_t plotEtBins[]={35,50,100,150,200,250,300,400,600,800};
	Double_t plotPtBins[]={25,30,35,40,45,50,55,60,65,70,75,80, 85,90,95,100,105,110,115,120,125,130, 135,140,146,152,158,164,170,177,184,192, 200,208,216,224,232,240,250,260,275,290, 305,325,345,370,400};
	TH1D *p_phoEt_data     = new TH1D("p_phoEt_data","",70,100,800);
	TH1D *p_mcEt = new TH1D("p_mcEt","",100,0,400);
	TH1D *p_phoEta_data    = new TH1D("p_phoEta_data","",60,-3,3);
	TH1D *p_phoPhi_data    = new TH1D("p_phoPhi_data","",64,-3.2, 3.2);
	TH1D *p_lepPt_data     = new TH1D("p_lepPt_data","",9,plotEtBins);
	TH1D *p_lepEta_data    = new TH1D("p_lepEta_data","",60,-3,3);
	TH1D *p_lepPhi_data    = new TH1D("p_lepPhi_data","",64,-3.2, 3.2);
	TH1D *p_sigMET_data    = new TH1D("p_sigMET_data","",100,0,400);
	TH1D *p_sigMT_data     = new TH1D("p_sigMT_data","",100,0,400);
	TH1D *p_sigMETPhi_data = new TH1D("p_sigMETPhi_data","",64,-3.2, 3.2);
	TH1D *p_dPhiLepMET_data= new TH1D("p_dPhiLepMET_data","",64,-3.2, 3.2);
	TH1D *p_nVertex_data   = new TH1D("p_nVertex_data","",100,0,100);
	TH1D *p_HT_data        = new TH1D("p_HT_data","",100,0,1000);
	TH1D *p_nJet_data      = new TH1D("p_nJet_data","",10,-0.5,9.5);
	TH1D *p_JetPt         = new TH1D("p_JetPt","p_JetPt",9,plotEtBins);

	TH1D *p_phoEt_ZG     = new TH1D("p_phoEt_ZG","",70,100,800);
	TH1D *p_phoEta_ZG    = new TH1D("p_phoEta_ZG","",60,-3,3);
	TH1D *p_phoPhi_ZG    = new TH1D("p_phoPhi_ZG","",64,-3.2, 3.2);
	TH1D *p_lepPt_ZG     = new TH1D("p_lepPt_ZG","",9,plotEtBins);
	TH1D *p_lepEta_ZG    = new TH1D("p_lepEta_ZG","",60,-3,3);
	TH1D *p_lepPhi_ZG    = new TH1D("p_lepPhi_ZG","",64,-3.2, 3.2);
	TH1D *p_sigMET_ZG    = new TH1D("p_sigMET_ZG","",100,0,400);
	TH1D *p_sigMT_ZG     = new TH1D("p_sigMT_ZG","",100,0,400);
	TH1D *p_sigMETPhi_ZG = new TH1D("p_sigMETPhi_ZG","",64,-3.2, 3.2);
	TH1D *p_dPhiLepMET_ZG= new TH1D("p_dPhiLepMET_ZG","",64,-3.2, 3.2);
	TH1D *p_nVertex_ZG   = new TH1D("p_nVertex_ZG","",100,0,100);
	TH1D *p_HT_ZG        = new TH1D("p_HT_ZG","",100,0,1000);
	TH1D *p_nJet_ZG      = new TH1D("p_nJet_ZG","",10,-0.5,9.5);
	TH1D *p_JetPt_ZG     = new TH1D("p_JetPt_ZG","p_JetPt_ZG",9,plotEtBins);

//************ Signal Tree **********************//
  TChain *tree = new TChain("mgTree");
  tree->Add("/uscms_data/d3/mengleis/Sep13/resTree_VGamma_WG_Pt130.root");
	float MCweight(0);
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
	float JetPt(0);

	
  tree->SetBranchAddress("MCweight",  &MCweight);
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
  tree->SetBranchAddress("HT",        &HT);
  tree->SetBranchAddress("nJet",      &nJet);
	tree->SetBranchAddress("ISRJetPt",     &JetPt);
	
  for(unsigned ievt(0); ievt<tree->GetEntries(); ++ievt){//loop on entries
		tree->GetEntry(ievt);

		double weight = 1; 

		if(fabs(phoEta) > 1.4442)continue;
		if( phoEt < 140)continue;

		p_phoEt_data->Fill(phoEt, weight);	
		p_lepPt_data->Fill(lepPt, weight);
		p_lepEta_data->Fill(lepEta, weight);
		p_lepPhi_data->Fill(lepPhi, weight);
		p_sigMET_data->Fill(sigMET, weight);
		p_sigMT_data->Fill(sigMT, weight);
		p_sigMETPhi_data->Fill(sigMETPhi, weight);
		p_dPhiLepMET_data->Fill(dPhiLepMET, weight);
		p_nVertex_data->Fill(nVertex, weight);
		p_HT_data->Fill(HT, weight);
		p_nJet_data->Fill(nJet, weight);
		p_JetPt->Fill(JetPt, weight);
	}//loop on  events




//************ Signal Tree **********************//
  TChain *ZGtree = new TChain("mgTree");
  ZGtree->Add("/uscms_data/d3/mengleis/Sep13/resTree_VGamma_ZG.root");
	float ZG_MCweight(0);
  float ZG_phoEt(0);
  float ZG_phoEta(0);
  float ZG_phoPhi(0);
  float ZG_lepPt(0);
  float ZG_lepEta(0);
  float ZG_lepPhi(0);
  float ZG_sigMT(0);
  float ZG_sigMET(0);
  float ZG_sigMETPhi(0);
  float ZG_dPhiLepMET(0);
  int   ZG_nVertex(0);
  float ZG_HT(0);
  float ZG_nJet(0);
	float ZG_JetPt(0);

  ZGtree->SetBranchAddress("MCweight",  &ZG_MCweight);
  ZGtree->SetBranchAddress("phoEt",     &ZG_phoEt);
  ZGtree->SetBranchAddress("phoEta",    &ZG_phoEta);
  ZGtree->SetBranchAddress("phoPhi",    &ZG_phoPhi);
  ZGtree->SetBranchAddress("lepPt",     &ZG_lepPt);
  ZGtree->SetBranchAddress("lepEta",    &ZG_lepEta);
  ZGtree->SetBranchAddress("lepPhi",    &ZG_lepPhi);
  ZGtree->SetBranchAddress("sigMT",     &ZG_sigMT);
  ZGtree->SetBranchAddress("sigMET",    &ZG_sigMET);
  ZGtree->SetBranchAddress("sigMETPhi", &ZG_sigMETPhi);
  ZGtree->SetBranchAddress("dPhiLepMET",&ZG_dPhiLepMET);
  ZGtree->SetBranchAddress("nVertex",   &ZG_nVertex);
  ZGtree->SetBranchAddress("HT",        &ZG_HT);
  ZGtree->SetBranchAddress("nJet",      &ZG_nJet);
	ZGtree->SetBranchAddress("ISRJetPt",     &ZG_JetPt);

  for(unsigned ievt(0); ievt<ZGtree->GetEntries(); ++ievt){//loop on entries
		ZGtree->GetEntry(ievt);

		double weight = 1; 
		if(fabs(ZG_phoEta) > 1.4442)continue;
		if( ZG_phoEt < 140)continue;

		p_phoEt_ZG->Fill(ZG_phoEt, weight);
		p_lepPt_ZG->Fill(ZG_lepPt, weight);
		p_lepEta_ZG->Fill(ZG_lepEta, weight);
		p_lepPhi_ZG->Fill(ZG_lepPhi, weight);
		p_sigMET_ZG->Fill(ZG_sigMET, weight);
		p_sigMT_ZG->Fill(ZG_sigMT, weight);
		p_sigMETPhi_ZG->Fill(ZG_sigMETPhi, weight);
		p_dPhiLepMET_ZG->Fill(ZG_dPhiLepMET, weight);
		p_nVertex_ZG->Fill(ZG_nVertex, weight);
		p_HT_ZG->Fill(ZG_HT, weight);
		p_nJet_ZG->Fill(ZG_nJet, weight);
		p_JetPt_ZG->Fill(ZG_JetPt, weight);
	}//loop on  events

	TCanvas *can_phoEt     = new TCanvas("can_phoEt",       "can_phoEt",600,600); 
	TCanvas *can_lepPt     = new TCanvas("can_lepPt",       "can_lepPt",600,600); 
	TCanvas *can_lepEta    = new TCanvas("can_lepEta",      "can_lepEta",600,600); 
	TCanvas *can_trail    = new TCanvas("can_trail",      "can_trail",600,600); 
	TCanvas *can_sigMET    = new TCanvas("can_sigMET",      "can_sigMET",600,600); 
	TCanvas *can_sigMT     = new TCanvas("can_sigMT",       "can_sigMT",600,600); 
	TCanvas *can_sigMETPhi = new TCanvas("can_sigMETPhi",   "can_sigMETPhi",600,600); 
	TCanvas *can_dPhiLepMET= new TCanvas("can_dPhiLepMET",  "can_dPhiLepMET",600,600); 
	TCanvas *can_nVertex   = new TCanvas("can_nVertex",     "can_nVertex",600,600); 
	TCanvas *can_dRPhoLep  = new TCanvas("can_dRPhoLep",    "can_dRPhoLep",600,600); 
	TCanvas *can_HT        = new TCanvas("can_HT",          "can_HT",600,600); 
	TCanvas *can_nJet      = new TCanvas("can_nJet",        "can_nJet",600,600); 
	TCanvas *can_invmass   = new TCanvas("can_invmass",     "can_invmass", 600,600);    
	TCanvas *can_llmass   = new TCanvas("can_llmass",     "can_llmass", 600,600);   

	float scalefactor = p_lepPt_data->Integral(1,10)/p_lepPt_ZG->Integral(1,10);
	can_phoEt->cd();
	gPad->SetLogy();
	p_phoEt_data->Draw();
	p_phoEt_ZG->SetLineColor(kRed);
	p_phoEt_ZG->Scale(scalefactor);
	p_phoEt_ZG->Draw("same");

 
	can_lepPt->cd();
	gPad->SetLogy();
	p_lepPt_data->Draw();
	p_lepPt_ZG->SetLineColor(kRed);
	p_lepPt_ZG->Scale(scalefactor);
	p_lepPt_ZG->Draw("same");

	can_lepEta->cd();
	p_lepEta_data->Draw();
	p_lepEta_ZG->SetLineColor(kRed);
	p_lepEta_ZG->Scale(scalefactor);
	p_lepEta_ZG->Draw("same");
	for(unsigned ibin(1); ibin < 60; ibin++){
		if(p_lepEta_data->GetBinContent(ibin) > 0)std::cout << p_lepEta_ZG->GetBinCenter(ibin) << " " << p_lepEta_ZG->GetBinContent(ibin)/p_lepEta_data->GetBinContent(ibin) << std::endl;
	} 

	can_sigMET->cd();
	p_sigMET_data->Draw();
	p_sigMET_ZG->SetLineColor(kRed);
	p_sigMET_ZG->Scale(scalefactor);
	p_sigMET_ZG->Draw("same");

	can_sigMT->cd();
	p_sigMT_data->Draw();
	p_sigMT_ZG->SetLineColor(kRed);
	p_sigMT_ZG->Scale(scalefactor);
	p_sigMT_ZG->Draw("same");

	can_sigMETPhi->cd();
	p_sigMETPhi_data->Draw();
	p_sigMETPhi_ZG->SetLineColor(kRed);
	p_sigMETPhi_ZG->Scale(scalefactor);
	p_sigMETPhi_ZG->Draw("same");

	can_dPhiLepMET->cd();
	p_dPhiLepMET_data->Draw();
	p_dPhiLepMET_ZG->SetLineColor(kRed);
	p_dPhiLepMET_ZG->Scale(scalefactor);
	p_dPhiLepMET_ZG->Draw("same");

	can_nVertex->cd();
	p_nVertex_data->Draw();
	p_nVertex_ZG->SetLineColor(kRed);
	p_nVertex_ZG->Scale(scalefactor);
	p_nVertex_ZG->Draw("same");

	can_HT->cd();
	gPad->SetLogy();
	p_HT_data->Draw();
	p_HT_ZG->SetLineColor(kRed);
	p_HT_ZG->Scale(scalefactor);
	p_HT_ZG->Draw("same");

	can_nJet->cd();
	p_JetPt->Draw();
	p_JetPt_ZG->Scale(scalefactor);
	p_JetPt_ZG->SetLineColor(kRed);
	p_JetPt_ZG->Draw("same");


}


