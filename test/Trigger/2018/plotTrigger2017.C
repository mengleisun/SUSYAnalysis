#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include "TH1F.h"
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TF1.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "../../../include/tdrstyle.C"

void plotTrigger2017(){//main  
	setTDRStyle();  
	gStyle->SetTitleXOffset(1); 
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0.5);
  Double_t xaxis2d[] = {25,40,70,200};
  Double_t yaxis2d[] = {10,20,70,200};
  Double_t plotPtBins[]={10,15,20,25,30,35,40,50,60,70,100,200};
	TEfficiency *p_HLT_pho = new TEfficiency("p_L1_pho","HLT/L1 efficiency;p_{T} (GeV);HLT/L1 efficiency",11, plotPtBins);
	TEfficiency *p_HLT_mu = new TEfficiency("p_L1_mu","HLT/L1 efficiency;p_{T} (GeV);HLT/L1 efficiency",11, plotPtBins);
	TProfile2D *p_L1eff_Z  = new TProfile2D("p_L1eff_Z", "L1 efficiency;#gamma p_{T} (GeV); #mu p_{T} (GeV)",3,xaxis2d,3,yaxis2d);

	p_L1eff_Z->GetXaxis()->SetMoreLogLabels(); 
	p_L1eff_Z->GetYaxis()->SetMoreLogLabels(); 
	

	TChain *Ztree = new TChain("ZTree","ZTree");
	Ztree->Add("../plot_MuonTrigger2018.root");
	float Z_phoEt(0);
	float Z_phoEta(0);
	float Z_muPt(0);
	int   Z_phofireL1;
	int   Z_mufireL1;
	int   Z_phofireHLT;
	int   Z_mufireHLT;

	Ztree->SetBranchAddress("phoEt",     &Z_phoEt);
	Ztree->SetBranchAddress("phoEta",    &Z_phoEta);
	Ztree->SetBranchAddress("muPt",      &Z_muPt);
	Ztree->SetBranchAddress("phofireL1", &Z_phofireL1);
	Ztree->SetBranchAddress("mufireL1",  &Z_mufireL1); 
	Ztree->SetBranchAddress("phofireHLT",&Z_phofireHLT);
	Ztree->SetBranchAddress("mufireHLT", &Z_mufireHLT);

	double pass(0);
	double total(0);
	for(unsigned ievt(0); ievt < Ztree->GetEntries(); ievt++){
		Ztree->GetEntry(ievt);
		if(Z_phoEt > 200)Z_phoEt = 199;
		if(Z_muPt > 200)Z_muPt = 199;

		if(fabs(Z_phoEta) > 2.5)continue;
		if(Z_phofireL1 > 0 && Z_mufireL1 > 0)
			p_L1eff_Z->Fill(Z_phoEt, Z_muPt, 1);
		else p_L1eff_Z->Fill(Z_phoEt, Z_muPt, 0);

		if(Z_phofireL1 > 0 && Z_mufireL1 > 0){
			p_HLT_pho->Fill(Z_phofireHLT > 0, Z_phoEt);
			p_HLT_mu->Fill(Z_mufireHLT > 0, Z_muPt);
		}
		if(Z_phoEt > 40 && Z_muPt > 20){
			total += 1;
			if(Z_phofireL1  > 0 && Z_mufireL1  > 0)pass+=1;
		}
	}
	std::cout << "eff: " << pass/total << std::endl;

	TCanvas *canZ = new TCanvas("canZ","",600,600);
	canZ->SetLeftMargin(0.2);
	canZ->SetRightMargin(0.2);
	canZ->SetTopMargin(0.1);
	canZ->SetBottomMargin(0.15);
	canZ->cd();
	gPad->SetLogx();
	gPad->SetLogy();
	gStyle->SetPaintTextFormat("4.2f");
	Int_t PaletteColors[] = {9, kBlue, kBlue-4,kCyan, kTeal, kGreen,kSpring, 5, 2};
	gStyle->SetPalette(9, PaletteColors);
	p_L1eff_Z->Draw("E colz text");
	canZ->SaveAs("mgTrigger_efficiency_2017.pdf");


	TProfile *p_mutree_phoL1 = new TProfile("p_mutree_phoL1","p_mutree_phoHLT",11, plotPtBins);
	TProfile *p_mutree_phoHLT = new TProfile("p_mutree_phoHLT","p_mutree_phoHLT",11, plotPtBins);
	TChain *mutree = new TChain("mgTree","mgTree");
	mutree->Add("../plot_MuonTrigger2018.root");
	float mu_phoEt(0);
	float mu_phoEta(0);
	float mu_muPt(0);
	int   mu_phofireL1;
	int   mu_mufireL1;
	int   mu_phofireHLT;
	int   mu_mufireHLT;

	mutree->SetBranchAddress("phoEt",     &mu_phoEt);
	mutree->SetBranchAddress("phoEta",    &mu_phoEta);
	mutree->SetBranchAddress("muPt",      &mu_muPt);
	mutree->SetBranchAddress("phofireL1", &mu_phofireL1);
	mutree->SetBranchAddress("mufireL1",  &mu_mufireL1); 
	mutree->SetBranchAddress("phofireHLT",&mu_phofireHLT);
	mutree->SetBranchAddress("mufireHLT", &mu_mufireHLT);

	for(unsigned ievt(0); ievt < mutree->GetEntries(); ievt++){
		mutree->GetEntry(ievt);
		if(mu_phoEt > 200)mu_phoEt = 199;
		if(mu_muPt > 200)mu_muPt = 199;

		if(fabs(mu_phoEta) > 1.4)continue;
		if(mu_muPt < 40)continue;
		if(mu_phofireL1 > 0)p_mutree_phoL1->Fill(mu_phoEt, 1);
		else p_mutree_phoL1->Fill(mu_phoEt, 0);

		//if(mu_phofireL1 > 0 && mu_mufireL1 > 0){
			if(mu_phofireHLT > 0)p_mutree_phoHLT->Fill(mu_phoEt, 1);
			else p_mutree_phoHLT->Fill(mu_phoEt, 0);
		//}
	}
	TCanvas *canmu = new TCanvas("canmu","",600,600);
	canmu->cd();
	TLegend *leg=new TLegend(0.5,0.65,0.87,0.8);
	p_HLT_pho->SetLineColor(kRed);
	p_HLT_pho->SetMarkerColor(kRed);
	leg->AddEntry(p_HLT_pho,"photon leg");
	leg->AddEntry(p_HLT_mu,"muon leg");
	p_HLT_pho->Draw("AP");
	p_HLT_mu->Draw("P same");
	leg->Draw("same");
	canmu->SaveAs("HLTeff.pdf");
}

