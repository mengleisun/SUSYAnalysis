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
#include "../include/tdrstyle.C"

void plotTrigger2017(){//main  
  gStyle->SetOptStat(0);
	setTDRStyle();   
  gStyle->SetErrorX(0.5);
  Double_t xaxis2d[] = {35,40,45,60,200};
  Double_t yaxis2d[] = {25,30,45,60,200};
  Double_t plotPtBins[]={10,15,20,25,30,35,40,50,60,70,100,200};
	TProfile *p_L1_pho = new TProfile("p_L1_pho","p_L1_pho",11, plotPtBins);
	TProfile2D *p_HLTeff_Z  = new TProfile2D("p_HLTeff_Z", "#mu#gamma trigger efficiency;#gamma p_{T} (GeV); #mu p_{T} (GeV)",4,xaxis2d,4,yaxis2d);
	TProfile2D *p_HLTeff_mg = new TProfile2D("p_HLTeff_mg","MET dataset",4,xaxis2d,4,yaxis2d);
	TProfile2D *p_HLTeff_DY = new TProfile2D("p_HLTeff_DY","DY efficiency",4,xaxis2d,4,yaxis2d);
	//TH2F 			 *p_crosseff = new TH2F("p_crosseff","MuonEG trigger efficiency; #gamma p_{T} (GeV); muon p_{T} (GeV)",4,xaxis2d,4,yaxis2d); 
	TH2F 			 *p_crosseff = new TH2F("p_crosseff","MuonEG trigger efficiency; #gamma p_{T} (GeV); muon p_{T} (GeV)",10,0,200,10,0,200);
	TH2F 			 *p_mgESF = new TH2F("p_mgESF","MuonEG trigger ESF; #gamma p_{T} (GeV); muon p_{T} (GeV)",4,xaxis2d,4,yaxis2d); 

	TChain *Ztree = new TChain("ZTree","ZTree");
	Ztree->Add("/uscms_data/d3/mengleis/FullStatusOct/plot_MuonTrigger2017.root");
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

	for(unsigned ievt(0); ievt < Ztree->GetEntries(); ievt++){
		Ztree->GetEntry(ievt);
		if(Z_phoEt > 200)Z_phoEt = 199;
		if(Z_muPt > 200)Z_muPt = 199;

		if( Z_phoEt > 25 ){
			if(Z_mufireHLT > 0)p_L1_pho->Fill(Z_muPt, 1);
			else p_L1_pho->Fill(Z_muPt, 0);
		}
//		if(Z_phofireL1 > 0 && Z_mufireL1 > 0){
//		if(Z_phofireHLT > 0 && Z_mufireHLT > 0)
//			p_HLTeff_Z->Fill(Z_phoEt, Z_muPt, 1);
//		else p_HLTeff_Z->Fill(Z_phoEt, Z_muPt, 0);}
//		if(Z_phofireL1 > 0 && Z_mufireL1 > 0)
//			p_HLTeff_Z->Fill(Z_phoEt, Z_muPt, 1);
//		else p_HLTeff_Z->Fill(Z_phoEt, Z_muPt, 0);
	}


	TCanvas *canmg = new TCanvas("canmg","",600,600);
	canmg->cd();
	p_HLTeff_mg->Draw("colz text");
  for(int binx(1); binx <= 4; binx++){
		for(int biny(1); biny <= 4; biny++){
			p_crosseff->SetBinContent(binx, biny, p_HLTeff_Z->GetBinContent(binx, biny));
		}
	}	

	TCanvas *canZ = new TCanvas("canZ","",600,600);
	canZ->cd();
//	gPad->SetLogx();
//	gPad->SetLogy();
//	gStyle->SetPaintTextFormat("4.2f");
	Int_t PaletteColors[] = {9, kBlue, kBlue-4,kCyan, kTeal, kGreen,kSpring, 5, 2};
	gStyle->SetPalette(9, PaletteColors);
	//p_HLTeff_Z->Draw("E colz text");
	//p_crosseff->Draw("E colz text");
	//p_crosseff->Draw();
	p_L1_pho->SetMarkerStyle(20);
	p_L1_pho->Draw("EP");
	canZ->SaveAs("mgTrigger_efficiency_2017.pdf");


}

