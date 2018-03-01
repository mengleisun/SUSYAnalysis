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

void plotMuMugammaR9(){//main  
  gStyle->SetOptStat(0);
	gStyle->SetPaintTextFormat("4.4f");
	setTDRStyle();   
  Double_t xaxis2d[] = {35,40,45,60,200};
  Double_t yaxis2d[] = {25,38,45,60,200};
	TProfile2D *p_HLTeff_Z  = new TProfile2D("p_HLTeff_Z", "#mu#gamma trigger efficiency;#gamma p_{T} (GeV); #mu p_{T} (GeV)",4,xaxis2d,4,-1.4442,1.4442);
	TProfile2D *p_HLTeff_mg = new TProfile2D("p_HLTeff_mg","MET dataset",4,xaxis2d,4,-1.4442,1.4442);
	TProfile2D *p_HLTeff_DY = new TProfile2D("p_HLTeff_DY","DY efficiency",4,xaxis2d,4,-1.4442,1.4442);
	TH2F 			 *p_crosseff = new TH2F("p_mgEff","MuonEG trigger efficiency; #gamma p_{T} (GeV); muon p_{T} (GeV)",4,xaxis2d,4,-1.4442,1.4442);
	TH2F 			 *p_mgESF = new TH2F("p_mgESF","MuonEG trigger ESF; #gamma p_{T} (GeV); muon p_{T} (GeV)",4,xaxis2d,4,-1.4442,1.4442);


	TChain *Ztree = new TChain("ZTree","ZTree");
	Ztree->Add("/uscms_data/d3/mengleis/FullStatusOct/plot_MuonTrigger_ReMiniAOD.root");
	float Z_phoEt(0);
	float Z_muPt(0);
	float Z_phoEta(0);
	float Z_dR(0);
	float Z_R9(0);

	Ztree->SetBranchAddress("phoEt",     &Z_phoEt);
	Ztree->SetBranchAddress("phoEta",    &Z_phoEta);
	Ztree->SetBranchAddress("muPt",      &Z_muPt);
	Ztree->SetBranchAddress("dR",        &Z_dR);
	Ztree->SetBranchAddress("phoR9",        &Z_R9);

	for(unsigned ievt(0); ievt < Ztree->GetEntries(); ievt++){
		Ztree->GetEntry(ievt);
		if(fabs(Z_phoEta) > 1.4442)continue;
		if(Z_phoEt > 200)Z_phoEt = 199;
		if(Z_muPt > 200)Z_muPt = 199;
		
		if( Z_R9 > 0.5)
			p_HLTeff_Z->Fill(Z_phoEt, Z_phoEta, 1);
		else p_HLTeff_Z->Fill(Z_phoEt, Z_phoEta, 0);
	}


	TCanvas *canmg = new TCanvas("canmg","",600,600);
	canmg->cd();
	p_HLTeff_mg->Draw("colz text");
  for(int binx(1); binx <= 4; binx++){
		for(int biny(1); biny <= 4; biny++){
			p_crosseff->SetBinContent(binx, biny, p_HLTeff_Z->GetBinContent(binx, biny));
			p_crosseff->SetBinError(binx, biny, 0); 
		}
	}	

	TChain *DYtree = new TChain("mgTree","mgTree");
	//DYtree->Add("/uscms_data/d3/mengleis/FullStatusOct/plot_MuonTrigger_DY.root");
	DYtree->Add("/uscms_data/d3/mengleis/plot_MuonTrigger_DY.root");
	float DY_phoEt(0);
	float DY_phoEta(0);
	float DY_muPt(0);
	float DY_muMiniIso(0);
	float DY_R9(0);

	DYtree->SetBranchAddress("phoEt",     &DY_phoEt);
	DYtree->SetBranchAddress("phoEta",    &DY_phoEta);
	DYtree->SetBranchAddress("muPt",      &DY_muPt);
	DYtree->SetBranchAddress("muMiniIso", &DY_muMiniIso);
	DYtree->SetBranchAddress("R9",        &DY_R9);

	for(unsigned ievt(0); ievt < DYtree->GetEntries(); ievt++){
		DYtree->GetEntry(ievt);
		if(DY_phoEt > 200)DY_phoEt = 199;
		if(DY_muPt > 200)DY_muPt = 199;
		if(DY_R9 > 0.5)
			p_HLTeff_DY->Fill(DY_phoEt, DY_phoEta, 1);
		else p_HLTeff_DY->Fill(DY_phoEt, DY_phoEta, 0);
	}
  for(int binx(1); binx <= 4; binx++){
		for(int biny(1); biny <= 4; biny++){
			p_mgESF->SetBinContent(binx, biny, p_crosseff->GetBinContent(binx, biny)/p_HLTeff_DY->GetBinContent(binx, biny));
			p_mgESF->SetBinError(binx, biny, p_crosseff->GetBinError(binx, biny)/p_HLTeff_DY->GetBinContent(binx, biny));
		}
	}	
	

	TCanvas *canESF = new TCanvas("canESF","",600,600);
	gPad->SetLogx();
	gPad->SetLogy();
	gStyle->SetPaintTextFormat("4.4f");
	p_mgESF->Draw("E colz text");
	//canESF->SaveAs("mgTrigger_ESF.pdf");
	
}

