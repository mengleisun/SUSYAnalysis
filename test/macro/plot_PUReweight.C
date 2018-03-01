#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TChain.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "../../include/tdrstyle.C"

void plot_PUReweight(){
	setTDRStyle();
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFillColor(0);
	
	TH1D *p_PU_data = new TH1D("p_PU_data",";N_{vtx};",100,0,100);
	TH1D *p_PU_MC = new TH1D("p_PU_MC","",100,0,100);
	TH1D *p_PU_MC_raw = new TH1D("p_PU_MC_raw","",100,0,100);

	TChain *sigtree = new TChain("signalTree");
	sigtree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_mgsignal_MuonEG_FullEcal.root");
	sigtree->Draw("nVertex >> p_PU_data");
 
	TChain *mctree;
  mctree = new TChain("mgTree","mgTree");
	mctree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_VGamma_DY.root");
  int   nVertex(0);
	float PUweight(1);
  mctree->SetBranchAddress("nVertex",   &nVertex);
	mctree->SetBranchAddress("PUweight",  &PUweight);
	for(unsigned ievt(0); ievt < mctree->GetEntries(); ievt++){
		mctree->GetEntry(ievt);
		p_PU_MC->Fill(nVertex,PUweight);
		p_PU_MC_raw->Fill(nVertex);
	}

	p_PU_data->Sumw2();
	p_PU_data->Scale(1.0/p_PU_data->Integral(1,100));
	p_PU_MC->Scale(1.0/p_PU_MC->Integral(1,100));
	p_PU_MC_raw->Scale(1.0/p_PU_MC_raw->Integral(1,100));
	
	gStyle->SetOptStat(0);
	TCanvas *c_PU = new TCanvas("PU", "PU",600,600);
	c_PU->cd();
	TPad *PU_pad1 = new TPad("PU_pad1", "PU_pad1", 0, 0.3, 1, 1.0);
	PU_pad1->SetBottomMargin(0);
	PU_pad1->Draw();  
	PU_pad1->cd(); 
	p_PU_data->SetMaximum(1.2*p_PU_data->GetBinContent(p_PU_data->GetMaximumBin()));
	p_PU_data->SetLineColor(1);
	p_PU_data->SetMarkerStyle(20);
	p_PU_data->GetXaxis()->SetTitle("N_{vtx}");
	p_PU_data->Draw("P");
	p_PU_MC->SetLineColor(kMagenta);
	p_PU_MC_raw->SetLineColor(kBlue);
	p_PU_MC_raw->Draw("hist same");
	p_PU_MC->Draw("hist same");
	TLegend *leg_PU =  new TLegend(0.5,0.65,0.9,0.9);
	leg_PU->SetFillStyle(0);
	leg_PU->AddEntry(p_PU_data,"Data");
	leg_PU->AddEntry(p_PU_MC_raw,"Simulation (no correction)");
	leg_PU->AddEntry(p_PU_MC,"Simulation (corrected)");
	leg_PU->Draw("same");
	prelim->Draw();
	lumitex->Draw();

	c_PU->cd();
	TPad *PU_pad2 = new TPad("PU_pad2", "PU_pad2", 0, 0.05, 1, 0.3);
	PU_pad2->SetTopMargin(0);
	PU_pad2->SetBottomMargin(0.3);
	PU_pad2->Draw();
	PU_pad2->cd();
  TLine *flatratio_PU = new TLine(0,1,100,1);
	TH1F *ratio_PU=(TH1F*)p_PU_data->Clone("transfer factor");
	ratio_PU->SetMarkerStyle(20);
	ratio_PU->SetLineColor(kBlack);
	ratio_PU->GetXaxis()->SetRangeUser(0,100);
	ratio_PU->GetYaxis()->SetRangeUser(0,2);
	ratio_PU->SetMinimum(0);
	ratio_PU->SetMaximum(2);
	ratio_PU->Divide(p_PU_MC);
	ratio_PU->SetTitle("");
	ratio_PU->GetYaxis()->SetTitle("Data/MC");
	ratio_PU->Draw();
	flatratio_PU->Draw("same");
	c_PU->SaveAs("PLOT_PUreweight.pdf");

}
