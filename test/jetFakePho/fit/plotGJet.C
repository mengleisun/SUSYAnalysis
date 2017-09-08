#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<ctime>
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TF1.h"
#include "TF3.h"
#include "TH2F.h"
#include "TMath.h"
#include "TSystem.h"
#include "TFractionFitter.h"
#include "TLatex.h"
#include "TStyle.h"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooAbsReal.h"
#include "RooConstVar.h"
#include "RooDataHist.h"
#include "RooKeysPdf.h"
#include "RooPlot.h"
#include "RooHistPdf.h"
#include "RooNumIntConfig.h"
#include "RooFFTConvPdf.h"
#include "RooPolynomial.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooFitResult.h"

#include "../../../include/tdrstyle.C"
#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_mcData.h"
#include "../../../include/analysis_tools.h"
#include "../../../include/analysis_fakes.h"

int
plotGJet(){

	setTDRStyle();   
//************ Signal Tree **********************//
	TChain *mctree = new TChain("egTree");
	mctree->Add("/uscms_data/d3/mengleis/FullStatusData/plot_hadron_GJet.root");
	float mc_phoEt(0);
	float mc_phoEta(0); 
	float mc_phoPhi(0); 
	float mc_phoSigma(0);
	float mc_phoChIso(0);
	std::vector<int> *mc_mcPID=0;
	std::vector<float> *mc_mcEta=0;
	std::vector<float> *mc_mcPhi=0;
	std::vector<float> *mc_mcPt=0;
	std::vector<int> *mc_mcMomPID=0; 

	mctree->SetBranchAddress("phoEt",     &mc_phoEt);
	mctree->SetBranchAddress("phoEta",    &mc_phoEta);
	mctree->SetBranchAddress("phoPhi",    &mc_phoPhi);
	mctree->SetBranchAddress("phoSigma",  &mc_phoSigma);
	mctree->SetBranchAddress("phoChIso",  &mc_phoChIso);
	mctree->SetBranchAddress("mcPID",   &mc_mcPID);
	mctree->SetBranchAddress("mcEta",   &mc_mcEta);
	mctree->SetBranchAddress("mcPhi",   &mc_mcPhi);
	mctree->SetBranchAddress("mcPt",    &mc_mcPt);
	mctree->SetBranchAddress("mcMomPID",&mc_mcMomPID);

	TH1F *p_total = new TH1F("p_total",";Iso_{h^{#pm}};Events",20,0,20);
	TH1F *p_fake = new TH1F("p_fake",";Iso_{h^{#pm}};Events",20,0,20);
	TH1F *p_totalsigma = new TH1F("p_totalsigma",";#sigma_{i#eta i#eta};",180,0.0013,0.0193);
	TH1F *p_sigsigma = new TH1F("p_sigsigma",";#sigma_{i#eta i#eta};",180,0.0013,0.0193);
	TH1F *p_sbsigma = new TH1F("p_sbsigma",";#sigma_{i#eta i#eta};",180,0.0013,0.0193);
	for(unsigned ievt(0); ievt < mctree->GetEntries(); ievt++){
		mctree->GetEntry(ievt);
	
		bool isFake(true);
		unsigned mcIndex(0);
		unsigned mcPhoIndex(0);
		float mindR(0.7);
		float phodR(0.7);
		bool hasMatch(false);
		for(unsigned ii(0); ii < mc_mcPID->size(); ii++){
			float dR = DeltaR(mc_phoEta, mc_phoPhi, (*mc_mcEta)[ii], (*mc_mcPhi)[ii]);
			float dE = fabs(mc_phoEt - (*mc_mcPt)[ii])/mc_phoEt;
			if((*mc_mcPID)[ii] == 22){
				phodR = DeltaR(mc_phoEta, mc_phoPhi, (*mc_mcEta)[ii], (*mc_mcPhi)[ii]);
				mcPhoIndex = ii;
			}
			if(dR < mindR && dE < 0.7){mindR = dR; mcIndex = ii; hasMatch = true;}
		}
		if(phodR < 0.1)mcIndex = mcPhoIndex;
		if(hasMatch)isFake = isHad(fabs((*mc_mcPID)[mcIndex]), fabs((*mc_mcMomPID)[mcIndex]));

		if(!isFake){
			p_totalsigma->Fill(mc_phoSigma);
			if(mc_phoSigma <= 0.0103)p_sigsigma->Fill(mc_phoSigma);
			else if(mc_phoSigma > 0.0103 && mc_phoSigma < 0.016)p_sbsigma->Fill(mc_phoSigma);
		}

		if(mc_phoSigma > 0.0103)continue;
		p_total->Fill(mc_phoChIso);
		if(isFake)p_fake->Fill(mc_phoChIso);
	}

	TCanvas *gjet = new TCanvas("gjet","",600,600);
	gjet->cd();
	gPad->SetLogy();	
	p_total->GetXaxis()->SetTitle("Iso_{h^{#pm}} (GeV)");
	p_total->GetXaxis()->SetTitleOffset(1);
	p_total->GetXaxis()->SetTitleSize(20);
	p_total->SetLineWidth(2);
	p_fake->SetLineWidth(2);
	p_total->Draw();
	p_fake->SetLineColor(kRed);
	p_fake->Draw("same");
	TLegend *leg= new TLegend(0.5,0.7,0.85,0.85);
	leg->AddEntry(p_total, "selected photons");
	leg->AddEntry(p_fake, "jet->#gamma fakes");
	leg->Draw("same");
	gjet->SaveAs("fake_in_GJet.pdf");

	TLatex* latex = new TLatex();
	TCanvas *sigma = new TCanvas("sigma","",600,600);
	sigma->cd();
	gPad->SetLogy();	
	p_totalsigma->GetXaxis()->SetTitleOffset(1);
	p_totalsigma->GetXaxis()->SetTitleSize(20);
	p_totalsigma->Scale(1/p_totalsigma->GetEntries());
	p_sigsigma->Scale(1/p_totalsigma->GetEntries());
	p_sbsigma->Scale(1/p_totalsigma->GetEntries());
	p_totalsigma->Draw();
	p_sigsigma->SetFillColor(kBlue);
	p_sbsigma->SetFillColor(kYellow);
	p_sigsigma->Draw("same");
	p_sbsigma->Draw("same");
	latex->DrawLatex(0.006,0.000005,"p_{signal}");
	latex->DrawLatex(0.011,0.000005,"p_{sideband}");
	sigma->SaveAs("iteration.pdf");

	return 1;
}
