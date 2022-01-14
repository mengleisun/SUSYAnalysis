#include "TChain.h"
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include "TSystem.h"
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TF1.h"
#include "TPad.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_rawData.h"
#include "../../../include/tdrstyle.C"
#include "../../../include/analysis_tools.h"

#define NBIN 17

void plotJetFake(){//main 

	setTDRStyle();
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	gStyle->SetErrorX(0.5);
	gStyle->SetTitleX(0.5);

	bool useMC(false);
	bool doCompare(false);
	//std::ifstream jetfake_file("../result/JetFakeRate-ChIso-DoubleEG-ReMiniAOD.txt");
	std::ifstream jetfake_file("JetFakeRate-DoubleEG-EB.txt");
	 
	gSystem->Load("/uscms/home/tmishra/work/CMSSW_10_2_22/src/SUSYAnalysis/lib/libAnaClasses.so");
	gStyle->SetOptStat(0);
	TGraphErrors *jetfrac = new TGraphErrors(NBIN);
	TGraphErrors *jetfracsys = new TGraphErrors(NBIN);
	TGraphErrors *jetfracTrue= new TGraphErrors(NBIN);
	TGraphErrors *ratio= new TGraphErrors(NBIN);
	TGraphErrors *jetfrac_nocorr = new TGraphErrors(NBIN);
	TGraphErrors *ratio_nocorr= new TGraphErrors(NBIN);	

	int i(0);
	float pt_lower(0), pt_upper(0);
	float fakerate(0), error(0), systematic(0); 
	float truerate(0);

	if(jetfake_file.is_open()){
		for(int i(0); i < NBIN; i++){ 
			if(!useMC)jetfake_file >> pt_lower >> pt_upper >> fakerate >> error >> systematic >> truerate;
			if(useMC)jetfake_file >>  pt_lower >> pt_upper >> fakerate >> error >> systematic >> truerate;
			jetfrac->SetPoint(i, (pt_lower+pt_upper)/2.0, fakerate);
                        cout<<fakerate<<endl;
			jetfrac->SetPointError(i, (pt_upper-pt_lower)/2.0, error);
			jetfracsys->SetPoint(i, (pt_lower+pt_upper)/2.0, fakerate);
			jetfracsys->SetPointError(i, 0, systematic);
			if(useMC)jetfracTrue->SetPoint(i, (pt_lower+pt_upper)/2.0, truerate);
			if(useMC)ratio->SetPoint(i, (pt_lower+pt_upper)/2.0, fakerate/truerate);
			if(useMC)ratio->SetPointError(i, (pt_upper-pt_lower)/2.0, error/truerate);
		}
		jetfake_file.close(); 
	}
	
	if(doCompare){
		std::ifstream jetfake_nocorrfile("../result/JetFakeRate-ChIso-GJet-nocorr.txt");
		if(jetfake_nocorrfile.is_open()){
			for(int i(0); i < NBIN; i++){
				if(!useMC)jetfake_nocorrfile >> pt_lower >> pt_upper >> fakerate >> error >> systematic;
				if(useMC)jetfake_nocorrfile >>  pt_lower >> pt_upper >> fakerate >> error >> systematic >> truerate;
				jetfrac_nocorr->SetPoint(i, (pt_lower+pt_upper)/2.0, fakerate);
				jetfrac_nocorr->SetPointError(i, (pt_upper-pt_lower)/2.0, error);
				if(useMC)ratio_nocorr->SetPoint(i, (pt_lower+pt_upper)/2.0, fakerate/truerate);
				if(useMC)ratio_nocorr->SetPointError(i, (pt_upper-pt_lower)/2.0, error/truerate);
			}
			jetfake_nocorrfile.close();
		}
	}

	jetfrac->SetMarkerStyle(20);
	jetfrac->SetMarkerColor(kBlack);
	jetfrac->SetLineColor(kBlack);
	jetfracsys->SetMarkerStyle(20);
	jetfracsys->SetMarkerColor(kBlack);
	jetfracsys->SetLineColor(kRed);
	jetfrac->SetLineWidth(2);
	jetfracsys->SetLineWidth(2);
	jetfrac->SetFillColor(0);
	jetfrac_nocorr->SetFillColor(0);
	jetfracTrue->SetFillColor(0);
	jetfracTrue->SetLineColor(kRed);
	jetfrac_nocorr->SetMarkerStyle(20);
	jetfrac_nocorr->SetMarkerColor(kBlue);
	jetfrac_nocorr->SetLineColor(kBlue);

	TCanvas *canvas = new TCanvas("Hadron Fraction"," fake rate",600,600);
	canvas->cd();
//	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
//	pad1->SetBottomMargin(0);
//	pad1->Draw();  
//	pad1->cd();  
	TH1F *dummy = new TH1F("Hadron Fraction","e#gamma channel;p_{T}(GeV);hadron fraction",17,30,200);
	dummy->GetXaxis()->SetTitle("p_{T} (GeV)");
	dummy->GetXaxis()->SetTitleOffset(1);
	dummy->SetMaximum(0.6);
	dummy->Draw();
	jetfrac->Draw("P same");
//	jetfracsys->Draw("P same");
	jetfracTrue->SetFillStyle(0);
	if(doCompare)jetfrac_nocorr->SetFillStyle(0);
	if(doCompare)jetfrac_nocorr->Draw("P same");
	if(useMC)jetfracTrue->Draw("P same");
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFillColor(0);
//	TLegend *leg =  new TLegend(0.4,0.75,0.85,0.9);
//	jetfrac->SetFillStyle(0);
//	leg->SetFillStyle(0);
//	leg->AddEntry(jetfrac, "GJet corrected fake rate");
//	if(doCompare)leg->AddEntry(jetfrac_nocorr, "GJet uncorrected fake rate");
//	if(useMC)leg->AddEntry(jetfracTrue,"MC truth");
//	leg->Draw("same");

//	canvas->cd();
//	TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
//	pad2->Draw();
//	pad2->cd();
//	pad2->SetTopMargin(0);
//	pad2->SetBottomMargin(0.3);
//	TH1F *dummy2 = new TH1F("dummy2",";p_{T}(GeV);estimated/true",17,30,200);
//	dummy2->SetMinimum(0);
//	dummy2->SetMaximum(5);
//	dummy2->SetMarkerStyle(20);
//	dummy2->Draw();
//  TLine *flatratio = new TLine(30,1,200,1);
//	ratio->Draw("P same");
//	if(doCompare)ratio_nocorr->SetMarkerStyle(20);
//	if(doCompare)ratio_nocorr->SetMarkerColor(kBlue);
//	if(doCompare)ratio_nocorr->Draw("P same");
//	flatratio->Draw("same");
	canvas->SaveAs("JetFakePho_DoubleEG_ReMiniAOD.pdf");
}

