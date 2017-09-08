#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include "TH1F.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TF1.h"
#include "TH2F.h"
#include "TMath.h"
#include "TSystem.h"
#include "TFractionFitter.h"
#include "TLatex.h"
#include "TGraphErrors.h"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooExponential.h"
#include "RooMCStudy.h"
#include "RooChi2MCSModule.h"
#include "RooPlot.h"
#include "TH1.h"
#include "RooFitResult.h"
#include "RooAbsReal.h"
#include "RooNumIntConfig.h"

#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_mcData.h"
#include "../../../include/analysis_tools.h"

void
plotFit(){

	gStyle->SetOptStat(0);
	TGraphErrors *p_frac = new TGraphErrors(4);
	TH1F *p_frac_0 = new TH1F("p_frac_0","",200,0,3);
	TH1F *p_frac_1 = new TH1F("p_frac_1","",200,0,3);
	TH1F *p_frac_2 = new TH1F("p_frac_2","",200,0,3);
	TH1F *p_frac_3 = new TH1F("p_frac_3","",200,0,3);

	std::ifstream vgammascalefile("VGamma_scalefactor_mg.txt");
	float fakescale, vgammascale, fakescaleerror, vgammascaleerror;
	float leplow, lephigh;

	float highest(0), lowest(1);
	float fittingerror(0), systematicerror(0), totalerror(0);

	for(unsigned i(0);  i < 4000; i++){
		vgammascalefile >> leplow >> lephigh >> fakescale >> fakescaleerror >> vgammascale >> vgammascaleerror;
		if(i == 0)p_frac->SetPoint(0, 25, fakescale);
		if(i == 0)fittingerror = fakescaleerror;
		p_frac_0->Fill(fakescale);
	}
	p_frac_0->Fit("gaus");
	systematicerror = p_frac_0->GetFunction("gaus")->GetParameter(2);
	if(p_frac_0->GetFunction("gaus")->GetParameter(1) + systematicerror > highest)highest=p_frac_0->GetFunction("gaus")->GetParameter(1) + systematicerror;
	if(p_frac_0->GetFunction("gaus")->GetParameter(1) - systematicerror < lowest)lowest=p_frac_0->GetFunction("gaus")->GetParameter(1) - systematicerror;
	totalerror = sqrt(systematicerror*systematicerror + fittingerror*fittingerror);
	p_frac->SetPointError(0,2.5,totalerror);

	for(unsigned i(0);  i < 1000; i++){
		vgammascalefile >> leplow >> lephigh >> fakescale >> fakescaleerror >> vgammascale >> vgammascaleerror;
		if(i == 0)p_frac->SetPoint(1, 40, vgammascale);
		if(i == 0)fittingerror = vgammascaleerror; 
		p_frac_1->Fill(vgammascale);
	}
	p_frac_1->Fit("gaus");
	systematicerror = p_frac_1->GetFunction("gaus")->GetParameter(2);
	if(p_frac_1->GetFunction("gaus")->GetParameter(1) + systematicerror > highest)highest=p_frac_1->GetFunction("gaus")->GetParameter(1) + systematicerror;
	if(p_frac_1->GetFunction("gaus")->GetParameter(1) - systematicerror < lowest)lowest=p_frac_1->GetFunction("gaus")->GetParameter(1) - systematicerror;
	totalerror = sqrt(systematicerror*systematicerror + fittingerror*fittingerror);
	p_frac->SetPointError(1,10,totalerror);


//	for(unsigned i(0);  i < 1000; i++){
//		vgammascalefile >> leplow >> lephigh >> fakescale >> fakescaleerror >> vgammascale >> vgammascaleerror;
//		if(i == 0)p_frac->SetPoint(2,60, fakescale);
//		if(i == 0)fittingerror = fakescaleerror; 
//		p_frac_2->Fill(fakescale);
//	}
//	p_frac_2->Fit("gaus");
//	systematicerror = p_frac_2->GetFunction("gaus")->GetParameter(2);
//	if(p_frac_2->GetFunction("gaus")->GetParameter(1) + systematicerror > highest)highest=p_frac_2->GetFunction("gaus")->GetParameter(1) + systematicerror;
//	if(p_frac_2->GetFunction("gaus")->GetParameter(1) - systematicerror < lowest)lowest=p_frac_2->GetFunction("gaus")->GetParameter(1) - systematicerror;
//	totalerror = sqrt(systematicerror*systematicerror + fittingerror*fittingerror);
//	p_frac->SetPointError(2,10,totalerror);
//
//	for(unsigned i(0);  i < 100; i++){
//		vgammascalefile >> leplow >> lephigh >> fakescale >> fakescaleerror >> vgammascale >> vgammascaleerror;
//		if(i == 0)p_frac->SetPoint(3,135, fakescale);
//		if(i == 0)fittingerror = fakescaleerror; 
//		p_frac_3->Fill(fakescale);
//	}
//	p_frac_3->Fit("gaus");
//	systematicerror = p_frac_3->GetFunction("gaus")->GetParameter(2);
//	if(p_frac_3->GetFunction("gaus")->GetParameter(1) + systematicerror > highest)highest=p_frac_3->GetFunction("gaus")->GetParameter(1) + systematicerror;
//	if(p_frac_3->GetFunction("gaus")->GetParameter(1) - systematicerror < lowest)lowest=p_frac_3->GetFunction("gaus")->GetParameter(1) - systematicerror;
//	totalerror = sqrt(systematicerror*systematicerror + fittingerror*fittingerror);
//	p_frac->SetPointError(3,65,totalerror);


//	TGraphErrors *p_syserror = new TGraphErrors(1);
//  p_syserror->SetPoint(0,100,(highest+lowest)/2);
//	p_syserror->SetPointError(0, 100, (highest-lowest)/2);
//
//	TCanvas *can = new TCanvas("scalefactor","",600,600);
//	TH2F *dummy = new TH2F("dummy","scalefactor for fake-lepton;muon Pt (GeV);f_{fakelepton}",4,0,200, 1,0,1);
//	//dummy->GetXaxis()->SetTitleSize(20);
//	dummy->Draw();
//	p_frac->SetMarkerStyle(20);
//	p_frac->Draw("EP same");
//	p_syserror->SetFillColor(2);
//  p_syserror->SetFillStyle(3244);
//  p_syserror->Draw("E2 same");	
//	p_frac->Draw("EP same");
//	can->SaveAs("fakescale_mg_ptdependence.pdf");

		
}
