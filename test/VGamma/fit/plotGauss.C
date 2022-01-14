// Distribution of the scale factors derived from 1000 toy MC experiments. plot 29, AN
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
plotGauss(){

	gStyle->SetOptStat(0);
	std::ifstream vgammascalefile("VGamma_scalefactor_mg_test.txt");
	float fakescale, vgammascale, fakescaleerror, vgammascaleerror;
	float leplow, lephigh;
	TH1F *p_frac_0 = new TH1F("p_frac_0","",100,1,1.5);

	float highest(0), lowest(1);
	float fittingerror(0), systematicerror(0), totalerror(0);

	double lowbound(2), highbound(0), norm(0);
	for(unsigned i(0);  i < 1000; i++){
		vgammascalefile >> leplow >> lephigh >> fakescale >> fakescaleerror >> vgammascale >> vgammascaleerror;
		if(i == 0)norm = vgammascale;
		if(i == 0)fittingerror = vgammascaleerror;
		p_frac_0->Fill(vgammascale);
		if(vgammascale < lowbound)lowbound = vgammascale;
		if(vgammascale > highbound)highbound = vgammascale;
	}
	p_frac_0->Fit("gaus");
	//systematicerror = 3*p_frac_0->GetFunction("gaus")->GetParameter(2);
	systematicerror = (highbound-norm) > (norm - lowbound)? (highbound-norm):(norm - lowbound); 
	if(p_frac_0->GetFunction("gaus")->GetParameter(1) + systematicerror > highest)highest=p_frac_0->GetFunction("gaus")->GetParameter(1) + systematicerror;
	if(p_frac_0->GetFunction("gaus")->GetParameter(1) - systematicerror < lowest)lowest=p_frac_0->GetFunction("gaus")->GetParameter(1) - systematicerror;
	totalerror = sqrt(systematicerror*systematicerror + fittingerror*fittingerror);

	std::cout << "highbound = " << highbound << std::endl;
	std::cout << "lowbound = " << lowbound << std::endl;
	std::cout << "error = " << totalerror << std::endl;


	TCanvas *canscale = new TCanvas("canscale","",600,600);
	canscale->cd();
	p_frac_0->Draw();	
	canscale->SaveAs("VGammaScale_test.pdf");	
}
