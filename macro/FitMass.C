#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include "TH1F.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TF1.h"
#include "TMath.h"
#include "TSystem.h"
#include "TChain.h"
#include "TTree.h"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooAbsReal.h"
#include "RooKeysPdf.h"
#include "RooConstVar.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooNumIntConfig.h"
#include "RooFFTConvPdf.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooHistPdf.h"
#include "RooMCStudy.h"
#include "RooChi2MCSModule.h"
#include "RooFitResult.h"
#include "RooExponential.h"
#include "RooCMSShape.h"
#include "RooDCBShape.h"
#include "RooUserPoly.h"
//#include "../include/RooCBExGaussShape.h"

void FitMass(){
  
	float fitrangehigh = 120;
	float fitrangelow  = 60;

//assume the histogram of invmass is named p_invmass
	
//************* Construct RooFit Models *********************************************************//
  RooRealVar mass_axis("invmass","invmass",fitrangelow, fitrangehigh);
  TCanvas *c_fitMass = new TCanvas("c_fitMass", "", 600, 600);
  RooDataHist datahist_data("both", "", mass_axis, p_invmass);

  c_fitMass->cd();
	double inilambda = (log(p_invmass->GetBinContent(2)) - log(p_invmass->GetBinContent(48)))/(p_invmass->GetBinCenter(2) - p_invmass->GetBinCenter(48)); 
	RooRealVar lambda("lambda", "slope", inilambda, -100, 100.);
	RooExponential *expo = new RooExponential("expo", "exponential PDF", mass_axis, lambda);

  RooRealVar m0( "m0", "m0", 91.188, 80,100);
  RooRealVar width( "width", "width", 2.495, 0, 15);
  RooRealVar mean("mean", "" ,0.,-1,1);
  RooRealVar sigma("sigma", "", 1.0, 0, 2);
  RooRealVar alpha("alpha", "", 1.0, 0.0, 20.0);
  RooRealVar n("n","", 1.0, 0.0, 20.0);
  RooRealVar alpha2("2ndalpha","", 1.0, 0.0, 20.0);
  RooRealVar n2("2ndn", "", 1.0, 0.0, 20.0);
  RooBreitWigner bw("bw", "", mass_axis, m0, width);
  RooDCBShape *cb = new RooDCBShape("cb","cb", mass_axis, mean, sigma, alpha, n, alpha2, n2);
  RooFFTConvPdf *signalRes = new RooFFTConvPdf("pdf", "pdf",mass_axis, bw, *cb);
  double iniSig = 0.8*p_invmass->Integral(1,100);
  double iniBkg = 60*(p_invmass->GetBinContent(10)+p_invmass->GetBinContent(50));
	if(iniBkg > p_invmass->GetEntries()/2)iniBkg = p_invmass->GetEntries()/2;
  RooRealVar nSig("nSig", "", iniSig, 0, p_invmass->GetEntries()*1.2);
  RooRealVar nBkg("nBkg", "", iniBkg, 0, p_invmass->GetEntries());
  RooAddPdf *model;
  model = new RooAddPdf("model", "", RooArgList(*expo, *signalRes),RooArgList(nBkg, nSig));   

  RooPlot* mass_Frame = mass_axis.frame(RooFit::Title(""),RooFit::Bins(100));
  model->fitTo(datahist_data);
  model->plotOn(mass_Frame, RooFit::Components(*expo),
				 RooFit::LineStyle(kDashed),
				 RooFit::Normalization(1.0, RooAbsReal::RelativeExpected));
  mass_Frame->Draw();

  mass_axis.setRange("signal",60,120);
  RooAbsReal* igx_sig;
	igx_sig = signalRes->createIntegral(mass_axis,RooFit::NormSet(mass_axis),RooFit::Range("signal"));
	double norminalmean = igx_sig->getVal()*(nSig.getVal());
	double norminalrms  = igx_sig->getVal()*(nSig.getError());

}
