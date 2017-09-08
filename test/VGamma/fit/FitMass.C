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
#include "TRandom3.h"
#include "TPad.h"
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

int
FitMass(int ih, int metlow, int methigh, int leplow, int lephigh, int isocut){
	gStyle->SetOptStat(0);
  RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooIntegrator1D").setRealValue("maxSteps",50000); 
	std::ostringstream histname;
	ofstream myfile;
	myfile.open("VGamma_scalefactor_mg.txt", std::ios_base::app | std::ios_base::out);

	std::ostringstream filename;
	filename.str("");
	filename << "../../Background/controlTree_mg_signal_ReMiniAOD" << "met" << metlow << "_" << methigh << "_pt" << leplow << "_" << lephigh << ".root";
	TFile *file_sig = TFile::Open(filename.str().c_str());
	filename.str("");
	filename << "../../Background/controlTree_mg_eleBkg_ReMiniAOD" << "met" << metlow << "_" << methigh << "_pt" << leplow << "_" << lephigh << ".root";
	TFile *file_ele = TFile::Open(filename.str().c_str());
	filename.str("");
	filename << "../../Background/controlTree_mg_jetbkg_ReMiniAOD" << "met" << metlow << "_" << methigh << "_pt" << leplow << "_" << lephigh << ".root";
	TFile *file_jet = TFile::Open(filename.str().c_str());
	filename.str("");
	filename << "../../Background/controlTree_mg_rareBkg_ReMiniAOD" << "met" << metlow << "_" << methigh << "_pt" << leplow << "_" << lephigh << ".root";
	TFile *file_rare = TFile::Open(filename.str().c_str());
	filename.str("");
	filename << "../../Background/controlTree_mg_qcd_ReMiniAOD" << "met" << metlow << "_" << methigh << "_pt" << leplow << "_" << lephigh << "_iso" << isocut << ".root";
	TFile *file_qcd = TFile::Open(filename.str().c_str());
	filename.str("");
	filename << "../../Background/controlTree_mg_VGBkg_ReMiniAOD" << "met" << metlow << "_" << methigh << "_pt" << leplow << "_" << lephigh << ".root";
	TFile *file_VG = TFile::Open(filename.str().c_str());

	TH1F  *p_sig = (TH1F*)file_sig->Get("p_Mt");
	TH1F  *p_rare = (TH1F*)file_rare->Get("p_Mt");
	histname.str("");
	if(ih == 0)histname << "p_Mt";
	else histname << "toy_eledPhiEleMET_" << ih;
	TH1F  *p_ele = (TH1F*)file_ele->Get(histname.str().c_str());
	TH1F  *p_jet = (TH1F*)file_jet->Get(histname.str().c_str());

	TH1F *p_rawsig = (TH1F*)p_sig->Clone("p_rawsig");
	TH1F *p_target = (TH1F*)p_sig->Clone("fit_target");
	p_target->Add(p_rare, -1);
	p_target->Add(p_ele, -1);
	p_target->Add(p_jet, -1);
	p_target->Sumw2();
  TH1F *p_proxy = (TH1F*)file_qcd->Get("p_Mt");
  TH1F *p_MC = (TH1F*)file_VG->Get("p_Mt"); 
	gRandom = new TRandom3(0);
	gRandom->SetSeed(0);
	double radomMC = -1+ gRandom->Rndm()*2.0;
	if(ih == 0)radomMC = 0;
	for(int ibin(1); ibin < p_MC->GetSize(); ibin++)p_MC->SetBinContent(ibin, p_MC->GetBinContent(ibin)+radomMC*p_MC->GetBinError(ibin));
	p_MC->Sumw2();

  p_proxy->GetXaxis()->SetTitle("#Delta#phi");
  p_target->GetXaxis()->SetTitle("#Delta#phi");

	p_rawsig->Rebin(4);
	p_target->Rebin(4);
	p_proxy->Rebin(4);
	p_MC->Rebin(4);	

  RooRealVar dphi("dphi","",0,150);
  RooDataHist* h_target = new RooDataHist("h_target","h_target",dphi,p_target);
  RooDataHist* h_proxy  = new RooDataHist("h_proxy", "h_proxy", dphi,p_proxy );
  RooDataHist* h_MC     = new RooDataHist("h_MC",    "h_MC",    dphi,p_MC);
  RooHistPdf*  pdf_proxy= new RooHistPdf("pdf_proxy","pdf_proxy",dphi, *h_proxy);
  RooHistPdf*  pdf_MC   = new RooHistPdf("pdf_MC",   "pdf_MC",   dphi, *h_MC);
  RooRealVar nFake("nFake","nFake", 0.5*p_target->Integral(1,25), 0, p_target->Integral(1,25)*2);
  RooRealVar nMC("nMC","nMC",0.5*p_target->Integral(1,25),0, p_target->Integral(1,25)*2);
  RooRealVar fakefrac("fakefrac","fakefrac",0.5,0,1.0);
  RooAddPdf model2("model2","",RooArgList(*pdf_proxy, *pdf_MC),RooArgList(nFake,nMC));
  RooAddPdf model("model","",RooArgList(*pdf_proxy, *pdf_MC),fakefrac);

  RooPlot* frame = dphi.frame(RooFit::Title("#Delta#phi(l,MET) (40 GeV < MET < 70 GeV)"));
  TCanvas *can=new TCanvas("can","",600,600);
  can->cd(1);
//	model2.fitTo(*h_target,RooFit::SumW2Error(kTRUE),RooFit::Save());
  RooFitResult* result = model.fitTo(*h_target,RooFit::SumW2Error(kTRUE),RooFit::Save());
  h_target->plotOn(frame);
  model.plotOn(frame,
		   RooFit::FillColor(kBlue-4));
  model.plotOn(frame, RooFit::Components(*pdf_proxy),
               RooFit::LineStyle(kDashed),
               RooFit::LineColor(kRed),
               RooFit::Normalization(1.0));
  model.plotOn(frame, RooFit::Components(*pdf_MC),
               RooFit::LineStyle(kDashed),
               RooFit::LineColor(kGreen),
               RooFit::Normalization(1.0));
  frame->Draw();
	std::ostringstream figurename;
  figurename.str("");
  figurename << "fit_MT_mg_met" << metlow << "_" << methigh << "_pt" <<  leplow << "_" << lephigh <<  "_iso" << isocut << ".pdf";
  can->SaveAs(figurename.str().c_str());


  return 1;
}
