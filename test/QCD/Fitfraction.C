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

#include "../../include/analysis_rawData.h"
#include "../../include/analysis_photon.h"
#include "../../include/analysis_muon.h"
#include "../../include/analysis_ele.h"
#include "../../include/analysis_mcData.h"
#include "../../include/analysis_tools.h"

int
Fitfraction(int ih, int leplow, int lephigh, int isocut){
	gStyle->SetOptStat(0);
  RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooIntegrator1D").setRealValue("maxSteps",50000); 
	std::ostringstream histname;
	ofstream myfile;
	myfile.open("VGamma_scalefactor_eg.txt", std::ios_base::app | std::ios_base::out);

	//TFile *file_sig = TFile::Open("../../Background/control/controlTree_eg_signal_ReMiniAOD.root");
	//TFile *file_ele = TFile::Open("../../Background/control/controlTree_eg_eleBkg_ReMiniAOD.root");
	//TFile *file_jet = TFile::Open("../../Background/control/controlTree_eg_jetbkg_ReMiniAOD.root");
	//TFile *file_rare = TFile::Open("../../Background/control/controlTree_eg_rareBkg_ReMiniAOD.root");
	//TFile *file_qcd = TFile::Open("../../Background/control/controlTree_eg_qcd_ReMiniAOD.root");
	//TFile *file_VG = TFile::Open("../../Background/control/controlTree_eg_VGBkg_ReMiniAOD.root");
	std::ostringstream filename;
	filename.str("");
	filename << "controlTree_eg_signal_ReMiniAOD" << leplow << "_" << lephigh << ".root";
	TFile *file_sig = TFile::Open(filename.str().c_str());
	filename.str("");
	filename << "controlTree_eg_qcd_ReMiniAOD_" << leplow << "_" << lephigh << "_met40" << "_iso" << isocut << ".root";
	TFile *file_qcd = TFile::Open(filename.str().c_str());
	filename.str("");
	filename << "controlTree_eg_VGBkg_ReMiniAOD" << leplow << "_" << lephigh << ".root";
	TFile *file_VG = TFile::Open(filename.str().c_str());
	//TFile *file_sig = TFile::Open("../../Background/control/controlTree_eg_signal_ReMiniAOD.root");
	//TFile *file_ele = TFile::Open("../../Background/control/controlTree_eg_eleBkg_ReMiniAOD.root");
	//TFile *file_jet = TFile::Open("../../Background/control/controlTree_eg_jetbkg_ReMiniAOD.root");
	//TFile *file_rare = TFile::Open("../../Background/control/controlTree_eg_rareBkg_ReMiniAOD.root");
	//TFile *file_qcd = TFile::Open("../../Background/control/controlTree_eg_qcd_ReMiniAOD.root");
	//TFile *file_VG = TFile::Open("../../Background/control/controlTree_eg_VGBkg_ReMiniAOD.root");

	TH1F  *p_sig = (TH1F*)file_sig->Get("p_dPhiEleMET");
	TH1F *p_target = (TH1F*)p_sig->Clone("fit_target");
	p_target->Sumw2();
  TH1F *p_proxy = (TH1F*)file_qcd->Get("p_dPhiEleMET");
  TH1F *p_MC = (TH1F*)file_VG->Get("p_dPhiEleMET"); 
	p_MC->Sumw2();

  p_proxy->GetXaxis()->SetTitle("#Delta#phi");
  p_target->GetXaxis()->SetTitle("#Delta#phi");

	p_target->Rebin(2);
	p_proxy->Rebin(2);
	p_MC->Rebin(2);	


  RooRealVar dphi("dphi","",0,3.2);
  RooDataHist* h_target = new RooDataHist("h_target","h_target",dphi,p_target);
  RooDataHist* h_proxy  = new RooDataHist("h_proxy", "h_proxy", dphi,p_proxy );
  RooDataHist* h_MC     = new RooDataHist("h_MC",    "h_MC",    dphi,p_MC);
  RooHistPdf*  pdf_proxy= new RooHistPdf("pdf_proxy","pdf_proxy",dphi, *h_proxy);
  RooHistPdf*  pdf_MC   = new RooHistPdf("pdf_MC",   "pdf_MC",   dphi, *h_MC);
  RooRealVar nFake("nFake","nFake", 20000, 0, p_target->GetEntries());
  RooRealVar nMC("nMC","nMC",20000,0, p_target->GetEntries());
  RooRealVar fakefrac("fakefrac","fakefrac",0.2,0,1.0);
  //RooAddPdf model("model","",RooArgList(*pdf_proxy, *pdf_MC),RooArgList(nFake,nMC));
  RooAddPdf model("model","",RooArgList(*pdf_proxy, *pdf_MC),fakefrac);

  RooPlot* frame = dphi.frame(RooFit::Title("#Delta#phi(l,MET) (40 GeV < MET < 70 GeV)"));
  TCanvas *can=new TCanvas("can","",600,600);
  can->cd(1);
  RooFitResult* result = model.fitTo(*h_target,RooFit::SumW2Error(kTRUE),RooFit::Save());
  h_target->plotOn(frame);
  model.plotOn(frame,
		   RooFit::FillColor(kBlue-4));
  model.plotOn(frame, RooFit::Components(*pdf_proxy),
               RooFit::LineStyle(kDashed),
               RooFit::LineColor(kRed),
               RooFit::Normalization(1.0));
  frame->Draw();

	TH1F *p_combine = new TH1F("p_combine","",16,0,3.2);
	TH1F *p_qcd  = new TH1F("p_qcd","",16,0,3.2);
	TH1F *norm_MC = (TH1F*)p_MC->Clone("norm_MC");
	TH1F *norm_proxy = (TH1F*)p_proxy->Clone("norm_proxy");
	norm_MC->Scale(1.0/norm_MC->Integral(1, norm_MC->GetSize()) );
	norm_proxy->Scale(1.0/norm_proxy->Integral(1, norm_proxy->GetSize()));
	p_combine->Add(norm_proxy, fakefrac.getVal());
	p_combine->Add(norm_MC, 1-fakefrac.getVal());
	p_combine->Scale(p_target->Integral(1, p_target->GetSize()));
	p_qcd->Add(norm_proxy, fakefrac.getVal()*p_target->Integral(1, p_target->GetSize()));
	TH1F *fitratio = new TH1F("fitratio",";#Delta#phi(l, E_{T}^{miss});fit/data", 16,0,3.2);
	TGraphErrors *fitratio_error = new TGraphErrors(16);
  double staterror = fakefrac.getError()/fakefrac.getVal();
	for(int ibin(1); ibin <= 16; ibin++)fitratio->SetBinContent(ibin, p_combine->GetBinContent(ibin)/p_target->GetBinContent(ibin));
	for(int ibin(1); ibin <= 16; ibin++)fitratio->SetBinError(ibin, p_target->GetBinError(ibin)/p_target->GetBinContent(ibin));
	for(int ibin(0); ibin < 16; ibin++)fitratio_error->SetPoint(ibin, -3.2+ibin*0.1, 1);
	for(int ibin(0); ibin < 16; ibin++)fitratio_error->SetPointError(ibin, 0.05, staterror);


//	RooDataHist* h_combine = new RooDataHist("h_combine","h_combine",dphi, p_combine);
//	h_combine->plotOn(frame,
//       RooFit::MarkerColor(kGreen-4));
  TCanvas *canres=new TCanvas("canres","",600,600);
  canres->cd(1);
	TPad *canpt_pad1 = new TPad("canpt_pad1", "pad1", 0, 0.3, 1, 1.0);
	canpt_pad1->SetBottomMargin(0);
	canpt_pad1->Draw();          
	canpt_pad1->cd();   
	p_target->SetTitle("#Delta#phi(e, E_{T}^{miss})"); 
	//p_target->SetMaximum(500);
	p_target->SetMinimum(0);
	p_target->SetMarkerStyle(20);
	p_target->Draw("EP");
	p_combine->SetLineColor(kBlue);
	p_combine->Draw("hist same");
	p_qcd->SetLineColor(kRed);
	p_qcd->Draw("hist same");
	TLegend *leg=new TLegend(0.6,0.7,0.85,0.85);
	leg->AddEntry(p_target, "observed");
	leg->AddEntry(p_combine, "fit result");
	leg->AddEntry(p_qcd, "fake leptons");
//	leg->Draw("same");
	 
	canres->cd();   
	TPad *canpt_pad2 = new TPad("canpt_pad2", "pad2", 0, 0.05, 1, 0.3);
	canpt_pad2->SetTopMargin(0);
	canpt_pad2->SetBottomMargin(0.3);
	canpt_pad2->Draw();
	canpt_pad2->cd(); 	
	TH1F *dummy_ptratio = new TH1F("dummy_ptratio",";#Delta#phi(l, E_{T}^{miss});fit/data",16,0,3.2);
	dummy_ptratio->SetMaximum(2);
	dummy_ptratio->SetMinimum(0);
	dummy_ptratio->Draw();
	fitratio->SetMarkerStyle(16);
	fitratio->Draw("EP same");
	fitratio_error->SetFillStyle(3005);
	fitratio_error->Draw("E2 same");	

	std::ostringstream plotname;
	plotname << "fit_dPhi_eg_" << leplow << "_" << lephigh << "_met40" << "_iso" << isocut << ".pdf";
  if(ih == 0)canres->SaveAs(plotname.str().c_str());

  double nMCtotal(0);
  for(unsigned i(1); i<=16; i++)nMCtotal+= p_MC->GetBinContent(i);
  std::cout << "factor for fake: " << fakefrac.getVal()*p_target->Integral(1, p_target->GetSize()) << "/" << p_proxy->Integral(1, p_proxy->GetSize()) << std::endl;
  std::cout << "factor for MC: " << (1-fakefrac.getVal())*p_target->Integral(1, p_target->GetSize()) << "/" << p_MC->Integral(1, p_MC->GetSize())  << std::endl;
  std::cout << "fakefrac=" << fakefrac.getVal()<< std::endl;
	myfile << leplow << " " << lephigh << " " << fakefrac.getVal()*p_target->Integral(1, p_target->GetSize())/p_proxy->Integral(1, p_proxy->GetSize()) << " " << (1-fakefrac.getVal())*p_target->Integral(1, p_target->GetSize())/p_MC->Integral(1, p_MC->GetSize()) << " " << fakefrac.getError()*p_target->Integral(1, p_target->GetSize())/p_MC->Integral(1, p_MC->GetSize()) << std::endl;
	myfile.close();
  return 1;
}
