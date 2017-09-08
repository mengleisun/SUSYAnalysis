#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include "TH1D.h"
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
FitfractionWGZG(int ih,int metlow, int methigh, int leplow, int lephigh, int isocut){
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
	filename << "../../Background/controlTree_mg_VGBkg_wgt" << "met" << metlow << "_" << methigh << "_pt" << leplow << "_" << lephigh << ".root";
	TFile *file_VG = TFile::Open(filename.str().c_str());

	TH1D  *p_sig = (TH1D*)file_sig->Get("p_dPhiEleMET");
	TH1D  *p_rare = (TH1D*)file_rare->Get("p_dPhiEleMET");
	histname.str("");
	if(ih == 0)histname << "p_dPhiEleMET";
	else histname << "toy_eledPhiEleMET_" << ih;
	TH1D  *p_ele = (TH1D*)file_ele->Get(histname.str().c_str());
	TH1D  *p_jet = (TH1D*)file_jet->Get(histname.str().c_str());

	TH1D *p_rawsig = (TH1D*)p_sig->Clone("p_rawsig");
	TH1D *p_target = (TH1D*)p_sig->Clone("fit_target");
	p_target->Add(p_rare, -1);
	p_target->Add(p_ele, -1);
	p_target->Add(p_jet, -1);
	p_target->Sumw2();
  TH1D *p_proxy = (TH1D*)file_qcd->Get("p_dPhiEleMET");
  TH1D *p_MC_WG = (TH1D*)file_VG->Get("p_dPhiEleMET_WG"); 
  TH1D *p_MC_ZG = (TH1D*)file_VG->Get("p_dPhiEleMET_ZG"); 
	gRandom = new TRandom3(0);
	gRandom->SetSeed(0);
	double radomMC = -1+ gRandom->Rndm()*2.0;
	if(ih == 0)radomMC = 0;
	for(int ibin(1); ibin < p_MC_WG->GetSize(); ibin++)p_MC_WG->SetBinContent(ibin, p_MC_WG->GetBinContent(ibin)+radomMC*p_MC_WG->GetBinError(ibin));
	for(int ibin(1); ibin < p_MC_ZG->GetSize(); ibin++)p_MC_ZG->SetBinContent(ibin, p_MC_ZG->GetBinContent(ibin)+radomMC*p_MC_ZG->GetBinError(ibin));
	p_MC_WG->Sumw2();
	p_MC_ZG->Sumw2();

	TCanvas *mccan=new TCanvas("mccan","",1200,600);
	mccan->Divide(2);
	mccan->cd(1);
	p_MC_WG->Draw();
	mccan->cd(2);
	p_MC_ZG->Draw();

  p_proxy->GetXaxis()->SetTitle("#Delta#phi");
  p_target->GetXaxis()->SetTitle("#Delta#phi");

//	p_rawsig->Rebin(2);
//	p_target->Rebin(2);
//	p_proxy->Rebin(2);
//	p_MC->Rebin(2);	

  RooRealVar dphi("dphi","",0,3.2);
  RooDataHist* h_target = new RooDataHist("h_target","h_target",dphi,p_target);
  RooDataHist* h_proxy  = new RooDataHist("h_proxy", "h_proxy", dphi,p_proxy );
  RooDataHist* h_MC_WG  = new RooDataHist("h_MC_WG", "h_MC_WG", dphi,p_MC_WG);
  RooDataHist* h_MC_ZG  = new RooDataHist("h_MC_ZG", "h_MC_ZG", dphi,p_MC_ZG);
  RooHistPdf*  pdf_proxy= new RooHistPdf("pdf_proxy","pdf_proxy",dphi, *h_proxy);
  RooHistPdf*  pdf_MC_WG= new RooHistPdf("pdf_MC_WG","pdf_MC_WG",dphi, *h_MC_WG);
  RooHistPdf*  pdf_MC_ZG= new RooHistPdf("pdf_MC_ZG","pdf_MC_ZG",dphi, *h_MC_ZG);
  RooRealVar fakefrac("fakefrac","fakefrac",0.3,0,1.0);
  RooRealVar mcfrac("mcfrac","mcfrac",0.7,0,1.0);
	RooAddPdf* pdf_MC_VG = new RooAddPdf("pdf_MC_VG", "", RooArgList(*pdf_MC_WG, *pdf_MC_ZG), mcfrac);
  RooAddPdf model("model","",RooArgList(*pdf_proxy, *pdf_MC_VG),RooArgList(fakefrac));

  RooPlot* frame = dphi.frame(RooFit::Title("#Delta#phi(l,MET) (40 GeV < MET < 70 GeV)"));
  TCanvas *can=new TCanvas("can","",600,600);
  can->cd(1);
	model.fitTo(*h_target,RooFit::SumW2Error(kTRUE),RooFit::Save());
  RooFitResult* result = model.fitTo(*h_target,RooFit::SumW2Error(kTRUE),RooFit::Save());
  h_target->plotOn(frame);
  model.plotOn(frame,
		   RooFit::FillColor(kBlue-4));
  model.plotOn(frame, RooFit::Components(*pdf_proxy),
               RooFit::LineStyle(kDashed),
               RooFit::LineColor(kRed),
               RooFit::Normalization(1.0));
  model.plotOn(frame, RooFit::Components(*pdf_MC_WG),
               RooFit::LineStyle(kDashed),
               RooFit::LineColor(kGreen),
               RooFit::Normalization(1.0));
  model.plotOn(frame, RooFit::Components(*pdf_MC_ZG),
               RooFit::LineStyle(kDashed),
               RooFit::LineColor(kMagenta),
               RooFit::Normalization(1.0));
  frame->Draw();
////    	std::ostringstream figurename;
////      figurename.str("");
////      figurename << "fit_lepPt_mg_met" << metlow << "_" << methigh << "_pt" <<  leplow << "_" << lephigh <<  "_iso" << isocut << ".pdf";
////      can->SaveAs(figurename.str().c_str());
////    
////    	TH1D *p_combine = new TH1D("p_combine","",32,0,3.2);
////    	TH1D *p_qcd  = new TH1D("p_qcd","",32,0,3.2);
////    	TH1D *p_VGAMMA = new TH1D("p_VGAMMA","",32,0,3.2);
////    	TH1D *norm_MC = (TH1D*)p_MC->Clone("norm_MC");
////    	TH1D *norm_proxy = (TH1D*)p_proxy->Clone("norm_proxy");
////    	norm_MC->Scale(1.0/norm_MC->Integral(1, norm_MC->GetSize()) );
////    	norm_proxy->Scale(1.0/norm_proxy->Integral(1, norm_proxy->GetSize()));
////    	p_combine->Add(norm_proxy, fakefrac.getVal());
////    	p_combine->Add(norm_MC, 1-fakefrac.getVal());
////    	p_combine->Scale(p_target->Integral(1, p_target->GetSize()));
////    	for(int ibin(1); ibin <= 32; ibin++)p_combine->SetBinError(ibin, fakefrac.getError()*p_combine->GetBinContent(ibin));
////    	p_qcd->Add(norm_proxy, fakefrac.getVal()*p_target->Integral(1, p_target->GetSize()));
////    	p_VGAMMA->Add(norm_MC, (1-fakefrac.getVal())*p_target->Integral(1, p_target->GetSize()));
////    	TH1D *fitratio = new TH1D("fitratio",";#Delta#phi(l, E_{T}^{miss});fit/data", 32,0,3.2);
////    	TGraphErrors *fitratio_error = new TGraphErrors(32);
////    	TGraphErrors *p_combine_error = new TGraphErrors(32);
////      double staterror = fakefrac.getError()/fakefrac.getVal();
////    	for(int ibin(1); ibin <= 32; ibin++)fitratio->SetBinContent(ibin, p_combine->GetBinContent(ibin)/p_target->GetBinContent(ibin));
////    	for(int ibin(1); ibin <= 32; ibin++)fitratio->SetBinError(ibin, p_target->GetBinError(ibin)/p_target->GetBinContent(ibin));
////    	for(int ibin(0); ibin < 32; ibin++)fitratio_error->SetPoint(ibin, 0.05+ibin*0.1, 1);
////    	for(int ibin(0); ibin < 32; ibin++)fitratio_error->SetPointError(ibin, 0.05, staterror);
////    	for(int ibin(0); ibin < 32; ibin++)p_combine_error->SetPoint(ibin, -0.05+ibin*0.1, p_combine->GetBinContent(ibin));
////    	for(int ibin(0); ibin < 32; ibin++)p_combine_error->SetPointError(ibin, 0.05, staterror*p_combine->GetBinContent(ibin));
////    
////    
////    //	RooDataHist* h_combine = new RooDataHist("h_combine","h_combine",dphi, p_combine);
////    //	h_combine->plotOn(frame,
////    //       RooFit::MarkerColor(kGreen-4));
////      TCanvas *canres=new TCanvas("canres","",600,600);
////      canres->cd(1);
////    	TPad *canpt_pad1 = new TPad("canpt_pad1", "pad1", 0, 0.3, 1, 1.0);
////    	canpt_pad1->SetBottomMargin(0);
////    	canpt_pad1->Draw();          
////    	canpt_pad1->cd();   
////    	p_target->SetTitle("#Delta#phi(#mu, E_{T}^{miss})"); 
////    	//p_target->SetMaximum(500);
////    	p_target->SetMinimum(0);
////    	p_target->SetMarkerStyle(20);
////    	p_target->Draw("EP");
////    	p_combine_error->SetFillColor(kBlue);
////    	p_combine_error->SetFillStyle(3001);
////    	p_combine_error->Draw("E2");
////    	p_target->Draw("EP same");
////    	p_combine->SetLineColor(kBlue);
////    	p_combine->SetLineWidth(2);
////    	p_combine->Draw("hist same");
////    	p_qcd->SetLineColor(kRed);
////    	p_qcd->Draw("hist same");
////    	p_VGAMMA->SetLineColor(kGreen);
////    	p_VGAMMA->Draw("hist same");
////    	
////    	TLegend *leg=new TLegend(0.6,0.1,0.85,0.3);
////    	leg->AddEntry(p_target, "observed");
////    	leg->AddEntry(p_combine, "fit result");
////    	leg->AddEntry(p_qcd, "fake leptons");
////    	leg->Draw("same");
////    	 
////    	canres->cd();   
////    	TPad *canpt_pad2 = new TPad("canpt_pad2", "pad2", 0, 0.05, 1, 0.3);
////    	canpt_pad2->SetTopMargin(0);
////    	canpt_pad2->SetBottomMargin(0.3);
////    	canpt_pad2->Draw();
////    	canpt_pad2->cd(); 	
////    	TH1D *dummy_ptratio = new TH1D("dummy_ptratio",";#Delta#phi(l, E_{T}^{miss});fit/data",32,0,3.2);
////    	dummy_ptratio->SetMaximum(2);
////    	dummy_ptratio->SetMinimum(0);
////    	dummy_ptratio->GetXaxis()->SetLabelFont(63);
////    	dummy_ptratio->GetXaxis()->SetLabelSize(14);
////    	dummy_ptratio->GetYaxis()->SetLabelFont(63);
////    	dummy_ptratio->GetYaxis()->SetLabelSize(14);
////    	dummy_ptratio->Draw();
////      TLine *flatratio = new TLine(0,1,3.2,1);
////    	flatratio->Draw("same");
////    	fitratio->SetMarkerStyle(20);
////    	fitratio->Draw("EP same");
////    	fitratio_error->SetFillStyle(3005);
////    	fitratio_error->Draw("E2 same");	
////    
////    	std::ostringstream plotname;
////    	plotname << "fit_dPhi_mg_" << leplow << "_" << lephigh << "_met40" << "_iso" << isocut << ".pdf";
////      if(ih == 0)canres->SaveAs(plotname.str().c_str());
////    
////      double nMCtotal(0);
////      for(unsigned i(1); i<=32; i++)nMCtotal+= p_MC->GetBinContent(i);
////      std::cout << "factor for fake: " << fakefrac.getVal()*p_target->Integral(1, p_target->GetSize()) << "/" << p_proxy->Integral(1, p_proxy->GetSize()) << std::endl;
////      std::cout << "factor for MC: " << (1-fakefrac.getVal())*p_target->Integral(1, p_target->GetSize()) << "/" << p_MC->Integral(1, p_MC->GetSize())  << std::endl;
////      std::cout << "fakefrac=" << fakefrac.getVal()<< std::endl;
////    
////    	
////    	float fakescale = fakefrac.getVal()*p_target->Integral(1, p_target->GetSize())/p_proxy->Integral(1, p_proxy->GetSize());
////    	float fakescaleerror = fakefrac.getError()*p_target->Integral(1, p_target->GetSize())/p_proxy->Integral(1, p_proxy->GetSize());
////    	float vgammascale = (1-fakefrac.getVal())*p_target->Integral(1, p_target->GetSize())/p_MC->Integral(1, p_MC->GetSize());
////    	float vgammascaleerror = fakefrac.getError()*p_target->Integral(1, p_target->GetSize())/p_MC->Integral(1, p_MC->GetSize()); 
////    	myfile << leplow << " " << lephigh << " " << fakescale << " " << fakescaleerror  << " " << vgammascale << " " << vgammascaleerror << std::endl; 
////    	myfile.close();
////    
////    
////            TH1D *p_LepPt_qcd = (TH1D*)file_qcd->Get("p_LepPt");
////            TH1D *p_LepPt_sig = (TH1D*)file_sig->Get("p_LepPt");
////            TH1D *p_LepPt_ele = (TH1D*)file_ele->Get("p_LepPt");
////            TH1D *p_LepPt_jet = (TH1D*)file_jet->Get("p_LepPt");
////            TH1D *p_LepPt_rare = (TH1D*)file_rare->Get("p_LepPt");
////            TH1D *p_LepPt_MC = (TH1D*)file_VG->Get("p_LepPt");
////            p_LepPt_sig->Add(p_LepPt_ele, -1);
////            p_LepPt_sig->Add(p_LepPt_jet, -1);
////            p_LepPt_sig->Add(p_LepPt_rare, -1);
////    				p_LepPt_sig->Add(p_LepPt_MC, -1*vgammascale);
////    
////    				p_LepPt_sig->Rebin(2);
////    				p_LepPt_qcd->Rebin(2);
////    
////            TCanvas *canscale = new TCanvas("canscale","", 600,600);
////            canscale->cd();
////    				TPad *canscale_pad1 = new TPad("canscale_pad1", "pad1", 0, 0.3, 1, 1.0);
////    				canscale_pad1->SetBottomMargin(0);
////    				canscale_pad1->Draw();          
////    				canscale_pad1->cd();   
////            gPad->SetLogy();
////            p_LepPt_sig->SetMarkerStyle(20);
////    				p_LepPt_qcd->SetMinimum(1);
////    			  p_LepPt_qcd->SetLineColor(kRed);
////            p_LepPt_qcd->Draw();
////            p_LepPt_sig->Draw("EP same");
////    
////    				TH1D *p_scale = (TH1D*)p_LepPt_sig->Clone("p_scale");
////    				p_scale->Divide(p_LepPt_qcd);
////    				canscale->cd();
////    				TPad *canscale_pad2 = new TPad("canscale_pad2", "pad2", 0, 0, 1, 0.3);
////            canscale_pad2->SetBottomMargin(0.05);
////            canscale_pad2->Draw();          
////            canscale_pad2->cd();
////    				p_scale->Draw(); 
////    //        plotname.str("");
////    //        plotname << "fit_lepPt_mg_met" << metlow << "_" << methigh << "_pt" <<  leplow << "_" << lephigh <<  "_iso" << isocut << ".pdf";
////    //        canscale->SaveAs(plotname.str().c_str());
////    //

  return 1;
}
