#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TF3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TChain.h"
#include "TSystem.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TGraphErrors.h"
void plot_bkg(){//main  

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");
	
//	TFile *file_sig = TFile::Open("bkgTree_eg_signal_ReMiniAOD.root");
//	TFile *file_ele = TFile::Open("bkgTree_eg_eleBkg_ReMiniAOD.root");
//	TFile *file_jet = TFile::Open("bkgTree_eg_jetbkg_ReMiniAOD.root");
//	TFile *file_qcd = TFile::Open("bkgTree_eg_qcd_ReMiniAOD.root");
//	TFile *file_VG  = TFile::Open("bkgTree_eg_VGBkg_ReMiniAOD.root");
//	TFile *file_rare= TFile::Open("bkgTree_eg_rareBkg_ReMiniAOD.root");
//

	TFile *file_sig = TFile::Open("controlTree_eg_signal_ReMiniAOD0_70.root");
	TFile *file_ele = TFile::Open("controlTree_eg_eleBkg_ReMiniAOD0_70.root");
	TFile *file_jet = TFile::Open("controlTree_eg_jetbkg_ReMiniAOD0_70.root");
	TFile *file_qcd = TFile::Open("controlTree_eg_qcd_ReMiniAOD_0_70_met40_iso4.root");
	TFile *file_VG  = TFile::Open("controlTree_eg_VGBkg_ReMiniAOD0_70.root");
	TFile *file_rare= TFile::Open("controlTree_eg_rareBkg_ReMiniAOD0_70.root");

	Double_t plotEtBins[]={35,40,45,50,55,60,65,70,75,80, 85,90,95,100,105,110,115,120,125,130, 135,140,146,152,158,164,170,177,184,192, 200,208,216,224,232,240,250,260,275,290, 305,325,345,370,400};
	Double_t plotPtBins[]={25,30,35,40,45,50,55,60,65,70,75,80, 85,90,95,100,105,110,115,120,125,130, 135,140,146,152,158,164,170,177,184,192, 200,208,216,224,232,240,250,260,275,290, 305,325,345,370,400};
	TGraphErrors *error_PhoEt = new TGraphErrors(44);
	TGraphErrors *error_LepPt = new TGraphErrors(46);
	TGraphErrors *error_MET = new TGraphErrors(100); 
	TGraphErrors *error_Mt = new TGraphErrors(200); 
	TGraphErrors *error_dPhiEleMET = new TGraphErrors(64);

	TH1F *p_allPhoEt = (TH1F*)file_sig->Get("p_PhoEt");
	TH1F *p_allLepPt = (TH1F*)file_sig->Get("p_LepPt");
	TH1F *p_allMET   = (TH1F*)file_sig->Get("p_MET");
	TH1F *p_allMt    = (TH1F*)file_sig->Get("p_Mt");
	TH1F *p_alldPhi  = (TH1F*)file_sig->Get("p_dPhiEleMET");
	TH1F *p_allPU    = (TH1F*)file_sig->Get("p_PU");
	
	
	TH1F *p_elePhoEt = (TH1F*)file_ele->Get("p_PhoEt");
	TH1F *p_eleLepPt = (TH1F*)file_ele->Get("p_LepPt");
	TH1F *p_eleMET   = (TH1F*)file_ele->Get("p_MET");
	TH1F *p_eleMt    = (TH1F*)file_ele->Get("p_Mt");
	TH1F *p_eledPhi  = (TH1F*)file_ele->Get("p_dPhiEleMET");
	TH1F *p_elePU    = (TH1F*)file_ele->Get("p_PU");

	TH1F *p_jetPhoEt = (TH1F*)file_jet->Get("p_PhoEt");
	TH1F *p_jetLepPt = (TH1F*)file_jet->Get("p_LepPt");
	TH1F *p_jetMET   = (TH1F*)file_jet->Get("p_MET");
	TH1F *p_jetMt    = (TH1F*)file_jet->Get("p_Mt");
	TH1F *p_jetdPhi  = (TH1F*)file_jet->Get("p_dPhiEleMET");
	TH1F *p_jetPU    = (TH1F*)file_jet->Get("p_PU");

	TH1F *p_qcdPhoEt = (TH1F*)file_qcd->Get("p_PhoEt");
	TH1F *p_qcdLepPt = (TH1F*)file_qcd->Get("p_LepPt");
	TH1F *p_qcdMET   = (TH1F*)file_qcd->Get("p_MET");
	TH1F *p_qcdMt    = (TH1F*)file_qcd->Get("p_Mt");
	TH1F *p_qcddPhi  = (TH1F*)file_qcd->Get("p_dPhiEleMET");
	TH1F *p_qcdPU    = (TH1F*)file_qcd->Get("p_PU");

	TH1F *p_VGPhoEt = (TH1F*)file_VG->Get("p_PhoEt");
	TH1F *p_VGLepPt = (TH1F*)file_VG->Get("p_LepPt");
	TH1F *p_VGMET   = (TH1F*)file_VG->Get("p_MET");
	TH1F *p_VGMt    = (TH1F*)file_VG->Get("p_Mt");
	TH1F *p_VGdPhi  = (TH1F*)file_VG->Get("p_dPhiEleMET");
	TH1F *p_VGPU    = (TH1F*)file_VG->Get("p_PU");

	TH1F *p_rarePhoEt = (TH1F*)file_rare->Get("p_PhoEt");
	TH1F *p_rareLepPt = (TH1F*)file_rare->Get("p_LepPt");
	TH1F *p_rareMET   = (TH1F*)file_rare->Get("p_MET");
	TH1F *p_rareMt    = (TH1F*)file_rare->Get("p_Mt");
	TH1F *p_raredPhi  = (TH1F*)file_rare->Get("p_dPhiEleMET");
	TH1F *p_rarePU    = (TH1F*)file_rare->Get("p_PU");


//	p_allMET->Rebin(5);
//	p_eleMET->Rebin(5);
//	p_jetMET->Rebin(5);
//	p_qcdMET->Rebin(5);
//	p_VGMET->Rebin(5);
//	p_rareMET->Rebin(5);

	p_allMt->Rebin(5);
	p_eleMt->Rebin(5);
	p_jetMt->Rebin(5);
	p_qcdMt->Rebin(5);
	p_VGMt->Rebin(5);
	p_rareMt->Rebin(5);

	gStyle->SetOptStat(0);
	TCanvas *c_pt = new TCanvas("Photon_Pt", "Photon P_{T}",600,600);
	c_pt->cd();
	TPad *pt_pad1 = new TPad("pt_pad1", "pt_pad1", 0, 0.3, 1, 1.0);
	pt_pad1->SetBottomMargin(0.1);
	pt_pad1->Draw();  
	pt_pad1->cd();  
	gPad->SetLogy();
	p_allPhoEt->SetMinimum(1);
	p_allPhoEt->GetXaxis()->SetRangeUser(35,300);
	p_allPhoEt->SetLineColor(kBlack);
	p_allPhoEt->SetMarkerStyle(20);
	p_allPhoEt->Draw("P");
	p_VGPhoEt->SetFillStyle(1001);
	p_VGPhoEt->SetLineColor(kMagenta);
	p_VGPhoEt->SetFillColor(kMagenta);
	p_rarePhoEt->SetFillStyle(1001);
	p_rarePhoEt->SetLineColor(kYellow-4);
	p_rarePhoEt->SetFillColor(kYellow-4);
	p_qcdPhoEt->SetFillStyle(1001);
	p_qcdPhoEt->SetLineColor(kBlue);
	p_qcdPhoEt->SetFillColor(kBlue);
	p_elePhoEt->SetFillStyle(1001);
	p_elePhoEt->SetLineColor(kRed);
	p_elePhoEt->SetFillColor(kRed);
	p_jetPhoEt->SetFillStyle(1001);
	p_jetPhoEt->SetLineColor(kGreen);
	p_jetPhoEt->SetFillColor(kGreen);
	for(int ibin(1); ibin < p_VGPhoEt->GetSize(); ibin++){
		error_PhoEt->SetPoint(ibin-1,p_VGPhoEt->GetBinCenter(ibin), p_VGPhoEt->GetBinContent(ibin));
		error_PhoEt->SetPointError(ibin-1,(p_VGPhoEt->GetBinLowEdge(ibin+1)-p_VGPhoEt->GetBinLowEdge(ibin))/2,p_VGPhoEt->GetBinError(ibin));
		std::cout << p_elePhoEt->GetBinError(ibin) << " " << p_jetPhoEt->GetBinError(ibin) << " " << p_qcdPhoEt->GetBinError(ibin) << " " << p_rarePhoEt->GetBinError(ibin)  << std::endl;
	}
	p_elePhoEt->Add(p_rarePhoEt); // ele 2nd
	p_jetPhoEt->Add(p_elePhoEt);  // jet 3rd
	p_qcdPhoEt->Add(p_jetPhoEt);  // qcd 4th
	p_VGPhoEt->Add(p_qcdPhoEt);   // VG  5th
	p_VGPhoEt->Sumw2();
	p_VGPhoEt->Draw("hist same");
	p_qcdPhoEt->Draw("hist same");
	p_jetPhoEt->Draw("hist same");
	p_elePhoEt->Draw("hist same");
	p_rarePhoEt->Draw("hist same");
  error_PhoEt->SetFillColor(2);
  error_PhoEt->SetFillStyle(3001);
	error_PhoEt->Draw("E2 same");
	TLegend *leg_pt =  new TLegend(0.6,0.75,0.9,0.9);
	leg_pt->SetFillStyle(0);
	gStyle->SetLegendBorderSize(1);
	gStyle->SetLegendFillColor(0);
	leg_pt->AddEntry(p_allPhoEt,"observed (MET < 100 GeV)");
	leg_pt->AddEntry(p_rarePhoEt,"t#bar{t}#gamma/WW#gamma/WZ#gamma");
	leg_pt->AddEntry(p_elePhoEt,"e->#gamma fake");
	leg_pt->AddEntry(p_jetPhoEt,"j->#gamma fake");
	leg_pt->AddEntry(p_qcdPhoEt,"j->e fake");
	leg_pt->AddEntry(p_VGPhoEt, "WG/ZG");
	leg_pt->Draw("same");
	p_allPhoEt->Draw("E same");

	c_pt->cd();
	TPad *pt_pad2 = new TPad("pt_pad2", "pt_pad2", 0, 0.05, 1, 0.25);
	pt_pad2->Draw();
	pt_pad2->cd();
  TLine *flatratio = new TLine(35,1,400,1);
	TH1F *ratio=(TH1F*)p_VGPhoEt->Clone("transfer factor");
	ratio->SetMarkerStyle(20);
	ratio->SetLineColor(kBlack);
	ratio->Divide(p_allPhoEt);
	ratio->SetTitle("");
	ratio->GetYaxis()->SetTitle("observed/bkg");
	ratio->GetXaxis()->SetLabelFont(63);
	ratio->GetXaxis()->SetLabelSize(14);
	ratio->GetYaxis()->SetLabelFont(63);
	ratio->GetYaxis()->SetLabelSize(14);
	ratio->Draw();
	flatratio->Draw("same");
	c_pt->SaveAs("BKG_eg_2016ReMiniAOD_pt.pdf");

// ******** MET ************************//
	gStyle->SetOptStat(0);
	TCanvas *c_met = new TCanvas("MET", "MET",600,600);
	c_met->cd();
	TPad *met_pad1 = new TPad("met_pad1", "met_pad1", 0, 0.3, 1, 1.0);
	met_pad1->SetBottomMargin(0.1);
	met_pad1->Draw();  
	met_pad1->cd();  
	gPad->SetLogy();
	p_allMET->SetMinimum(1);
	p_allMET->GetXaxis()->SetRangeUser(0,400);
	p_allMET->SetLineColor(kBlack);
	p_allMET->SetMarkerStyle(20);
	p_allMET->Draw("P");
	p_VGMET->SetFillStyle(1001);
	p_VGMET->SetLineColor(kMagenta);
	p_VGMET->SetFillColor(kMagenta);
	p_rareMET->SetFillStyle(1001);
	p_rareMET->SetLineColor(kYellow-4);
	p_rareMET->SetFillColor(kYellow-4);
	p_qcdMET->SetFillStyle(1001);
	p_qcdMET->SetLineColor(kBlue);
	p_qcdMET->SetFillColor(kBlue);
	p_eleMET->SetFillStyle(1001);
	p_eleMET->SetLineColor(kRed);
	p_eleMET->SetFillColor(kRed);
	p_jetMET->SetFillStyle(1001);
	p_jetMET->SetLineColor(kGreen);
	p_jetMET->SetFillColor(kGreen);
	p_eleMET->Add(p_rareMET); // ele 2nd
	p_jetMET->Add(p_eleMET);  // jet 3rd
	p_qcdMET->Add(p_jetMET);  // qcd 4th
	p_VGMET->Add(p_qcdMET);   // VG  5th
	for(int ibin(1); ibin < p_VGMET->GetSize(); ibin++){
		error_MET->SetPoint(ibin-1,p_VGMET->GetBinCenter(ibin), p_VGMET->GetBinContent(ibin));
		error_MET->SetPointError(ibin-1,(p_VGMET->GetBinLowEdge(ibin+1)-p_VGMET->GetBinLowEdge(ibin))/2,p_VGMET->GetBinError(ibin));
	}
	p_VGMET->Draw("hist same");
	p_qcdMET->Draw("hist same");
	p_jetMET->Draw("hist same");
	p_eleMET->Draw("hist same");
	p_rareMET->Draw("hist same");
  error_MET->SetFillColor(2);
  error_MET->SetFillStyle(3001);
	error_MET->Draw("E2 same");
	TLegend *leg_met =  new TLegend(0.6,0.75,0.9,0.9);
	leg_met->SetFillStyle(0);
	gStyle->SetLegendBorderSize(1);
	gStyle->SetLegendFillColor(0);
	leg_met->AddEntry(p_allMET,"observed (MET < 100 GeV)");
	leg_met->AddEntry(p_rareMET,"t#bar{t}#gamma/WW#gamma/WZ#gamma");
	leg_met->AddEntry(p_eleMET,"e->#gamma fake");
	leg_met->AddEntry(p_jetMET,"j->#gamma fake");
	leg_met->AddEntry(p_qcdMET,"j->e fake");
	leg_met->AddEntry(p_VGMET, "WG/ZG");
	leg_met->Draw("same");
	p_allMET->Draw("E same");

	c_met->cd();
	TPad *met_pad2 = new TPad("met_pad2", "met_pad2", 0, 0.05, 1, 0.25);
	met_pad2->Draw();
	met_pad2->cd();
  TLine *flatratio_met = new TLine(0,1,100,1);
	TH1F *ratio_met=(TH1F*)p_VGMET->Clone("transfer factor");
	ratio_met->GetXaxis()->SetRangeUser(0,400);
	ratio_met->SetLineColor(kBlack);
	ratio_met->SetMarkerStyle(20);
	ratio_met->Divide(p_allMET);
	ratio_met->SetTitle("");
	ratio_met->GetYaxis()->SetTitle("observed/bkg");
	ratio_met->GetXaxis()->SetLabelFont(63);
	ratio_met->GetXaxis()->SetLabelSize(14);
	ratio_met->GetYaxis()->SetLabelFont(63);
	ratio_met->GetYaxis()->SetLabelSize(14);
	ratio_met->Draw();
	flatratio_met->Draw("same");
	c_met->SaveAs("BKG_eg_2016ReMiniAOD_met.pdf");


// ******** LepPt ************************//
	gStyle->SetOptStat(0);
	TCanvas *c_leppt = new TCanvas("LepPt", "LepPt",600,600);
	c_leppt->cd();
	TPad *leppt_pad1 = new TPad("leppt_pad1", "leppt_pad1", 0, 0.3, 1, 1.0);
	leppt_pad1->SetBottomMargin(0.1);
	leppt_pad1->Draw();  
	leppt_pad1->cd();  
	gPad->SetLogy();
	p_allLepPt->SetMinimum(1);
	p_allLepPt->GetXaxis()->SetRangeUser(35,300);
	p_allLepPt->SetLineColor(kBlack);
	p_allLepPt->SetMarkerStyle(20);
	p_allLepPt->Draw("P");
	p_VGLepPt->SetFillStyle(1001);
	p_VGLepPt->SetLineColor(kMagenta);
	p_VGLepPt->SetFillColor(kMagenta);
	p_rareLepPt->SetFillStyle(1001);
	p_rareLepPt->SetLineColor(kYellow-4);
	p_rareLepPt->SetFillColor(kYellow-4);
	p_qcdLepPt->SetFillStyle(1001);
	p_qcdLepPt->SetLineColor(kBlue);
	p_qcdLepPt->SetFillColor(kBlue);
	p_eleLepPt->SetFillStyle(1001);
	p_eleLepPt->SetLineColor(kRed);
	p_eleLepPt->SetFillColor(kRed);
	p_jetLepPt->SetFillStyle(1001);
	p_jetLepPt->SetLineColor(kGreen);
	p_jetLepPt->SetFillColor(kGreen);
	p_eleLepPt->Add(p_rareLepPt); // ele 2nd
	p_jetLepPt->Add(p_eleLepPt);  // jet 3rd
	p_qcdLepPt->Add(p_jetLepPt);  // qcd 4th
	p_VGLepPt->Add(p_qcdLepPt);   // VG  5th
	for(int ibin(1); ibin < p_VGLepPt->GetSize(); ibin++){
		error_LepPt->SetPoint(ibin-1,p_VGLepPt->GetBinCenter(ibin), p_VGLepPt->GetBinContent(ibin));
		error_LepPt->SetPointError(ibin-1,(p_VGLepPt->GetBinLowEdge(ibin+1)-p_VGLepPt->GetBinLowEdge(ibin))/2,p_VGLepPt->GetBinError(ibin));
	}
	p_VGLepPt->Draw("hist same");
	p_qcdLepPt->Draw("hist same");
	p_jetLepPt->Draw("hist same");
	p_eleLepPt->Draw("hist same");
	p_rareLepPt->Draw("hist same");
  error_LepPt->SetFillColor(2);
  error_LepPt->SetFillStyle(3001);
	error_LepPt->Draw("E2 same");
	TLegend *leg_leppt =  new TLegend(0.6,0.75,0.9,0.9);
	leg_leppt->SetFillStyle(0);
	gStyle->SetLegendBorderSize(1);
	gStyle->SetLegendFillColor(0);
	leg_leppt->AddEntry(p_allLepPt,"observed (LepPt < 100 GeV)");
	leg_leppt->AddEntry(p_rareLepPt,"t#bar{t}#gamma/WW#gamma/WZ#gamma");
	leg_leppt->AddEntry(p_eleLepPt,"e->#gamma fake");
	leg_leppt->AddEntry(p_jetLepPt,"j->#gamma fake");
	leg_leppt->AddEntry(p_qcdLepPt,"j->e fake");
	leg_leppt->AddEntry(p_VGLepPt, "WG/ZG");
	leg_leppt->Draw("same");
	p_allLepPt->Draw("E same");

	c_leppt->cd();
	TPad *leppt_pad2 = new TPad("leppt_pad2", "leppt_pad2", 0, 0.05, 1, 0.25);
	leppt_pad2->Draw();
	leppt_pad2->cd();
  TLine *flatratio_leppt = new TLine(25,1,400,1);
	TH1F *ratio_leppt=(TH1F*)p_VGLepPt->Clone("transfer factor");
	ratio_leppt->SetMarkerStyle(20);
	ratio_leppt->SetLineColor(kBlack);
	ratio_leppt->Divide(p_allLepPt);
	ratio_leppt->SetTitle("");
	ratio_leppt->GetYaxis()->SetTitle("observed/bkg");
	ratio_leppt->GetXaxis()->SetLabelFont(63);
	ratio_leppt->GetXaxis()->SetLabelSize(14);
	ratio_leppt->GetYaxis()->SetLabelFont(63);
	ratio_leppt->GetYaxis()->SetLabelSize(14);
	ratio_leppt->Draw();
	flatratio_leppt->Draw("same");
	c_leppt->SaveAs("BKG_eg_2016ReMiniAOD_leppt.pdf");
// ******** Mt ************************//
	gStyle->SetOptStat(0);
	TCanvas *c_mt = new TCanvas("Mt", "Mt",600,600);
	c_mt->cd();
	TPad *mt_pad1 = new TPad("mt_pad1", "mt_pad1", 0, 0.3, 1, 1.0);
	mt_pad1->SetBottomMargin(0.1);
	mt_pad1->Draw();  
	mt_pad1->cd();  
	gPad->SetLogy();
	p_allMt->SetMinimum(1);
  p_allMt->SetMaximum(1000000);
	p_allMt->GetXaxis()->SetRangeUser(0,400);
	p_allMt->SetLineColor(kBlack);
	p_allMt->SetMarkerStyle(20);
	p_allMt->Draw("P");
	p_VGMt->SetFillStyle(1001);
	p_VGMt->SetLineColor(kMagenta);
	p_VGMt->SetFillColor(kMagenta);
	p_rareMt->SetFillStyle(1001);
	p_rareMt->SetLineColor(kYellow-4);
	p_rareMt->SetFillColor(kYellow-4);
	p_qcdMt->SetFillStyle(1001);
	p_qcdMt->SetLineColor(kBlue);
	p_qcdMt->SetFillColor(kBlue);
	p_eleMt->SetFillStyle(1001);
	p_eleMt->SetLineColor(kRed);
	p_eleMt->SetFillColor(kRed);
	p_jetMt->SetFillStyle(1001);
	p_jetMt->SetLineColor(kGreen);
	p_jetMt->SetFillColor(kGreen);
	p_eleMt->Add(p_rareMt); // ele 2nd
	p_jetMt->Add(p_eleMt);  // jet 3rd
	p_qcdMt->Add(p_jetMt);  // qcd 4th
	p_VGMt->Add(p_qcdMt);   // VG  5th
	for(int ibin(1); ibin < p_VGMt->GetSize(); ibin++){
		error_Mt->SetPoint(ibin-1,p_VGMt->GetBinCenter(ibin), p_VGMt->GetBinContent(ibin));
		error_Mt->SetPointError(ibin-1,(p_VGMt->GetBinLowEdge(ibin+1)-p_VGMt->GetBinLowEdge(ibin))/2,p_VGMt->GetBinError(ibin));
	}
	p_VGMt->Draw("hist same");
	p_qcdMt->Draw("hist same");
	p_jetMt->Draw("hist same");
	p_eleMt->Draw("hist same");
	p_rareMt->Draw("hist same");
  error_Mt->SetFillColor(2);
  error_Mt->SetFillStyle(3001);
	error_Mt->Draw("E2 same");
	TLegend *leg_mt =  new TLegend(0.6,0.75,0.9,0.9);
	leg_mt->SetFillStyle(0);
	gStyle->SetLegendBorderSize(1);
	gStyle->SetLegendFillColor(0);
	leg_mt->AddEntry(p_allMt,"observed (Mt < 100 GeV)");
	leg_mt->AddEntry(p_rareMt,"t#bar{t}#gamma/WW#gamma/WZ#gamma");
	leg_mt->AddEntry(p_eleMt,"e->#gamma fake");
	leg_mt->AddEntry(p_jetMt,"j->#gamma fake");
	leg_mt->AddEntry(p_qcdMt,"j->e fake");
	leg_mt->AddEntry(p_VGMt, "WG/ZG");
	leg_mt->Draw("same");
	p_allMt->Draw("E same");

	c_mt->cd();
	TPad *mt_pad2 = new TPad("mt_pad2", "mt_pad2", 0, 0.05, 1, 0.25);
	mt_pad2->Draw();
	mt_pad2->cd();
  TLine *flatratio_mt = new TLine(0,1,400,1);
	TH1F *ratio_mt=(TH1F*)p_VGMt->Clone("transfer factor");
	ratio_mt->SetMarkerStyle(20);
	ratio_mt->SetLineColor(kBlack);
	ratio_mt->GetXaxis()->SetRangeUser(0,400);
	ratio_mt->SetMinimum(0);
	ratio_mt->SetMaximum(2);
	ratio_mt->Divide(p_allMt);
	ratio_mt->SetTitle("");
	ratio_mt->GetYaxis()->SetTitle("observed/bkg");
	ratio_mt->GetXaxis()->SetLabelFont(63);
	ratio_mt->GetXaxis()->SetLabelSize(14);
	ratio_mt->GetYaxis()->SetLabelFont(63);
	ratio_mt->GetYaxis()->SetLabelSize(14);
	ratio_mt->Draw();
	flatratio_mt->Draw("same");
	c_mt->SaveAs("BKG_eg_2016ReMiniAOD_mt.pdf");



// ******** dPhi ************************//
	gStyle->SetOptStat(0);
	TCanvas *c_dphi = new TCanvas("dPhi", "dPhi",600,600);
	c_dphi->cd();
	TPad *dphi_pad1 = new TPad("dphi_pad1", "dphi_pad1", 0, 0.3, 1, 1.0);
	dphi_pad1->SetBottomMargin(0.1);
	dphi_pad1->Draw();  
	dphi_pad1->cd();  
	p_alldPhi->SetMinimum(1);
	p_alldPhi->SetLineColor(kBlack);
	p_alldPhi->SetMarkerStyle(20);
	p_alldPhi->Draw("P");
	p_VGdPhi->SetFillStyle(1001);
	p_VGdPhi->SetLineColor(kMagenta);
	p_VGdPhi->SetFillColor(kMagenta);
	p_raredPhi->SetFillStyle(1001);
	p_raredPhi->SetLineColor(kYellow-4);
	p_raredPhi->SetFillColor(kYellow-4);
	p_qcddPhi->SetFillStyle(1001);
	p_qcddPhi->SetLineColor(kBlue);
	p_qcddPhi->SetFillColor(kBlue);
	p_eledPhi->SetFillStyle(1001);
	p_eledPhi->SetLineColor(kRed);
	p_eledPhi->SetFillColor(kRed);
	p_jetdPhi->SetFillStyle(1001);
	p_jetdPhi->SetLineColor(kGreen);
	p_jetdPhi->SetFillColor(kGreen);
	p_eledPhi->Add(p_raredPhi); // ele 2nd
	p_jetdPhi->Add(p_eledPhi);  // jet 3rd
//	p_qcddPhi->Add(p_jetdPhi);  // qcd 4th
//	p_VGdPhi->Add(p_qcddPhi);   // VG  5th
//	p_VGdPhi->Draw("hist same");
//	p_qcddPhi->Draw("hist same");
	p_jetdPhi->Draw("hist same");
	p_eledPhi->Draw("hist same");
	p_raredPhi->Draw("hist same");
	TLegend *leg_dphi =  new TLegend(0.6,0.75,0.9,0.9);
	leg_dphi->SetFillStyle(0);
	gStyle->SetLegendBorderSize(1);
	gStyle->SetLegendFillColor(0);
	leg_dphi->AddEntry(p_alldPhi,"observed (dPhi < 100 GeV)");
	leg_dphi->AddEntry(p_raredPhi,"t#bar{t}#gamma/WW#gamma/WZ#gamma");
	leg_dphi->AddEntry(p_eledPhi,"e->#gamma fake");
	leg_dphi->AddEntry(p_jetdPhi,"j->#gamma fake");
	leg_dphi->AddEntry(p_qcddPhi,"j->e fake");
	leg_dphi->AddEntry(p_VGdPhi, "WG/ZG");
	leg_dphi->Draw("same");
	p_alldPhi->Draw("E same");

	c_dphi->cd();
	TPad *dphi_pad2 = new TPad("dphi_pad2", "dphi_pad2", 0, 0.05, 1, 0.25);
	dphi_pad2->Draw();
	dphi_pad2->cd();
  TLine *flatratio_dphi = new TLine(0,1,400,1);
	TH1F *ratio_dphi=(TH1F*)p_VGdPhi->Clone("transfer factor");
	ratio_dphi->SetMarkerStyle(20);
	ratio_dphi->SetLineColor(kBlack);
	ratio_dphi->GetXaxis()->SetRangeUser(0,400);
	ratio_dphi->SetMinimum(0);
	ratio_dphi->SetMaximum(2);
	ratio_dphi->Divide(p_alldPhi);
	ratio_dphi->SetTitle("");
	ratio_dphi->GetYaxis()->SetTitle("observed/bkg");
	ratio_dphi->GetXaxis()->SetLabelFont(63);
	ratio_dphi->GetXaxis()->SetLabelSize(14);
	ratio_dphi->GetYaxis()->SetLabelFont(63);
	ratio_dphi->GetYaxis()->SetLabelSize(14);
	ratio_dphi->Draw();
	flatratio_dphi->Draw("same");
	c_dphi->SaveAs("BKG_eg_2016ReMiniAOD_dphi.pdf");
}


