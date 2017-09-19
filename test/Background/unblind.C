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
void unblind(){//main  

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

	TFile *file_sig = TFile::Open("signalTree_mg_signal.root");
	TFile *file_ele = TFile::Open("signalTree_mg_eleBkg.root");
	TFile *file_jet = TFile::Open("signalTree_mg_jetbkg.root");
	TFile *file_qcd = TFile::Open("signalTree_mg_qcd.root");
	TFile *file_VG  = TFile::Open("signalTree_mg_VGBkg.root");
	TFile *file_rare= TFile::Open("signalTree_mg_rareBkg.root");

	TFile *file_t5 = TFile::Open("signalTree_T5WG.root");
	TFile *file_tchi=TFile::Open("signalTree_TChiWG.root");
	TH1D *p_t5wg_MET_valid_1000_500 = (TH1D*)file_t5->Get("p_t5wg_MET_valid_1000_500");
	TH1D *p_t5wg_MET_valid_1800_1000= (TH1D*)file_t5->Get("p_t5wg_MET_valid_1800_1000");
	TH1D *p_t5wg_MET_valid_1800_200 = (TH1D*)file_t5->Get("p_t5wg_MET_valid_1800_200");
	TH1D *p_t5wg_MET_signal_1000_500= (TH1D*)file_t5->Get("p_t5wg_MET_signal_1000_500");
	TH1D *p_t5wg_MET_signal_1800_1000=(TH1D*)file_t5->Get("p_t5wg_MET_signal_1800_1000");
	TH1D *p_t5wg_MET_signal_1800_200= (TH1D*)file_t5->Get("p_t5wg_MET_signal_1800_200");
	TH1D *p_t5wg_HT_signal_1000_500 = (TH1D*)file_t5->Get("p_t5wg_HT_signal_1000_500");
	TH1D *p_t5wg_HT_signal_1800_1000= (TH1D*)file_t5->Get("p_t5wg_HT_signal_1800_1000");
	TH1D *p_t5wg_HT_signal_1800_200 = (TH1D*)file_t5->Get("p_t5wg_HT_signal_1800_200");
	TH1D *p_t5wg_PhoEt_signal_1000_500 = (TH1D*)file_t5->Get("p_t5wg_PhoEt_signal_1000_500_1");
	TH1D *p_t5wg_PhoEt_signal_1800_1000= (TH1D*)file_t5->Get("p_t5wg_PhoEt_signal_1800_1000_1");
	TH1D *p_t5wg_PhoEt_signal_1800_200 = (TH1D*)file_t5->Get("p_t5wg_PhoEt_signal_1800_200_1");

	TH1D *p_tchiwg_MET_valid_500  = (TH1D*)file_tchi->Get("p_tchiwg_MET_valid_500");
	TH1D *p_tchiwg_MET_valid_800  = (TH1D*)file_tchi->Get("p_tchiwg_MET_valid_800");
	TH1D *p_tchiwg_MET_valid_1000 = (TH1D*)file_tchi->Get("p_tchiwg_MET_valid_1000");
	TH1D *p_tchiwg_MET_signal_500 = (TH1D*)file_tchi->Get("p_tchiwg_MET_signal_500");
	TH1D *p_tchiwg_MET_signal_800 = (TH1D*)file_tchi->Get("p_tchiwg_MET_signal_800");
	TH1D *p_tchiwg_MET_signal_1000= (TH1D*)file_tchi->Get("p_tchiwg_MET_signal_1000");
	TH1D *p_tchiwg_HT_signal_500  = (TH1D*)file_tchi->Get("p_tchiwg_HT_signal_500");
	TH1D *p_tchiwg_HT_signal_800  = (TH1D*)file_tchi->Get("p_tchiwg_HT_signal_800");
	TH1D *p_tchiwg_HT_signal_1000 = (TH1D*)file_tchi->Get("p_tchiwg_HT_signal_1000");

	TGraphErrors *error_PhoEt = new TGraphErrors(44);
	TGraphErrors *error_LepPt = new TGraphErrors(46);
	TGraphErrors *error_MET = new TGraphErrors(100); 
	TGraphErrors *error_Mt = new TGraphErrors(200); 
	TGraphErrors *error_HT = new TGraphErrors(64);
	TGraphErrors *ratioerror_PhoEt = new TGraphErrors(44);
	TGraphErrors *ratioerror_LepPt = new TGraphErrors(46);
	TGraphErrors *ratioerror_MET = new TGraphErrors(100); 
	TGraphErrors *ratioerror_Mt = new TGraphErrors(200); 
	TGraphErrors *ratioerror_HT = new TGraphErrors(64);

	TH1F *p_allPhoEt = (TH1F*)file_sig->Get("p_PhoEt");
	TH1F *p_allLepPt = (TH1F*)file_sig->Get("p_LepPt");
	TH1F *p_allMET   = (TH1F*)file_sig->Get("p_MET");
	TH1F *p_allMt    = (TH1F*)file_sig->Get("p_Mt");
	TH1F *p_allHT  = (TH1F*)file_sig->Get("p_HT");
	TH1F *p_allPU    = (TH1F*)file_sig->Get("p_PU");

	p_allPhoEt->SetMarkerColor(1);
	p_allLepPt->SetMarkerColor(1);
	p_allMET->SetMarkerColor(1);
	p_allMt->SetMarkerColor(1);  
	p_allHT->SetMarkerColor(1);
	p_allPU->SetMarkerColor(1);
	
	TH1F *p_elePhoEt = (TH1F*)file_ele->Get("p_PhoEt");
	TH1F *p_eleLepPt = (TH1F*)file_ele->Get("p_LepPt");
	TH1F *p_eleMET   = (TH1F*)file_ele->Get("p_MET");
	TH1F *p_eleMt    = (TH1F*)file_ele->Get("p_Mt");
	TH1F *p_eleHT  = (TH1F*)file_ele->Get("p_HT");
	TH1F *p_elePU    = (TH1F*)file_ele->Get("p_PU");

	TH1F *p_jetPhoEt = (TH1F*)file_jet->Get("p_PhoEt");
	TH1F *p_jetLepPt = (TH1F*)file_jet->Get("p_LepPt");
	TH1F *p_jetMET   = (TH1F*)file_jet->Get("p_MET");
	TH1F *p_jetMt    = (TH1F*)file_jet->Get("p_Mt");
	TH1F *p_jetHT  = (TH1F*)file_jet->Get("p_HT");
	TH1F *p_jetPU    = (TH1F*)file_jet->Get("p_PU");

	TH1F *p_qcdPhoEt = (TH1F*)file_qcd->Get("p_PhoEt");
	TH1F *p_qcdLepPt = (TH1F*)file_qcd->Get("p_LepPt");
	TH1F *p_qcdMET   = (TH1F*)file_qcd->Get("p_MET");
	TH1F *p_qcdMt    = (TH1F*)file_qcd->Get("p_Mt");
	TH1F *p_qcdHT  = (TH1F*)file_qcd->Get("p_HT");
	TH1F *p_qcdPU    = (TH1F*)file_qcd->Get("p_PU");

	TH1F *p_VGPhoEt = (TH1F*)file_VG->Get("p_PhoEt");
	TH1F *p_VGLepPt = (TH1F*)file_VG->Get("p_LepPt");
	TH1F *p_VGMET   = (TH1F*)file_VG->Get("p_MET");
	TH1F *p_VGMt    = (TH1F*)file_VG->Get("p_Mt");
	TH1F *p_VGHT  = (TH1F*)file_VG->Get("p_HT");
	TH1F *p_VGPU    = (TH1F*)file_VG->Get("p_PU");

	TH1F *p_rarePhoEt = (TH1F*)file_rare->Get("p_PhoEt");
	TH1F *p_rareLepPt = (TH1F*)file_rare->Get("p_LepPt");
	TH1F *p_rareMET   = (TH1F*)file_rare->Get("p_MET");
	TH1F *p_rareMt    = (TH1F*)file_rare->Get("p_Mt");
	TH1F *p_rareHT  = (TH1F*)file_rare->Get("p_HT");
	TH1F *p_rarePU    = (TH1F*)file_rare->Get("p_PU");


//	p_elePhoEt->Reset();
//	p_eleLepPt->Reset();
//	p_eleMET->Reset();   
//	p_eleMt->Reset();    
//	p_eleHT->Reset();  
//	p_elePU->Reset();    
//	p_jetPhoEt->Reset();
//	p_jetLepPt->Reset();
//	p_jetMET->Reset();   
//	p_jetMt->Reset();    
//	p_jetHT->Reset();  
//	p_jetPU->Reset();    
//	p_qcdPhoEt->Reset();
//	p_qcdLepPt->Reset();
//	p_qcdMET->Reset();   
//	p_qcdMt->Reset();    
//	p_qcdHT->Reset();  
//	p_qcdPU->Reset();
//	p_rarePhoEt->Reset();
//	p_rareLepPt->Reset();
//	p_rareMET->Reset();   
//	p_rareMt->Reset();    
//	p_rareHT->Reset();  
//	p_rarePU->Reset();

	gStyle->SetOptStat(0);
	TCanvas *c_pt = new TCanvas("Photon_Pt", "Photon P_{T}",600,600);
	c_pt->cd();
	TPad *pt_pad1 = new TPad("pt_pad1", "pt_pad1", 0, 0.3, 1, 1.0);
	pt_pad1->SetBottomMargin(0.1);
	pt_pad1->Draw();  
	pt_pad1->cd();  
	gPad->SetLogy();
	p_allPhoEt->SetMinimum(0.01);
	p_allPhoEt->GetXaxis()->SetRangeUser(35,300);
	p_allPhoEt->SetLineColor(1);
	p_allPhoEt->SetMarkerStyle(20);
//	p_allPhoEt->Draw("P");
	TH1D *dumEt = (TH1D*)p_allPhoEt->Clone("dumEt");
	dumEt->Reset();
	dumEt->SetMaximum(p_allPhoEt->Integral(1,10));
	dumEt->Draw();
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
	p_elePhoEt->Add(p_rarePhoEt); // ele 2nd
	p_jetPhoEt->Add(p_elePhoEt);  // jet 3rd
	p_qcdPhoEt->Add(p_jetPhoEt);  // qcd 4th
	p_VGPhoEt->Add(p_qcdPhoEt);   // VG  5th
	p_VGPhoEt->Sumw2();
	for(int ibin(1); ibin < p_VGPhoEt->GetSize(); ibin++){
		error_PhoEt->SetPoint(ibin-1,p_VGPhoEt->GetBinCenter(ibin), p_VGPhoEt->GetBinContent(ibin));
		float prederror = p_VGPhoEt->GetBinError(ibin);
	//	prederror += p_elePhoEt->GetBinError(ibin);
	//	prederror += p_jetPhoEt->GetBinError(ibin);
	//	prederror += p_qcdPhoEt->GetBinError(ibin);
	//	prederror += p_rarePhoEt->GetBinError(ibin);
		error_PhoEt->SetPointError(ibin-1,(p_VGPhoEt->GetBinLowEdge(ibin+1)-p_VGPhoEt->GetBinLowEdge(ibin))/2,prederror);
		std::cout << p_elePhoEt->GetBinError(ibin) << " " << p_jetPhoEt->GetBinError(ibin) << " " << p_qcdPhoEt->GetBinError(ibin) << " " << p_rarePhoEt->GetBinError(ibin)  << std::endl;
		ratioerror_PhoEt->SetPoint(ibin-1,p_VGPhoEt->GetBinCenter(ibin), p_allPhoEt->GetBinContent(ibin)/p_VGPhoEt->GetBinContent(ibin));
		ratioerror_PhoEt->SetPointError(ibin-1,(p_VGPhoEt->GetBinLowEdge(ibin+1)-p_VGPhoEt->GetBinLowEdge(ibin))/2, prederror/p_VGPhoEt->GetBinContent(ibin)); 
	}
	p_VGPhoEt->Draw("hist same");
	p_qcdPhoEt->Draw("hist same");
	p_jetPhoEt->Draw("hist same");
	p_elePhoEt->Draw("hist same");
	p_rarePhoEt->Draw("hist same");
  error_PhoEt->SetFillColor(2);
  error_PhoEt->SetFillStyle(3001);
	//error_PhoEt->Draw("E2 same");
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
//	p_allPhoEt->Draw("E same");

	p_t5wg_PhoEt_signal_1000_500->SetLineColor(kBlack);
  p_t5wg_PhoEt_signal_1800_1000->SetLineColor(8);
  p_t5wg_PhoEt_signal_1800_200->SetLineColor(kCyan);
	p_t5wg_PhoEt_signal_1000_500->SetLineWidth(4);
  p_t5wg_PhoEt_signal_1800_1000->SetLineWidth(4);
  p_t5wg_PhoEt_signal_1800_200->SetLineWidth(4);
	p_t5wg_PhoEt_signal_1000_500->Draw("same");
  p_t5wg_PhoEt_signal_1800_1000->Draw("same");
  p_t5wg_PhoEt_signal_1800_200->Draw("same");

	c_pt->cd();
	TPad *pt_pad2 = new TPad("pt_pad2", "pt_pad2", 0, 0.05, 1, 0.25);
	pt_pad2->Draw();
	pt_pad2->cd();
  TLine *flatratio = new TLine(35,1,400,1);
	TH1F *ratio=(TH1F*)p_allPhoEt->Clone("transfer factor");
	ratio->SetMinimum(0);
	ratio->SetMarkerStyle(20);
	ratio->SetLineColor(kBlack);
	ratio->Divide(p_VGPhoEt);
	ratio->SetTitle("");
	ratio->GetYaxis()->SetTitle("observed/bkg");
	ratio->GetXaxis()->SetLabelFont(63);
	ratio->GetXaxis()->SetLabelSize(14);
	ratio->GetYaxis()->SetLabelFont(63);
	ratio->GetYaxis()->SetLabelSize(14);
	ratio->Draw();
	ratioerror_PhoEt->SetFillColor(2);
	ratioerror_PhoEt->SetFillStyle(3001);
	ratioerror_PhoEt->Draw("E2 same");
	ratio->Draw("same");
	flatratio->Draw("same");
	c_pt->SaveAs("SIG_mg_2016ReMiniAOD_pt.pdf");

// ******** MET ************************//
	gStyle->SetOptStat(0);
	TCanvas *c_met = new TCanvas("MET", "MET",600,600);
	c_met->cd();
//	TPad *met_pad1 = new TPad("met_pad1", "met_pad1", 0, 0.3, 1, 1.0);
//	met_pad1->SetBottomMargin(0.1);
//	met_pad1->Draw();  
//	met_pad1->cd();  
	gPad->SetLogy();
	p_allMET->GetYaxis()->SetRangeUser(1,1000000);
	//p_allMET->GetYaxis()->SetRangeUser(0.01, 1000);
	p_allMET->SetMinimum(0.001);
	p_allMET->GetXaxis()->SetRangeUser(0,1000);
	p_allMET->SetLineColor(1);
	p_allMET->SetMarkerStyle(20);
//	p_allMET->Draw("P");
	TH1D *dumMET = (TH1D*)p_allMET->Clone("dumEt");
	dumMET->Reset();
	dumMET->SetMaximum(10*p_allMET->Integral(1,10));
	dumMET->Draw();
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
		float prederror = p_VGMET->GetBinError(ibin);
	//	prederror += p_eleMET->GetBinError(ibin);
	//	prederror += p_jetMET->GetBinError(ibin);
	//	prederror += p_qcdMET->GetBinError(ibin);
	//	prederror += p_rareMET->GetBinError(ibin);
		error_MET->SetPoint(ibin-1,p_VGMET->GetBinCenter(ibin), p_VGMET->GetBinContent(ibin));
		error_MET->SetPointError(ibin-1,(p_VGMET->GetBinLowEdge(ibin+1)-p_VGMET->GetBinLowEdge(ibin))/2,prederror);
		ratioerror_MET->SetPoint(ibin-1,p_VGMET->GetBinCenter(ibin), p_allMET->GetBinContent(ibin)/p_VGMET->GetBinContent(ibin));
		ratioerror_MET->SetPointError(ibin-1,(p_VGMET->GetBinLowEdge(ibin+1)-p_VGMET->GetBinLowEdge(ibin))/2, prederror/p_VGMET->GetBinContent(ibin)); 
	}
	p_VGMET->Draw("hist same");
	p_qcdMET->Draw("hist same");
	p_jetMET->Draw("hist same");
	p_eleMET->Draw("hist same");
	p_rareMET->Draw("hist same");
  error_MET->SetFillColor(2);
  error_MET->SetFillStyle(3001);
	//error_MET->Draw("E2 same");
	//
	p_t5wg_MET_signal_1000_500->SetLineColor(kBlack);
  p_t5wg_MET_signal_1800_1000->SetLineColor(8);
  p_t5wg_MET_signal_1800_200->SetLineColor(kCyan);
	p_t5wg_MET_signal_1000_500->SetLineWidth(4);
  p_t5wg_MET_signal_1800_1000->SetLineWidth(4);
  p_t5wg_MET_signal_1800_200->SetLineWidth(4);
	p_t5wg_MET_signal_1000_500->Draw("same");
  p_t5wg_MET_signal_1800_1000->Draw("same");
  p_t5wg_MET_signal_1800_200->Draw("same");
	p_tchiwg_MET_signal_500->SetLineColor(kBlack);
	p_tchiwg_MET_signal_800->SetLineColor(8);
	p_tchiwg_MET_signal_1000->SetLineColor(kCyan);
	p_tchiwg_MET_signal_500->SetLineStyle(2);
	p_tchiwg_MET_signal_800->SetLineStyle(2);
	p_tchiwg_MET_signal_1000->SetLineStyle(2);
	p_tchiwg_MET_signal_500->SetLineWidth(4);
	p_tchiwg_MET_signal_800->SetLineWidth(4);
	p_tchiwg_MET_signal_1000->SetLineWidth(4);
	p_tchiwg_MET_signal_500->Draw("same");
	p_tchiwg_MET_signal_800->Draw("same");
	p_tchiwg_MET_signal_1000->Draw("same");

	TLegend *leg_met =  new TLegend(0.6,0.75,0.9,0.9);
	leg_met->SetFillStyle(0);
	gStyle->SetLegendBorderSize(1);
	gStyle->SetLegendFillColor(0);
	leg_met->AddEntry(p_t5wg_MET_signal_1000_500, "T5WG,1000,500");
  leg_met->AddEntry(p_t5wg_MET_signal_1800_1000,"T5WG,1800,1000");
  leg_met->AddEntry(p_t5wg_MET_signal_1800_200, "T5WG,1800,200");
	leg_met->AddEntry(p_tchiwg_MET_signal_500, "TChiWG, 500");
	leg_met->AddEntry(p_tchiwg_MET_signal_800, "TChiWG, 800");
	leg_met->AddEntry(p_tchiwg_MET_signal_1000,"TChiWG, 1000");
	leg_met->Draw("same");
//	leg_met->AddEntry(p_allMET,"observed (MET < 100 GeV)");
//	leg_met->AddEntry(p_rareMET,"t#bar{t}#gamma/WW#gamma/WZ#gamma");
//	leg_met->AddEntry(p_eleMET,"e->#gamma fake");
//	leg_met->AddEntry(p_jetMET,"j->#gamma fake");
//	leg_met->AddEntry(p_qcdMET,"j->e fake");
//	leg_met->AddEntry(p_VGMET, "WG/ZG");
//	leg_met->Draw("same");
	//p_allMET->Draw("E same");

//	c_met->cd();
//	TPad *met_pad2 = new TPad("met_pad2", "met_pad2", 0, 0.05, 1, 0.25);
//	met_pad2->Draw();
//	met_pad2->cd();
//  TLine *flatratio_met = new TLine(0,1,1000,1);
//	TH1F *ratio_met=(TH1F*)p_allMET->Clone("transfer factor");
//	ratio_met->GetXaxis()->SetRangeUser(0,1000);
//	ratio_met->SetLineColor(kBlack);
//	ratio_met->SetMarkerStyle(20);
//	ratio_met->Divide(p_VGMET);
//	ratio_met->SetTitle("");
//	ratio_met->GetYaxis()->SetTitle("observed/bkg");
//	ratio_met->GetYaxis()->SetRangeUser(0,2);
//	ratio_met->GetXaxis()->SetLabelFont(63);
//	ratio_met->GetXaxis()->SetLabelSize(14);
//	ratio_met->GetYaxis()->SetLabelFont(63);
//	ratio_met->GetYaxis()->SetLabelSize(14);
//	ratio_met->Draw();
//	ratioerror_MET->SetFillColor(2);
//	ratioerror_MET->SetFillStyle(3001);
//	ratioerror_MET->Draw("E2 same");
//	ratio_met->Draw("same");
//	flatratio_met->Draw("same");
	c_met->SaveAs("SIG_mg_2016ReMiniAOD_met.pdf");


// ******** LepPt ************************//
	gStyle->SetOptStat(0);
	TCanvas *c_leppt = new TCanvas("LepPt", "LepPt",600,600);
	c_leppt->cd();
	TPad *leppt_pad1 = new TPad("leppt_pad1", "leppt_pad1", 0, 0.3, 1, 1.0);
	leppt_pad1->SetBottomMargin(0.1);
	leppt_pad1->Draw();  
	leppt_pad1->cd();  
	gPad->SetLogy();
	p_allLepPt->SetMinimum(0.001);
	p_allLepPt->GetXaxis()->SetRangeUser(35,300);
	p_allLepPt->SetLineColor(1);
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
		float prederror = p_VGLepPt->GetBinError(ibin);
//		prederror += p_eleLepPt->GetBinError(ibin);
//		prederror += p_jetLepPt->GetBinError(ibin);
//		prederror += p_qcdLepPt->GetBinError(ibin);
//		prederror += p_rareLepPt->GetBinError(ibin);
		error_LepPt->SetPointError(ibin-1,(p_VGLepPt->GetBinLowEdge(ibin+1)-p_VGLepPt->GetBinLowEdge(ibin))/2,prederror);
		ratioerror_LepPt->SetPoint(ibin-1,p_VGLepPt->GetBinCenter(ibin), p_allLepPt->GetBinContent(ibin)/p_VGLepPt->GetBinContent(ibin));
		ratioerror_LepPt->SetPointError(ibin-1,(p_VGLepPt->GetBinLowEdge(ibin+1)-p_VGLepPt->GetBinLowEdge(ibin))/2, prederror/p_VGLepPt->GetBinContent(ibin)); 
	}
	p_VGLepPt->Draw("hist same");
	p_qcdLepPt->Draw("hist same");
	p_jetLepPt->Draw("hist same");
	p_eleLepPt->Draw("hist same");
	p_rareLepPt->Draw("hist same");
  error_LepPt->SetFillColor(2);
  error_LepPt->SetFillStyle(3001);
	//error_LepPt->Draw("E2 same");
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
	TH1F *ratio_leppt=(TH1F*)p_allLepPt->Clone("transfer factor");
	ratio_leppt->SetMarkerStyle(20);
	ratio_leppt->SetLineColor(kBlack);
	ratio_leppt->Divide(p_VGLepPt);
	ratio_leppt->SetTitle("");
	ratio_leppt->GetYaxis()->SetTitle("observed/bkg");
	ratio_leppt->GetYaxis()->SetRangeUser(0,2);
	ratio_leppt->GetXaxis()->SetLabelFont(63);
	ratio_leppt->GetXaxis()->SetLabelSize(14);
	ratio_leppt->GetYaxis()->SetLabelFont(63);
	ratio_leppt->GetYaxis()->SetLabelSize(14);
	ratio_leppt->Draw();
	ratioerror_LepPt->SetFillColor(2);
	ratioerror_LepPt->SetFillStyle(3001);
	ratioerror_LepPt->Draw("E2 same");
	ratio_leppt->Draw("same");
	flatratio_leppt->Draw("same");
	c_leppt->SaveAs("SIG_mg_2016ReMiniAOD_leppt.pdf");
// ******** Mt ************************//
	gStyle->SetOptStat(0);
	TCanvas *c_mt = new TCanvas("Mt", "Mt",600,600);
	c_mt->cd();
	TPad *mt_pad1 = new TPad("mt_pad1", "mt_pad1", 0, 0.3, 1, 1.0);
	mt_pad1->SetBottomMargin(0.1);
	mt_pad1->Draw();  
	mt_pad1->cd();  
	gPad->SetLogy();
	p_allMt->SetMinimum(0.001);
  p_allMt->SetMaximum(1000000);
	p_allMt->GetXaxis()->SetRangeUser(0,1000);
	p_allMt->SetLineColor(1);
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
		float prederror = p_VGMt->GetBinError(ibin);
//		prederror += p_eleMt->GetBinError(ibin);
//		prederror += p_jetMt->GetBinError(ibin);
//		prederror += p_qcdMt->GetBinError(ibin);
//		prederror += p_rareMt->GetBinError(ibin);
		error_Mt->SetPointError(ibin-1,(p_VGMt->GetBinLowEdge(ibin+1)-p_VGMt->GetBinLowEdge(ibin))/2,prederror);
		ratioerror_Mt->SetPoint(ibin-1,p_VGMt->GetBinCenter(ibin), p_allMt->GetBinContent(ibin)/p_VGMt->GetBinContent(ibin));
		ratioerror_Mt->SetPointError(ibin-1,(p_VGMt->GetBinLowEdge(ibin+1)-p_VGMt->GetBinLowEdge(ibin))/2, prederror/p_VGMt->GetBinContent(ibin)); 
	}
	p_VGMt->Draw("hist same");
	p_qcdMt->Draw("hist same");
	p_jetMt->Draw("hist same");
	p_eleMt->Draw("hist same");
	p_rareMt->Draw("hist same");
  error_Mt->SetFillColor(2);
  error_Mt->SetFillStyle(3001);
	//error_Mt->Draw("E2 same");
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
  TLine *flatratio_mt = new TLine(0,1,100,1);
	TH1F *ratio_mt=(TH1F*)p_allMt->Clone("transfer factor");
	ratio_mt->SetMarkerStyle(20);
	ratio_mt->SetLineColor(kBlack);
	ratio_mt->GetXaxis()->SetRangeUser(0,100);
	ratio_mt->GetYaxis()->SetRangeUser(0,2);
	ratio_mt->SetMinimum(0);
	ratio_mt->SetMaximum(2);
	ratio_mt->Divide(p_VGMt);
	ratio_mt->SetTitle("");
	ratio_mt->GetYaxis()->SetTitle("observed/bkg");
	ratio_mt->GetXaxis()->SetLabelFont(63);
	ratio_mt->GetXaxis()->SetLabelSize(14);
	ratio_mt->GetYaxis()->SetLabelFont(63);
	ratio_mt->GetYaxis()->SetLabelSize(14);
	ratio_mt->Draw();
	ratioerror_Mt->SetFillColor(2);
	ratioerror_Mt->SetFillStyle(3001);
	ratioerror_Mt->Draw("E2 same");
	ratio_mt->Draw("same");
	flatratio_mt->Draw("same");
	c_mt->SaveAs("SIG_mg_2016ReMiniAOD_mt.pdf");


TCanvas *check = new TCanvas("check","",1200,800);
check->Divide(3,2);
check->cd(1);
p_allHT->Draw();
check->cd(2);
p_VGHT->Draw();
check->cd(3);
p_rareHT->Draw();
check->cd(4);
p_qcdHT->Draw();
check->cd(5);
p_eleHT->Draw();

// ******** HT ************************//
	gStyle->SetOptStat(0);
	TCanvas *c_HT = new TCanvas("HT", "HT",600,600);
	c_HT->cd();
	gPad->SetLogy();
//	TPad *HT_pad1 = new TPad("HT_pad1", "HT_pad1", 0, 0.3, 1, 1.0);
//	HT_pad1->SetBottomMargin(0.1);
//	HT_pad1->Draw();  
//	HT_pad1->cd();  
	p_allHT->SetMinimum(1);
	p_allHT->SetLineColor(1);
	p_allHT->SetMarkerStyle(20);
//	p_allHT->Draw("P");
	TH1D *dumHT = (TH1D*)p_allHT->Clone("dumEt");
	dumHT->Reset();
	dumHT->SetMaximum(10*p_allHT->Integral(1,10));
	dumHT->SetMinimum(0.01);
	dumHT->Draw();
	p_VGHT->SetFillStyle(1001);
	p_VGHT->SetLineColor(kMagenta);
	p_VGHT->SetFillColor(kMagenta);
	p_rareHT->SetFillStyle(1001);
	p_rareHT->SetLineColor(kYellow-4);
	p_rareHT->SetFillColor(kYellow-4);
	p_qcdHT->SetFillStyle(1001);
	p_qcdHT->SetLineColor(kBlue);
	p_qcdHT->SetFillColor(kBlue);
	p_eleHT->SetFillStyle(1001);
	p_eleHT->SetLineColor(kRed);
	p_eleHT->SetFillColor(kRed);
	p_jetHT->SetFillStyle(1001);
	p_jetHT->SetLineColor(kGreen);
	p_jetHT->SetFillColor(kGreen);
	p_eleHT->Add(p_rareHT); // ele 2nd
	p_jetHT->Add(p_eleHT);  // jet 3rd
	p_qcdHT->Add(p_jetHT);  // qcd 4th
	p_VGHT->Add(p_qcdHT);   // VG  5th
	p_VGHT->Draw("hist same");
	p_qcdHT->Draw("hist same");
	p_jetHT->Draw("hist same");
	p_eleHT->Draw("hist same");
	p_rareHT->Draw("hist same");


	p_t5wg_HT_signal_1000_500->SetLineColor(kBlack);
  p_t5wg_HT_signal_1800_1000->SetLineColor(8);
  p_t5wg_HT_signal_1800_200->SetLineColor(kCyan);
	p_t5wg_HT_signal_1000_500->SetLineWidth(4);
  p_t5wg_HT_signal_1800_1000->SetLineWidth(4);
  p_t5wg_HT_signal_1800_200->SetLineWidth(4);
	p_t5wg_HT_signal_1000_500->Draw("same");
  p_t5wg_HT_signal_1800_1000->Draw("same");
  p_t5wg_HT_signal_1800_200->Draw("same");
	p_tchiwg_HT_signal_500->SetLineColor(kBlack);
	p_tchiwg_HT_signal_800->SetLineColor(8);
	p_tchiwg_HT_signal_1000->SetLineColor(kCyan);
	p_tchiwg_HT_signal_500->SetLineStyle(2);
	p_tchiwg_HT_signal_800->SetLineStyle(2);
	p_tchiwg_HT_signal_1000->SetLineStyle(2);
	p_tchiwg_HT_signal_500->SetLineWidth(4);
	p_tchiwg_HT_signal_800->SetLineWidth(4);
	p_tchiwg_HT_signal_1000->SetLineWidth(4);
	p_tchiwg_HT_signal_500->Draw("same");
	p_tchiwg_HT_signal_800->Draw("same");
	p_tchiwg_HT_signal_1000->Draw("same");

	TLegend *leg_HT =  new TLegend(0.6,0.75,0.9,0.9);
	leg_HT->SetFillStyle(0);
	gStyle->SetLegendBorderSize(1);
	gStyle->SetLegendFillColor(0);
//	leg_HT->AddEntry(p_allHT,"observed (HT < 100 GeV)");
//	leg_HT->AddEntry(p_rareHT,"t#bar{t}#gamma/WW#gamma/WZ#gamma");
//	leg_HT->AddEntry(p_eleHT,"e->#gamma fake");
//	leg_HT->AddEntry(p_jetHT,"j->#gamma fake");
//	leg_HT->AddEntry(p_qcdHT,"j->e fake");
//	leg_HT->AddEntry(p_VGHT, "WG/ZG");
	leg_met->Draw("same");
//	p_allHT->Draw("E same");

//	c_HT->cd();
//	TPad *HT_pad2 = new TPad("HT_pad2", "HT_pad2", 0, 0.05, 1, 0.25);
//	HT_pad2->Draw();
//	HT_pad2->cd();
//  TLine *flatratio_HT = new TLine(0,1,400,1);
//	TH1F *ratio_HT=(TH1F*)p_VGHT->Clone("transfer factor");
//	ratio_HT->SetMarkerStyle(20);
//	ratio_HT->SetLineColor(kBlack);
//	ratio_HT->GetXaxis()->SetRangeUser(0,400);
//	ratio_HT->SetMinimum(0);
//	ratio_HT->SetMaximum(2);
//	ratio_HT->Divide(p_allHT);
//	ratio_HT->SetTitle("");
//	ratio_HT->GetYaxis()->SetTitle("observed/bkg");
//	ratio_HT->GetXaxis()->SetLabelFont(63);
//	ratio_HT->GetXaxis()->SetLabelSize(14);
//	ratio_HT->GetYaxis()->SetLabelFont(63);
//	ratio_HT->GetYaxis()->SetLabelSize(14);
//	ratio_HT->Draw();
//	flatratio_HT->Draw("same");
	c_HT->SaveAs("SIG_mg_2016ReMiniAOD_HT.pdf");

}


