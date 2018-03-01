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
#include "../../include/tdrstyle.C"

int channel = 2; // 1 = eg, 2 = mg
int plottype = 2; // 1 = bkg, 2 = valid

void plot_bkg(){//main  

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");
	setTDRStyle();

	std::ostringstream signame;  signame.str("");
	std::ostringstream VGname;   VGname.str("");

	if(channel == 1){
			signame << "validTree_egamma_VGBkg.root";
			VGname  << "NLOvalidTree_egamma_VGBkg.root";
	}
	else if(channel == 2){
			signame << "validTree_mg_VGBkg.root";
			VGname  << "NLOvalidTree_mg_VGBkg.root";
	}

	TFile *file_sig = TFile::Open(signame.str().c_str());
	TFile *file_VG  = TFile::Open(VGname.str().c_str());

	TH1F *p_allPhoEt = (TH1F*)file_sig->Get("p_PhoEt");
	TH1F *p_allLepPt = (TH1F*)file_sig->Get("p_LepPt");
	TH1F *p_allMET   = (TH1F*)file_sig->Get("p_MET");
	TH1F *p_allMt    = (TH1F*)file_sig->Get("p_Mt");
	TH1F *p_allHT  = (TH1F*)file_sig->Get("p_HT");
	TH1F *p_allISRJet    = (TH1F*)file_sig->Get("p_ISRJet");

	p_allPhoEt->SetMarkerColor(1);
	p_allLepPt->SetMarkerColor(1);
	p_allMET->SetMarkerColor(1);
	p_allMt->SetMarkerColor(1);  
	p_allHT->SetMarkerColor(1);
	p_allISRJet->SetMarkerColor(1);

	TH1F *p_VGPhoEt = (TH1F*)file_VG->Get("p_PhoEt");
	TH1F *p_VGLepPt = (TH1F*)file_VG->Get("p_LepPt");
	TH1F *p_VGMET   = (TH1F*)file_VG->Get("p_MET");
	TH1F *p_VGMt    = (TH1F*)file_VG->Get("p_Mt");
	TH1F *p_VGHT  = (TH1F*)file_VG->Get("p_HT");
	TH1F *p_VGISRJet    = (TH1F*)file_VG->Get("p_ISRJet");

	std::ostringstream etplot;  etplot.str("");
	std::ostringstream ptplot;  ptplot.str("");
	std::ostringstream metplot; metplot.str("");
	std::ostringstream mtplot;  mtplot.str("");
	std::ostringstream htplot;  htplot.str("");
	
	if(channel == 1){
			etplot << "NLOvsLO_egamma_2016ReMiniAOD_pt.pdf";
			ptplot << "NLOvsLO_egamma_2016ReMiniAOD_leppt.pdf";		
			metplot << "NLOvsLO_egamma_2016ReMiniAOD_met.pdf";
			mtplot << "NLOvsLO_egamma_2016ReMiniAOD_mt.pdf";
			htplot << "NLOvsLO_egamma_2016ReMiniAOD_ht.pdf";	
	}
	else if(channel == 2){
			etplot << "NLOvsLO_mg_2016ReMiniAOD_pt.pdf";
			ptplot << "NLOvsLO_mg_2016ReMiniAOD_leppt.pdf";		
			metplot << "NLOvsLO_mg_2016ReMiniAOD_met.pdf";
			mtplot << "NLOvsLO_mg_2016ReMiniAOD_mt.pdf";
			htplot << "NLOvsLO_mg_2016ReMiniAOD_ht.pdf";	
	}
	// ******** Mt ************************//
	gStyle->SetOptStat(0);
	TCanvas *c_mt = new TCanvas("Mt", "Mt",600,600);
	setCanvas(c_mt); 
	c_mt->cd();
	TPad *mt_pad1 = new TPad("mt_pad1", "mt_pad1", 0, 0.3, 1, 1.0);
	setTopPad(mt_pad1); 
	mt_pad1->SetBottomMargin(0);
	mt_pad1->Draw();  
	mt_pad1->cd();  
	gPad->SetLogy();
	p_allMt->SetMinimum(0.5);
  p_allMt->SetMaximum(10*p_allMt->GetBinContent(p_allMt->GetMaximumBin()));
	p_allMt->GetXaxis()->SetRangeUser(0,1000);
	p_allMt->SetLineColor(1);
	p_allMt->SetMarkerStyle(20);
	p_allMt->Draw("P");
	p_VGMt->SetLineColor(kMagenta);
	p_VGMt->Draw("hist same");
	TLegend *leg_mt =  new TLegend(0.4,0.6,0.9,0.9);
	leg_mt->SetNColumns(2);
	leg_mt->SetFillStyle(0);
	leg_mt->SetBorderSize(0);
	leg_mt->SetFillColor(0);
	leg_mt->AddEntry(p_allMt,"LO W+#gamma");
	leg_mt->AddEntry(p_VGMt, "NLO W+#gamma");
	leg_mt->Draw("same");
	p_allMt->Draw("E same");

	c_mt->cd();
	TPad *mt_pad2 = new TPad("mt_pad2", "mt_pad2", 0, 0, 1, 0.3);
	mt_pad2->SetTopMargin(0);
	mt_pad2->SetBottomMargin(0.4);
	mt_pad2->Draw();
	mt_pad2->cd();
  TLine *flatratio_mt = new TLine(0,1,1000,1);
	TH1F *ratio_mt=(TH1F*)p_allMt->Clone("transfer factor");
	ratio_mt->SetMarkerStyle(20);
	ratio_mt->SetLineColor(kBlack);
	ratio_mt->GetXaxis()->SetRangeUser(0,800);
	ratio_mt->GetYaxis()->SetRangeUser(0,2.1);
	ratio_mt->SetMinimum(0);
	ratio_mt->SetMaximum(2);
	ratio_mt->Divide(p_VGMt);
	ratio_mt->SetTitle("");
	ratio_mt->GetYaxis()->SetTitle("LO/NLO");
	ratio_mt->GetXaxis()->SetLabelFont(63);
	ratio_mt->GetXaxis()->SetLabelSize(14);
	ratio_mt->GetYaxis()->SetLabelFont(63);
	ratio_mt->GetYaxis()->SetLabelSize(14);
	ratio_mt->Draw();
	ratio_mt->Draw("same");
	flatratio_mt->Draw("same");
	c_mt->SaveAs(mtplot.str().c_str());
 
  gStyle->SetOptStat(0);
	TCanvas *c_pt = new TCanvas("Photon_Pt", "Photon P_{T}",600,600);
	setCanvas(c_pt); 
	c_pt->cd();
	TPad *pt_pad1 = new TPad("pt_pad1", "pt_pad1", 0, 0.3, 1, 1.0);
	setTopPad(pt_pad1); 
	pt_pad1->SetBottomMargin(0);
	pt_pad1->Draw();  
	pt_pad1->cd();  
	gPad->SetLogy();
	p_allPhoEt->SetTitle("p_{T}^{#gamma}");
	p_allPhoEt->SetMaximum(10*p_allPhoEt->GetBinContent(p_allPhoEt->GetMaximumBin()));
	p_allPhoEt->SetMinimum(0.5);
	p_allPhoEt->GetXaxis()->SetRangeUser(35,800);
	p_allPhoEt->SetLineColor(1);
	p_allPhoEt->SetMarkerStyle(20);
	p_allPhoEt->Draw("P");
	p_VGPhoEt->SetLineColor(kMagenta);
	p_VGPhoEt->Sumw2();
	p_VGPhoEt->Draw("hist same");
	leg_mt->Draw("same");
	p_allPhoEt->Draw("E same");


	c_pt->cd();
	TPad *pt_pad2 = new TPad("pt_pad2", "pt_pad2", 0, 0, 1, 0.3);
	pt_pad2->SetTopMargin(0);
	pt_pad2->SetBottomMargin(0.4);
	pt_pad2->Draw();
	pt_pad2->cd();
  TLine *flatratio = new TLine(35,1,800,1);
	TH1F *ratio=(TH1F*)p_allPhoEt->Clone("transfer factor");
	ratio->SetMinimum(0);
	ratio->SetMaximum(2);
	ratio->SetMarkerStyle(20);
	ratio->SetLineColor(kBlack);
	ratio->Divide(p_VGPhoEt);
	ratio->SetTitle("");
	ratio->GetYaxis()->SetTitle("LO/NLO");
	ratio->GetXaxis()->SetLabelFont(63);
	ratio->GetXaxis()->SetLabelSize(14);
	ratio->GetYaxis()->SetLabelFont(63);
	ratio->GetYaxis()->SetLabelSize(14);
	ratio->Draw();
	ratio->Draw("same");
	flatratio->Draw("same");
	c_pt->SaveAs(etplot.str().c_str());

// ******** MET ************************//
	gStyle->SetOptStat(0);
	TCanvas *c_met = new TCanvas("MET", "MET",600,600);
	setCanvas(c_met); 
	c_met->cd();
	TPad *met_pad1 = new TPad("met_pad1", "met_pad1", 0, 0.3, 1, 1.0);
	setTopPad(met_pad1); 
	met_pad1->SetBottomMargin(0);
	met_pad1->Draw();  
	met_pad1->cd();  
	gPad->SetLogy();
	p_allMET->GetYaxis()->SetRangeUser(1, 10*p_allMET->GetBinContent(p_allMET->GetMaximumBin()));
	p_allMET->SetMinimum(0.5);
	p_allMET->GetXaxis()->SetRangeUser(0,1000);
	p_allMET->SetLineColor(1);
	p_allMET->SetMarkerStyle(20);
	p_allMET->Draw("P");
	p_VGMET->SetLineColor(kMagenta);
	p_VGMET->Draw("hist same");
	leg_mt->Draw("same");
	p_allMET->Draw("E same");

	std::cout << "scale " << p_allMET->Integral(1,20)/p_VGMET->Integral(1,20) << std::endl;
	for(unsigned i(1); i < 20; i++)std::cout << p_allMET->GetBinContent(i)/p_VGMET->GetBinContent(i) << std::endl;
	
	c_met->cd();
	TPad *met_pad2 = new TPad("met_pad2", "met_pad2", 0, 0, 1, 0.3);
	met_pad2->SetTopMargin(0);
	met_pad2->SetBottomMargin(0.4);
	met_pad2->Draw();
	met_pad2->cd();
  TLine *flatratio_met = new TLine(0,1,1000,1);
	TH1F *ratio_met=(TH1F*)p_allMET->Clone("transfer factor");
	ratio_met->GetXaxis()->SetRangeUser(0,1000);
	ratio_met->SetLineColor(kBlack);
	ratio_met->SetMarkerStyle(20);
	ratio_met->Divide(p_VGMET);
	ratio_met->SetTitle("");
	ratio_met->GetYaxis()->SetTitle("LO/NLO");
	ratio_met->GetYaxis()->SetRangeUser(0,2.1);
	ratio_met->GetXaxis()->SetLabelFont(63);
	ratio_met->GetXaxis()->SetLabelSize(14);
	ratio_met->GetYaxis()->SetLabelFont(63);
	ratio_met->GetYaxis()->SetLabelSize(14);
	ratio_met->Draw();
	ratio_met->Draw("same");
	flatratio_met->Draw("same");
	c_met->SaveAs(metplot.str().c_str());


// ******** LepPt ************************//
	gStyle->SetOptStat(0);
	TCanvas *c_leppt = new TCanvas("LepPt", "LepPt",600,600);
	setCanvas(c_leppt); 
	c_leppt->cd();
	TPad *leppt_pad1 = new TPad("leppt_pad1", "leppt_pad1", 0, 0.3, 1, 1.0);
	leppt_pad1->SetBottomMargin(0);
	leppt_pad1->Draw();  
	leppt_pad1->cd();  
	gPad->SetLogy();
	p_allLepPt->GetXaxis()->SetRangeUser(35,800);
	p_allLepPt->SetLineColor(1);
	p_allLepPt->SetMarkerStyle(20);
	p_allLepPt->Draw("P");
	p_VGLepPt->SetLineColor(kMagenta);
	p_VGLepPt->Draw("hist same");

	c_leppt->cd();
	TPad *leppt_pad2 = new TPad("leppt_pad2", "leppt_pad2", 0, 0.05, 1, 0.25);
	leppt_pad2->Draw();
	leppt_pad2->cd();
  TLine *flatratio_leppt = new TLine(25,1,800,1);
	TH1F *ratio_leppt=(TH1F*)p_allLepPt->Clone("transfer factor");
	ratio_leppt->SetMarkerStyle(20);
	ratio_leppt->SetLineColor(kBlack);
	ratio_leppt->Divide(p_VGLepPt);
	ratio_leppt->SetTitle("");
	ratio_leppt->GetYaxis()->SetTitle("LO/NLO");
	ratio_leppt->GetYaxis()->SetRangeUser(0,2.1);
	ratio_leppt->GetXaxis()->SetLabelFont(63);
	ratio_leppt->GetXaxis()->SetLabelSize(14);
	ratio_leppt->GetYaxis()->SetLabelFont(63);
	ratio_leppt->GetYaxis()->SetLabelSize(14);
	ratio_leppt->Draw();
	flatratio_leppt->Draw("same");
	c_leppt->SaveAs(ptplot.str().c_str());

// ******** HT ************************//
	gStyle->SetOptStat(0);
	TCanvas *c_HT = new TCanvas("HT", "HT",600,600);
	setCanvas(c_HT); 
	c_HT->cd();
	TPad *HT_pad1 = new TPad("HT_pad1", "HT_pad1", 0, 0.3, 1, 1.0);
	setTopPad(HT_pad1); 
	HT_pad1->SetBottomMargin(0);
	HT_pad1->Draw();  
	HT_pad1->cd();  
	gPad->SetLogy();
	p_allHT->SetMinimum(5);
	p_allHT->SetMaximum(10*p_allHT->GetBinContent(p_allHT->GetMaximumBin()));
	p_allHT->SetLineColor(1);
	p_allHT->SetMarkerStyle(20);
	p_allHT->Draw("P");
	p_VGHT->SetLineColor(kMagenta);
	p_VGHT->Draw("hist same");
	leg_mt->Draw("same");
	p_allHT->Draw("E same");

	c_HT->cd();
	TPad *HT_pad2 = new TPad("HT_pad2", "HT_pad2", 0, 0, 1, 0.3);
	HT_pad2->SetTopMargin(0);
	HT_pad2->SetBottomMargin(0.4);
	HT_pad2->Draw();
	HT_pad2->cd();
  TLine *flatratio_HT = new TLine(0,1,1000,1);
	TH1F *ratio_HT=(TH1F*)p_allHT->Clone("transfer factor");
	ratio_HT->SetMarkerStyle(20);
	ratio_HT->SetLineColor(kBlack);
	ratio_HT->GetXaxis()->SetRangeUser(0,1000);
	ratio_HT->GetYaxis()->SetRangeUser(0,2.1);
	ratio_HT->SetMinimum(0);
	ratio_HT->SetMaximum(2);
	ratio_HT->Divide(p_VGHT);
	ratio_HT->SetTitle("");
	ratio_HT->GetYaxis()->SetTitle("LO/NLO");
	ratio_HT->GetXaxis()->SetLabelFont(63);
	ratio_HT->GetXaxis()->SetLabelSize(14);
	ratio_HT->GetYaxis()->SetLabelFont(63);
	ratio_HT->GetYaxis()->SetLabelSize(14);
	ratio_HT->Draw();
	c_HT->SaveAs(htplot.str().c_str());


 //******** ISRJet ************************//
	p_allISRJet->Rebin(5);
	p_VGISRJet->Rebin(5);

	gStyle->SetOptStat(0);
	TCanvas *c_ISRJet = new TCanvas("ISRJet", "ISRJet",600,600);
	setCanvas(c_ISRJet); 
	c_ISRJet->cd();
	TPad *ISRJet_pad1 = new TPad("ISRJet_pad1", "ISRJet_pad1", 0, 0.3, 1, 1.0);
	setTopPad(ISRJet_pad1); 
	ISRJet_pad1->SetBottomMargin(0);
	ISRJet_pad1->Draw();  
	ISRJet_pad1->cd();  
	gPad->SetLogy();
	p_allISRJet->SetMinimum(0.1);
	p_allISRJet->SetMaximum(10*p_allISRJet->GetBinContent(p_allISRJet->GetMaximumBin()));
	p_allISRJet->SetLineColor(1);
	p_allISRJet->SetMarkerStyle(20);
	p_allISRJet->Draw("P");
	p_VGISRJet->SetLineColor(kMagenta);
	p_VGISRJet->Draw("hist same");
	leg_mt->Draw("same");
	p_allISRJet->Draw("E same");

	c_ISRJet->cd();
	TPad *ISRJet_pad2 = new TPad("ISRJet_pad2", "ISRJet_pad2", 0, 0, 1, 0.3);
	ISRJet_pad2->SetTopMargin(0);
	ISRJet_pad2->SetBottomMargin(0.4);
	ISRJet_pad2->Draw();
	ISRJet_pad2->cd();
  TLine *flatratio_ISRJet = new TLine(0,1,1000,1);
	TH1F *ratio_ISRJet=(TH1F*)p_allISRJet->Clone("transfer factor");
	ratio_ISRJet->SetMarkerStyle(20);
	ratio_ISRJet->SetLineColor(kBlack);
	ratio_ISRJet->GetXaxis()->SetRangeUser(0,1000);
	ratio_ISRJet->GetYaxis()->SetRangeUser(0,2.1);
	ratio_ISRJet->SetMinimum(0);
	ratio_ISRJet->SetMaximum(2);
	ratio_ISRJet->Divide(p_VGISRJet);
	ratio_ISRJet->SetTitle("");
	ratio_ISRJet->GetYaxis()->SetTitle("LO/NLO");
	ratio_ISRJet->GetXaxis()->SetLabelFont(63);
	ratio_ISRJet->GetXaxis()->SetLabelSize(14);
	ratio_ISRJet->GetYaxis()->SetLabelFont(63);
	ratio_ISRJet->GetYaxis()->SetLabelSize(14);
	ratio_ISRJet->Draw();
	flatratio_ISRJet->Draw("same");
	c_ISRJet->SaveAs(htplot.str().c_str());
}


