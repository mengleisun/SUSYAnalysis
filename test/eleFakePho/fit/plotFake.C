#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TSystem.h"
#include "TChain.h"
#include "TTree.h"
#include "../../../include/tdrstyle.C"
#include "../../../include/analysis_fakes.h"

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
#include "TFitResult.h"
#include "../../../include/RooCMSShape.h"
#include "../../../include/RooDCBShape.h"
#include "../../../include/RooUserPoly.h"

void plotFake(){//main 

	bool DoPlotPt(true);
	bool DoPlotPolEta(false);
	bool DoPlotEta(true);
	bool DoPlotVtx(true);
 
	gStyle->SetOptStat(0);
	int PtBins[]={30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,65,70,75,80,90,100,120,150,180};
	unsigned nPtBins = sizeof(PtBins)/sizeof(int);
	float EtaBins[30];
	for(int i(0); i<29; i++)EtaBins[i]=0.05*i;
	unsigned nEtaBins = sizeof(EtaBins)/sizeof(float);
	float VtxBins[]={0,4,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46};
	unsigned nVtxBins = sizeof(VtxBins)/sizeof(float);
	TGraphErrors *fr_bothcount_pt 		= new TGraphErrors(nPtBins-1);
	TGraphErrors *fr_bothcount_pt_pol = new TGraphErrors(nPtBins-1);
	TGraphErrors *fr_bothcount_eta		= new TGraphErrors(nEtaBins-1);
	TGraphErrors *fr_bothcount_eta_pol= new TGraphErrors(nEtaBins-1);
	TGraphErrors *fr_bothcount_vtx		= new TGraphErrors(nVtxBins-1);
	TGraphErrors *fr_bothcount_vtx_pol= new TGraphErrors(nVtxBins-1);

	std::string bintype;
	std::string numtype;
	float lowcut;
	float signal1, error1;
	float tmpsig, tmperror;
	double pt_den[nPtBins-1];
	double pt_denerror[nPtBins-1];
	double pt_num[nPtBins-1];
	double pt_numerror[nPtBins-1];

	double pt_den_pol[nPtBins-1];
	double pt_denerror_pol[nPtBins-1];
	double pt_num_pol[nPtBins-1];
	double pt_numerror_pol[nPtBins-1];

	double eta_den[nEtaBins-1];
	double eta_denerror[nEtaBins-1];
	double eta_num[nEtaBins-1];
	double eta_numerror[nEtaBins-1];
	double eta_den_pol[nEtaBins-1];
	double eta_denerror_pol[nEtaBins-1];
	double eta_num_pol[nEtaBins-1];
	double eta_numerror_pol[nEtaBins-1];

	double vtx_den[nVtxBins-1];
	double vtx_denerror[nVtxBins-1];
	double vtx_num[nVtxBins-1];
	double vtx_numerror[nVtxBins-1];

	std::ifstream Pt_file("EleFakeRate-2016ReReco-ker-pt.txt");
	std::ifstream Pt_file_pol("../ReMiniAODresult/EleFakeRate-2016ReReco-ker-pt.txt");
	std::ifstream Eta_file("../ReMiniAODresult/EleFakeRate-2016ReReco-ker-eta.txt");
	std::ifstream Eta_file_pol("EleFakeRate-2016ReReco-pol-eta.txt");
	std::ifstream Vtx_file("../ReMiniAODresult/EleFakeRate-2016ReReco-ker-vtx.txt");
	 
	if(DoPlotPt){
		if(Pt_file.is_open()){
			for(int i(0); i<nPtBins-1; i++){ 
				Pt_file >> bintype >> numtype >> lowcut >> signal1 >> error1;
				pt_den[i] = signal1;
				pt_denerror[i] = error1;
			}    
			for(int i(0); i<nPtBins-1; i++){ 
				Pt_file >> bintype >> numtype >> lowcut >> signal1 >> error1;
				pt_num[i] = signal1;
				pt_numerror[i] = error1;
			}    
		}
		if(Pt_file_pol.is_open()){
			for(int i(0); i<nPtBins-1; i++){ 
				Pt_file_pol >> bintype >> numtype >> lowcut >> signal1 >> error1;
				pt_den_pol[i] = signal1;
				pt_denerror_pol[i] = error1;
			}    
			for(int i(0); i<nPtBins-1; i++){ 
				Pt_file_pol >> bintype >> numtype >> lowcut >> signal1 >> error1;
				pt_num_pol[i] = signal1;
				pt_numerror_pol[i] = error1;
			}    
		}
		for(int i(0); i<nPtBins-1; i++){
			double syserror = (pt_num_pol[i]/(pt_den_pol[i]+pt_num_pol[i]) - pt_num[i]/(pt_den[i]+pt_num[i]) );
			double staterror = pt_num[i]*pt_denerror[i]/(pt_den[i]*pt_den[i])+pt_numerror[i]/pt_den[i];
			fr_bothcount_pt->SetPoint(i,(PtBins[i]+PtBins[i+1])/2.0,pt_num[i]/(pt_den[i]+pt_num[i]) );
			fr_bothcount_pt->SetPointError(i, (PtBins[i+1]-PtBins[i])/2.0, sqrt(syserror*syserror+staterror*staterror));
			fr_bothcount_pt_pol->SetPoint(i,(PtBins[i]+PtBins[i+1])/2.0,pt_num[i]/(pt_den[i]+pt_num[i]) );
			fr_bothcount_pt_pol->SetPointError(i, 0, syserror);
		}
	}

	if(DoPlotEta){
		if(Eta_file.is_open()){
			for(int i(0); i<nEtaBins-1; i++){
				Eta_file >> bintype >> numtype >> lowcut >> signal1 >> error1;
				eta_den[i] = signal1;
				eta_denerror[i] = error1;
			}
			for(int i(0); i<nEtaBins-1; i++){
				Eta_file >> bintype >> numtype >> lowcut >> signal1 >> error1;
				eta_num[i] = signal1;
				eta_numerror[i] = error1;
			}
		}
		for(int i(0); i<nEtaBins-1; i++){
			fr_bothcount_eta->SetPoint(i,(EtaBins[i]+EtaBins[i+1])/2.0,eta_num[i]/(eta_den[i]+eta_num[i]) );
			fr_bothcount_eta->SetPointError(i, (EtaBins[i+1]-EtaBins[i])/2.0, eta_num[i]*eta_denerror[i]/(eta_den[i]*eta_den[i])+eta_numerror[i]/eta_den[i]);
			std::cout << "eta " << EtaBins[i] << " " << eta_num[i]/(eta_den[i]+eta_num[i]) << std::endl;
		}
	}

	if(DoPlotPolEta){
		if(Eta_file_pol.is_open()){
			for(int i(0); i<nEtaBins-1; i++){
				Eta_file_pol >> bintype >> numtype >> lowcut >> signal1 >> error1;
				eta_den_pol[i] = signal1;
				eta_denerror_pol[i] = error1;
			}
			for(int i(0); i<nEtaBins-1; i++){
				Eta_file_pol >> bintype >> numtype >> lowcut >> signal1 >> error1;
				eta_num_pol[i] = signal1;
				eta_numerror_pol[i] = error1;
			}
		}
		for(int i(0); i<nEtaBins-1; i++){
			fr_bothcount_eta_pol->SetPoint(i,(EtaBins[i]+EtaBins[i+1])/2.0,eta_num_pol[i]/(eta_den_pol[i]+eta_num_pol[i]) );
			fr_bothcount_eta_pol->SetPointError(i, (EtaBins[i+1]-EtaBins[i])/2.0, eta_num_pol[i]*eta_denerror_pol[i]/(eta_den_pol[i]*eta_den_pol[i])+eta_numerror_pol[i]/eta_den_pol[i]);
		}
	}

	if(DoPlotVtx){
		if(Vtx_file.is_open()){
			for(int i(0); i<nVtxBins-1; i++){
				Vtx_file >> bintype >> numtype >> lowcut >> signal1 >> error1; 
				vtx_den[i] = signal1;
				vtx_denerror[i] = error1;
			}
			for(int i(0); i<nVtxBins-1; i++){
				Vtx_file >> bintype >> numtype >> lowcut >> signal1 >> error1; 
				vtx_num[i] = signal1;
				vtx_numerror[i] = error1;
			}
		}
		for(int i(0); i<nVtxBins-1; i++){
			fr_bothcount_vtx->SetPoint(i,(VtxBins[i]+VtxBins[i+1])/2.0,vtx_num[i]/(vtx_den[i]+vtx_num[i]) );
			fr_bothcount_vtx->SetPointError(i, (VtxBins[i+1]-VtxBins[i])/2.0, vtx_num[i]*vtx_denerror[i]/(vtx_den[i]*vtx_den[i])+vtx_numerror[i]/vtx_den[i]);
		}
	}

	setTDRStyle();
	gStyle->SetOptFit(0);
	fr_bothcount_pt->SetMarkerStyle(20);
	fr_bothcount_pt->SetMarkerColor(kBlue);
	fr_bothcount_pt->SetLineColor(kBlue);
	fr_bothcount_pt->SetLineWidth(2);
	fr_bothcount_pt->SetFillColor(0);
	fr_bothcount_eta->SetMarkerStyle(20);
	fr_bothcount_eta->SetMarkerColor(kBlue);
	fr_bothcount_eta->SetLineColor(kBlue);
	fr_bothcount_eta->SetLineWidth(2);
	fr_bothcount_eta->SetFillColor(0);
	fr_bothcount_vtx->SetMarkerStyle(20);
	fr_bothcount_vtx->SetMarkerColor(kBlue);
	fr_bothcount_vtx->SetLineColor(kBlue);
	fr_bothcount_vtx->SetLineWidth(2);
	fr_bothcount_vtx->SetFillColor(0);
	fr_bothcount_eta_pol->SetMarkerStyle(20);
	fr_bothcount_eta_pol->SetMarkerColor(kRed);
	fr_bothcount_eta_pol->SetLineColor(kRed);
	fr_bothcount_eta_pol->SetLineWidth(2);
	fr_bothcount_eta_pol->SetFillColor(0);

	TF1 *fakerate_pt;
	TF1 *fakerate_vtx = new TF1("vtxfake","pol1",0,44);
	if(DoPlotPt){
		TH1F *dummy_pt = new TH1F("",";Pt(GeV);fake rate",16,20,180);
		TCanvas *canvas_pt = new TCanvas("fake-rate vs. p_{T}","",600,600);
		canvas_pt->cd();
		dummy_pt->SetMaximum(0.04);
		dummy_pt->Draw();
		TF1 *f1 = new TF1("f1", fakerate_ptDependence,30,170,4);
		f1->SetParameters(0.00123, -0.002, -0.834, 0.002);
		f1->SetParNames("slope","constant","index","coeff");
		fr_bothcount_pt->Draw("P same");
	//	fr_bothcount_pt_pol->SetLineColor(kRed);
	//	fr_bothcount_pt_pol->SetLineWidth(2);
	//	fr_bothcount_pt_pol->Draw("E same");
		TLegend *leg=new TLegend(0.6,0.7,0.9,0.9);
		fr_bothcount_pt->Fit("f1","R");
		fakerate_pt = fr_bothcount_pt->GetFunction("f1");
		std::cout << "fitting result" << "  slope=" << fakerate_pt->GetParameter(0) << "  constant=" << fakerate_pt->GetParameter(1) << " index=" << fakerate_pt->GetParameter(2) << " coeff=" << fakerate_pt->GetParameter(3) << std::endl;
		canvas_pt->SaveAs("elefakepho_rate_pt.pdf");
	}

	if(DoPlotEta){
		TH1F *dummy_eta = new TH1F("",";#eta;fake rate",14,0,1.4);
		TCanvas *canvas_eta = new TCanvas("fake-rate vs. #eta","fake rate",600,600);
		canvas_eta->cd();
		dummy_eta->SetMaximum(0.08);
		dummy_eta->Draw();
		fr_bothcount_eta->Draw("L same");
		//fr_bothcount_eta_pol->Draw("L same");
		canvas_eta->SaveAs("elefakepho_rate_eta.pdf");
	}

	if(DoPlotVtx){
		TH1F *dummy_vtx = new TH1F("",";nVtx;fake rate",44,0,44);
		TCanvas *canvas_vtx = new TCanvas("fake-rate vs. nVtx"," fake rate",600,600);
		canvas_vtx->cd();
		dummy_vtx->SetMaximum(0.08);
		dummy_vtx->Draw();
		fr_bothcount_vtx->Draw("P same");
		fakerate_vtx->SetParameter(0,0.02);
		fakerate_vtx->SetParameter(1,1);
		fr_bothcount_vtx->Fit("vtxfake");
		fakerate_vtx->Draw("same");
		canvas_vtx->SaveAs("elefakepho_rate_vtx.pdf");
	}

}


