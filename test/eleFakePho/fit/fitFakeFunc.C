#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include "TH1D.h"
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
#include "TMatrixDSym.h"
#include "TRandom3.h"
#include "../../../include/tdrstyle.C"

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
#include "TVirtualFitter.h"
#include "RooMultiVarGaussian.h"
#include "../../../include/RooCMSShape.h"
#include "../../../include/RooDCBShape.h"
#include "../../../include/RooUserPoly.h"
#include "../../../include/analysis_fakes.h"
#define NTOY 100

bool doEB = false;

Double_t fakerate_ptDependence(Double_t *x, Double_t *par)
{
	double slope = par[0];
	double constant = par[1]; 
	double index = par[2];
	double coeff = 1.0; 
	
	double pt = TMath::Max(x[0],0.000001);
	
	double arg = 0;
	arg = slope*pt + constant; 
	double fitval = pow(arg, index)*coeff; 
	return fitval;
}

Double_t fakerate_etaDependence(Double_t *x, Double_t *par){
	Double_t weight_eta(0);

	double eta = x[0];
	if(eta >= 0 && eta < 1.4442){
		for(int ieta(0); ieta < 29; ieta++)
			if(eta > ieta*0.05 && eta <= (ieta+1)*0.05)weight_eta = etaRatesEB[ieta];
	}
	else if(eta > 1.56 && eta <= 2.5){
		for(int ieta(0); ieta < 94; ieta++)
			if(eta > 1.56 + ieta*0.01 && eta <= 1.56+(ieta+1)*0.01)weight_eta = etaRatesEE[ieta];
	}
	else weight_eta = 0;

	return weight_eta;
}
	

void fitFakeFunc(){//main 
 
	setTDRStyle();
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	gStyle->SetErrorX(0.5);
	gStyle->SetTitleX(0.5);

	int ptlowedge(30), pthighedge(200);
	int vtxlowedge(0), vtxhighedge(46);
	int graphPtBins = (int)((pthighedge - ptlowedge)/0.5);
	int graphVtxBins = (int)((vtxhighedge - vtxlowedge)/0.5);
	int PtBins[] ={30,35,40,45,50,55,60,65,70,75,80,90,100,120,150,180,200};
	unsigned nPtBins = sizeof(PtBins)/sizeof(int);
	unsigned nEtaBins(0);
	if(doEB)nEtaBins=30;
	else nEtaBins = 95;
	float EtaBins[nEtaBins];
	for(int i(0); i<nEtaBins; i++){
		if(doEB)EtaBins[i]=0.05*i;
		else EtaBins[i]=1.56 + 0.01*i;	
	}
	float VtxBins[]={0,4,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46};
	unsigned nVtxBins = sizeof(VtxBins)/sizeof(float);
	TGraphErrors *fr_bothcount_pt = new TGraphErrors(nPtBins-1);
	TGraphErrors *fr_bothcount_eta= new TGraphErrors(nEtaBins-1);
	TGraphErrors *fr_bothcount_vtx= new TGraphErrors(nVtxBins-1);
	TGraphErrors *fr_pt_sigmaband = new TGraphErrors(graphPtBins);
	TGraphErrors *fr_vtx_sigmaband = new TGraphErrors(graphVtxBins);
	TGraphErrors *fr_eta_sigmaband = new TGraphErrors(nEtaBins-1);
	TGraphErrors *fr_pt_ratio = new TGraphErrors(nVtxBins-1);
	TGraphErrors *fr_pt_ratioError = new TGraphErrors(nVtxBins-1);
	TGraphErrors *fr_vtx_ratio = new TGraphErrors(nVtxBins-1);
	TGraphErrors *fr_vtx_ratioError = new TGraphErrors(nVtxBins-1);

	std::string bintype;
	std::string numtype;
	float lowcut;
	float signal1, error1;
	float fitmean, fitrms;
	float DYmean, DYrms;
	float Polmean, Polrms;
	float tmpsig, tmperror;
	double pt_den[nPtBins-1];
	double pt_denerror[nPtBins-1];
	double pt_num[nPtBins-1];
	double pt_numerror[nPtBins-1];

	double eta_den[nEtaBins-1];
	double eta_denerror[nEtaBins-1];
	double eta_num[nEtaBins-1];
	double eta_numerror[nEtaBins-1];

	double vtx_den[nVtxBins-1];
	double vtx_denerror[nVtxBins-1];
	double vtx_num[nVtxBins-1];
	double vtx_numerror[nVtxBins-1];

//	std::ifstream Pt_file("/uscms_data/d3/mengleis/SUSYAnalysis/test/eleFakePho/Aug3/EleFakeRate-2016ReReco-Bw-ker-pt-60-120.txt");
//	std::ifstream Pt_DYfile("/uscms_data/d3/mengleis/SUSYAnalysis/test/eleFakePho/Aug3/EleFakeRate-2016ReReco-Bw-ker-pt-60-120.txt");
//	std::ifstream Pt_Polfile("/uscms_data/d3/mengleis/SUSYAnalysis/test/eleFakePho/Aug3/EleFakeRate-2016ReReco-Bw-ker-pt-60-120.txt");
//	std::ifstream Eta_file("/uscms_data/d3/mengleis/SUSYAnalysis/test/eleFakePho/Aug3/EleFakeRate-2016ReReco-Bw-ker-eta-60-120.txt");
//	std::ifstream Eta_DYfile("/uscms_data/d3/mengleis/SUSYAnalysis/test/eleFakePho/Aug3/EleFakeRate-2016ReReco-Bw-ker-eta-60-120.txt");
//	std::ifstream Eta_Polfile("/uscms_data/d3/mengleis/SUSYAnalysis/test/eleFakePho/Aug3/EleFakeRate-2016ReReco-Bw-ker-eta-60-120.txt");
//	std::ifstream Vtx_file("/uscms_data/d3/mengleis/SUSYAnalysis/test/eleFakePho/Aug3/EleFakeRate-2016ReReco-Bw-ker-vtx-60-120.txt");
//	std::ifstream Vtx_DYfile("/uscms_data/d3/mengleis/SUSYAnalysis/test/eleFakePho/Aug3/EleFakeRate-2016ReReco-Bw-ker-vtx-60-120.txt");
//	std::ifstream Vtx_Polfile("/uscms_data/d3/mengleis/SUSYAnalysis/test/eleFakePho/Aug3/EleFakeRate-2016ReReco-Bw-ker-vtx-60-120.txt");
	
	std::ifstream Pt_file("/uscms_data/d3/mengleis/SUSYAnalysis/test/eleFakePho/FullEcal/EleFakeRate-2016ReReco-Bw-ker-pt-40-140.txt");
	std::ifstream Pt_DYfile("/uscms_data/d3/mengleis/SUSYAnalysis/test/eleFakePho/FullEcal/EleFakeRate-2016ReReco-Bw-ker-pt-40-140.txt");
	std::ifstream Pt_Polfile("/uscms_data/d3/mengleis/SUSYAnalysis/test/eleFakePho/FullEcal/EleFakeRate-2016ReReco-Bw-ker-pt-40-140.txt");
	std::ifstream Eta_file("/uscms_data/d3/mengleis/SUSYAnalysis/test/eleFakePho/FullEcal/EleFakeRate-2016ReReco-Bw-ker-eta-40-140.txt");
	std::ifstream Eta_DYfile("/uscms_data/d3/mengleis/SUSYAnalysis/test/eleFakePho/FullEcal/EleFakeRate-2016ReReco-Bw-ker-eta-40-140.txt");
	std::ifstream Eta_Polfile("/uscms_data/d3/mengleis/SUSYAnalysis/test/eleFakePho/FullEcal/EleFakeRate-2016ReReco-Bw-ker-eta-40-140.txt");
	std::ifstream Vtx_file("/uscms_data/d3/mengleis/SUSYAnalysis/test/eleFakePho/FullEcal/EleFakeRate-2016ReReco-Bw-ker-vtx-40-140.txt");
	std::ifstream Vtx_DYfile("/uscms_data/d3/mengleis/SUSYAnalysis/test/eleFakePho/FullEcal/EleFakeRate-2016ReReco-Bw-ker-vtx-40-140.txt");
	std::ifstream Vtx_Polfile("/uscms_data/d3/mengleis/SUSYAnalysis/test/eleFakePho/FullEcal/EleFakeRate-2016ReReco-Bw-ker-vtx-40-140.txt");
//************************************ Fill and fit the pt dependence. Calculate ratio **********************************************************// 
	if(Pt_file.is_open() && Pt_DYfile.is_open() && Pt_Polfile.is_open()){
		for(unsigned i(0); i<nPtBins-1; i++){ 
			Pt_file >> bintype >> numtype >> lowcut >> signal1 >> error1 >> fitmean >> fitrms;
			pt_den[i] = fitmean;
			Pt_DYfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> DYmean >> DYrms;
			Pt_Polfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> Polmean >> Polrms;
			double sysdiff = fabs(DYmean - fitmean) > fabs(Polmean - fitmean)? fabs(DYmean - fitmean):fabs(Polmean - fitmean);
			pt_denerror[i] = sqrt(fitrms*fitrms + sysdiff*sysdiff);
		}    
		for(unsigned i(0); i<nPtBins-1; i++){ 
			Pt_file >> bintype >> numtype >> lowcut >> signal1 >> error1 >> fitmean >> fitrms;
			pt_num[i] = fitmean;
			Pt_DYfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> DYmean >> DYrms;
			Pt_Polfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> Polmean >> Polrms;
			double sysdiff = fabs(DYmean - fitmean) > fabs(Polmean - fitmean)? fabs(DYmean - fitmean):fabs(Polmean - fitmean);
			pt_numerror[i] = sqrt(fitrms*fitrms + sysdiff*sysdiff);
		}    
	}
	for(unsigned i(0); i<nPtBins-1; i++){
		double fakerate = pt_num[i]/pt_den[i];
		double error = sqrt(fakerate*fakerate*pt_denerror[i]*pt_denerror[i]/(pt_den[i]*pt_den[i])+ pt_numerror[i]*pt_numerror[i]/pt_den[i]/pt_den[i]);
		fr_bothcount_pt->SetPoint(i,(PtBins[i]+PtBins[i+1])/2.0, fakerate );
		fr_bothcount_pt->SetPointError(i, (PtBins[i+1]-PtBins[i])/2.0, error);
	}
	TCanvas *can=new TCanvas("can","",600,600);
	can->cd();
	fr_bothcount_pt->Draw("AP");
    	TF1 *f1 = new TF1("f1", fakerate_ptDependence,30,1000,3);
    	TF1 *fit_fakerate_pt;
    	TFitResultPtr result_fitpt;
    	TVirtualFitter::SetMaxIterations(1000000);
// 	for(int i(1); i<50; i++){
// 		for(int j(1); j <50; j++){
// 			std::cout << i << " " << j << std::endl;
// 		}
// 	}
	//*************** EB ******************//
////	f1->SetParameter(0, 20);
////	f1->SetParLimits(0, 0, 1000);
////	f1->SetParameter(1, 206);
////	f1->SetParameter(2, -1.92);
////	f1->SetParLimits(2, -4, 0);
////	f1->SetParNames("slope","constant","index");
////	result_fitpt = fr_bothcount_pt->Fit("f1","R S");
	//*************** EE ******************//
	f1->SetParameter(0, 1);
	f1->SetParLimits(0, 0, 1000);
	f1->SetParameter(1, 10);
	f1->SetParameter(2, -1.92);
	f1->SetParLimits(2, -4, 0);
	f1->SetParNames("slope","constant","index");
	fr_bothcount_pt->Fit("f1","R S");
    	result_fitpt = fr_bothcount_pt->Fit("f1","R S");
    
 	fit_fakerate_pt = fr_bothcount_pt->GetFunction("f1");
 	for(unsigned i(0); i<nPtBins-1; i++){
 		double fakerate = pt_num[i]/pt_den[i];
 		double error = sqrt(fakerate*fakerate*pt_denerror[i]*pt_denerror[i]/(pt_den[i]*pt_den[i])+ pt_numerror[i]*pt_numerror[i]/pt_den[i]/pt_den[i]);
 		fr_pt_ratio->SetPoint(i,(PtBins[i]+PtBins[i+1])/2.0, fit_fakerate_pt->Eval((PtBins[i]+PtBins[i+1])/2.0)/fakerate);
 		fr_pt_ratio->SetPointError(i, (PtBins[i+1]-PtBins[i])/2.0, error/fit_fakerate_pt->Eval((PtBins[i]+PtBins[i+1])/2.0));
 	}
 
 //********************************* Fill and fit the eta dependence ********************************************************************************//
 
 	if(Eta_file.is_open() && Eta_DYfile.is_open() && Eta_Polfile.is_open()){
 		for(unsigned i(0); i<nEtaBins-1; i++){ 
 			Eta_file >> bintype >> numtype >> lowcut >> signal1 >> error1 >> fitmean >> fitrms;
 			std::cout << lowcut << " " << signal1 << std::endl;
 			eta_den[i] = fitmean;
 			Eta_DYfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> DYmean >> DYrms;
 			Eta_Polfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> Polmean >> Polrms;
 			double sysdiff = fabs(DYmean - fitmean) > fabs(Polmean - fitmean)? fabs(DYmean - fitmean):fabs(Polmean - fitmean);
 			eta_denerror[i] = sqrt(fitrms*fitrms + sysdiff*sysdiff);
 		}    
 		for(unsigned i(0); i<nEtaBins-1; i++){ 
 			Eta_file >> bintype >> numtype >> lowcut >> signal1 >> error1 >> fitmean >> fitrms;
 			std::cout << lowcut << " " << signal1 << std::endl;
 			eta_num[i] = fitmean;
 			Eta_DYfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> DYmean >> DYrms;
 			Eta_Polfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> Polmean >> Polrms;
 			double sysdiff = fabs(DYmean - fitmean) > fabs(Polmean - fitmean)? fabs(DYmean - fitmean):fabs(Polmean - fitmean);
 			eta_numerror[i] = sqrt(fitrms*fitrms + sysdiff*sysdiff);
 		}    
 	}
 	for(unsigned i(0); i<nEtaBins-1; i++){
 		double fakerate = eta_num[i]/eta_den[i];
 		double error = sqrt(fakerate*fakerate*eta_denerror[i]*eta_denerror[i]/(eta_den[i]*eta_den[i])+ eta_numerror[i]*eta_numerror[i]/eta_den[i]/eta_den[i]);
 		std::cout << "eta " << EtaBins[i] << " " << fakerate << std::endl;
 		fr_bothcount_eta->SetPoint(i,(EtaBins[i]+EtaBins[i+1])/2.0, fakerate );
 		fr_bothcount_eta->SetPointError(i, (EtaBins[i+1]-EtaBins[i])/2.0, error);
 		fr_eta_sigmaband->SetPoint(i,(EtaBins[i]+EtaBins[i+1])/2.0, fakerate );
 		fr_eta_sigmaband->SetPointError(i, (EtaBins[i+1]-EtaBins[i])/2.0, error);
 	}
	TF1 *fit_fakerate_eta = new TF1("fit_fakerate_eta", fakerate_etaDependence, 0, 2.5, 0); 
 
 //******************************** Fill and fit vtx dependence **********************************************************************************//
 	if(Vtx_file.is_open() && Vtx_DYfile.is_open() && Vtx_Polfile.is_open()){
 		for(unsigned i(0); i<nVtxBins-1; i++){ 
 			Vtx_file >> bintype >> numtype >> lowcut >> signal1 >> error1 >> fitmean >> fitrms;
 			vtx_den[i] = fitmean;
 			Vtx_DYfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> DYmean >> DYrms;
 			Vtx_Polfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> Polmean >> Polrms;
 			double sysdiff = fabs(DYmean - fitmean) > fabs(Polmean - fitmean)? fabs(DYmean - fitmean):fabs(Polmean - fitmean);
 			vtx_denerror[i] = sqrt(fitrms*fitrms + sysdiff*sysdiff);
 		}    
 		for(unsigned i(0); i<nVtxBins-1; i++){ 
 			Vtx_file >> bintype >> numtype >> lowcut >> signal1 >> error1 >> fitmean >> fitrms;
 			vtx_num[i] = fitmean;
 			Vtx_DYfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> DYmean >> DYrms;
 			Vtx_Polfile >> bintype >> numtype >> lowcut >> signal1 >> error1 >> Polmean >> Polrms;
 			double sysdiff = fabs(DYmean - fitmean) > fabs(Polmean - fitmean)? fabs(DYmean - fitmean):fabs(Polmean - fitmean);
 			vtx_numerror[i] = sqrt(fitrms*fitrms + sysdiff*sysdiff);
 		}    
 	}
 	for(unsigned i(0); i<nVtxBins-1; i++){
 		double fakerate = vtx_num[i]/vtx_den[i];
 		double error = sqrt(fakerate*fakerate*vtx_denerror[i]*vtx_denerror[i]/(vtx_den[i]*vtx_den[i])+ vtx_numerror[i]*vtx_numerror[i]/vtx_den[i]/vtx_den[i]);
 		fr_bothcount_vtx->SetPoint(i,(VtxBins[i]+VtxBins[i+1])/2.0, fakerate );
 		fr_bothcount_vtx->SetPointError(i, (VtxBins[i+1]-VtxBins[i])/2.0, error);
 	}
 	TF1 *fit_fakerate_vtx = new TF1("vtxfake","pol1",0,44);
 	TFitResultPtr result_fitvtx;
 	fit_fakerate_vtx->SetParameter(0,0.02);
 	fit_fakerate_vtx->SetParameter(1,1);
 	result_fitvtx = fr_bothcount_vtx->Fit("vtxfake","S");
 	for(unsigned i(0); i<nVtxBins-1; i++){
 		double fakerate = vtx_num[i]/vtx_den[i];
 		double error = sqrt(fakerate*fakerate*vtx_denerror[i]*vtx_denerror[i]/(vtx_den[i]*vtx_den[i])+ vtx_numerror[i]*vtx_numerror[i]/vtx_den[i]/vtx_den[i]);
 		fr_vtx_ratio->SetPoint(i,(VtxBins[i]+VtxBins[i+1])/2.0, fit_fakerate_vtx->Eval((VtxBins[i]+VtxBins[i+1])/2.0)/fakerate);
 		fr_vtx_ratio->SetPointError(i, (VtxBins[i+1]-VtxBins[i])/2.0, error/fit_fakerate_vtx->Eval((VtxBins[i]+VtxBins[i+1])/2.0));
 	}
 
// *****************************************************************************************************************************//
 
 
	gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libRooFitClasses.so");
	
	TH1D* invmass_den = new TH1D("invmass_den", "invmass_den",100,40,140);
	TH1D* invmass_num = new TH1D("invmass_num", "invmass_num",100,40,140);
	TH1D* h_invmass_bg  = new TH1D("invmass_bg",  "invmass_bg", 100,40,140);
	TH1D* invmass_prednum = new TH1D("invmass_prednum", "invmass_prednum",100,40,140);
	
	TChain *etree = new TChain("FakeRateTree");
	etree->Add("/uscms_data/d3/mengleis/Sep1/plot_elefakepho-FullEcalTnP.root");
	float invmass=0; 
	float tagPt=0; 
	float probePt=0; 
	float probeEta=0;
	bool  vetovalue=0;
	bool  FSRveto = 0;
	int   nVertex=0;
	etree->SetBranchAddress("invmass",   &invmass); 
	etree->SetBranchAddress("tagPt",     &tagPt);
	etree->SetBranchAddress("probePt",   &probePt);
	etree->SetBranchAddress("probeEta",  &probeEta);
	etree->SetBranchAddress("vetovalue", &vetovalue);
	etree->SetBranchAddress("nVertex",   &nVertex);
	etree->SetBranchAddress("FSRveto",   &FSRveto);
	
	std::vector<float> etreeEt;
	std::vector<float> etreeVtx;
	std::vector<float> etreeEta;
	std::vector<float> etreeInvmass;
	etreeEt.clear();
	etreeVtx.clear(); 
	etreeEta.clear();
	etreeInvmass.clear();
	for(unsigned iEvt(0); iEvt < etree->GetEntries(); iEvt++){
		etree->GetEntry(iEvt);
		     	
		if(probePt < 35)continue;
		if(doEB && fabs(probeEta) > 1.4442)continue;
		else if(!doEB && (fabs(probeEta) < 1.56 || fabs(probeEta) > 2.1))continue;

		if(vetovalue== false)invmass_den->Fill(invmass);
		else if(vetovalue== true && FSRveto == true)invmass_num->Fill(invmass);
		
		float w_ele = fit_fakerate_pt->Eval(probePt)*fit_fakerate_vtx->Eval(nVertex)*fit_fakerate_eta->Eval(fabs(probeEta));
		if(w_ele > 0.5 || w_ele < 0)std::cout << "etree " << probePt << " " << probeEta << " " << w_ele << std::endl;
		if(vetovalue==false){
			invmass_prednum->Fill(invmass, w_ele);
			etreeEt.push_back(probePt);
			etreeVtx.push_back(nVertex);
			etreeEta.push_back(fabs(probeEta));
			etreeInvmass.push_back(invmass);
		} 
	}
	invmass_prednum->Sumw2();
 
	TChain *bgtree = new TChain("BGTree");
	bgtree->Add("/uscms_data/d3/mengleis/Sep1/plot_bgtemplate_FullEcal.root");
	float invmass_bg=0; 
	float probePt_bg=0; 
	float probeEta_bg=0; 
	bool  vetovalue_bg=0;
	bool  FSRveto_bg=0;	
	int   nVertex_bg=0;
	bgtree->SetBranchAddress("invmass",   &invmass_bg); 
	bgtree->SetBranchAddress("probePt",   &probePt_bg);
	bgtree->SetBranchAddress("probeEta",   &probeEta_bg);
	bgtree->SetBranchAddress("vetovalue", &vetovalue_bg);
	bgtree->SetBranchAddress("FSRveto",   &FSRveto_bg);
	bgtree->SetBranchAddress("nVertex",   &nVertex_bg);
	
	
	for(unsigned iEvt(0); iEvt < bgtree->GetEntries(); iEvt++){
		bgtree->GetEntry(iEvt);
	
		if(probePt_bg < 35)continue; 
		if(doEB && fabs(probeEta_bg) > 1.4442)continue;
		else if(!doEB && (fabs(probeEta_bg) < 1.56 || fabs(probeEta_bg) > 2.1))continue;
	     	
		if(vetovalue_bg== true && FSRveto_bg==true){h_invmass_bg->Fill(invmass_bg);}
	}
	
	RooRealVar mass_axis("invmass","invmass",70,110);
	TCanvas *c_fitMass = new TCanvas("c_fitMass", "", 600, 600);
	c_fitMass->cd();
 
	RooDataHist datahist_data("both", "", mass_axis, invmass_num);
	RooDataHist datahist_bg("bg","",mass_axis, h_invmass_bg);
	RooHistPdf  pdf_bg("pdf_bg","pdf_bg",mass_axis, datahist_bg);
 
	RooRealVar m0( "m0", "m0", 91.188);
	RooRealVar width( "width", "width", 2.495,0,4);
	RooRealVar mean("mean", "" ,0.);
	RooRealVar sigma("sigma", "",2.4 , 0.0, 15.0);
	RooRealVar alpha("alpha", "", 1.0, 0.0, 20.0);
	RooRealVar n("n","", 1.0, 0.0, 20.0);
	RooRealVar alpha2("2ndalpha","", 1.0, 0.0, 20.0);
	RooRealVar n2("2ndn", "", 1.0, 0.0, 20.0);
	RooBreitWigner bw("bw", "", mass_axis, m0, width);
	RooDCBShape *cb;
	cb = new RooDCBShape("cb","cb", mass_axis, mean, sigma, alpha, n, alpha2, n2);
	RooGaussian gauss("gs", "gs", mass_axis, mean, sigma);
	RooFFTConvPdf signalRes("pdf", "pdf",mass_axis, bw, *cb);
	double iniSig = 0.2*invmass_num->Integral(1,100);
	int    lowBinNumber = invmass_num->FindBin(70.0);
	int    highBinNumber= invmass_num->FindBin(110.0);
	double iniBkg = 40*(invmass_num->GetBinContent(lowBinNumber)+invmass_num->GetBinContent(highBinNumber));
	RooRealVar nSig("nSig", "", iniSig, 0, invmass_num->GetEntries()*1.2);
	RooRealVar nBkg("nBkg", "", iniBkg, 0.5*iniBkg, invmass_num->GetEntries());
	RooAddPdf *model = new RooAddPdf("model", "", RooArgList(pdf_bg, signalRes),RooArgList(nBkg, nSig));   
 
	RooPlot* mass_Frame = mass_axis.frame(RooFit::Title("totalNum"),RooFit::Bins(40));
	model->fitTo(datahist_data);
	datahist_data.plotOn(mass_Frame);
	model->plotOn(mass_Frame, RooFit::Components(pdf_bg),
		RooFit::LineStyle(kDashed),
		RooFit::Normalization(1.0, RooAbsReal::RelativeExpected));
	model->plotOn(mass_Frame,
		RooFit::Components(RooArgSet(pdf_bg, signalRes)),
		RooFit::LineStyle(kDotted),
		RooFit::Normalization(1.0, RooAbsReal::RelativeExpected));
	model->plotOn(mass_Frame, 
		RooFit::LineStyle(kSolid));
	mass_Frame->Draw();
	c_fitMass->SaveAs("fit_totalNum_new.pdf");
 
	mass_axis.setRange("signal",70,110);
	RooAbsReal* igx_sig = signalRes.createIntegral(mass_axis,RooFit::NormSet(mass_axis),RooFit::Range("signal"));
	std::cout << "num = " <<  igx_sig->getVal()*(nSig.getVal()) << " error = " << igx_sig->getVal()*(nSig.getError()) <<  std::endl;
 	std::cout << "predict = " << invmass_prednum->Integral(lowBinNumber,highBinNumber) << std::endl; 
 	std::cout << "scale factor = " << igx_sig->getVal()*(nSig.getVal())/invmass_prednum->Integral(lowBinNumber,highBinNumber) << std::endl;
 
	TCanvas *cancompare = new TCanvas("cancompare","",600,600);
	cancompare->cd();
	RooPlot* compare_Frame = mass_axis.frame(RooFit::Title("compare"),RooFit::Bins(40));
	model->plotOn(compare_Frame,
		RooFit::Components(signalRes),
		RooFit::LineStyle(kDotted),
		RooFit::Normalization(1.0, RooAbsReal::RelativeExpected));
	invmass_prednum->Scale(igx_sig->getVal()*(nSig.getVal())*1.0/invmass_prednum->Integral(lowBinNumber,highBinNumber));
	RooDataHist datahist_prednum("both", "", mass_axis, invmass_prednum);
	datahist_prednum.plotOn(compare_Frame, RooFit::MarkerColor(kRed));
	compare_Frame->Draw();
	cancompare->SaveAs("compare_predvsnum_new.pdf");
 
//*************************************   Toy MC **************************************************************************//
 
	double totalnumError = sqrt(pow(0.13*igx_sig->getVal()*(nSig.getVal()), 2) + pow(igx_sig->getVal()*(nSig.getError()), 2));
	gRandom = new TRandom3(0);
	gRandom->SetSeed(0);
	double random_totalnum[NTOY]; 
	for(unsigned ir(1); ir<NTOY; ir++)	
		random_totalnum[ir] = igx_sig->getVal()*(nSig.getVal()) + totalnumError*(-1+ gRandom->Rndm()*2.0);
     
     	TMatrixDSym cov_pt = result_fitpt->GetCovarianceMatrix(); 
     	result_fitpt->Print("V");     
     	TVectorD mean_pt(3) ;
     	mean_pt(0) = result_fitpt->Parameter(0); 
     	mean_pt(1) = result_fitpt->Parameter(1);
     	mean_pt(2) = result_fitpt->Parameter(2);
     	RooRealVar central_slope("central_slope","",mean_pt(0)-result_fitpt->ParError(0), mean_pt(0)+result_fitpt->ParError(0));
     	RooRealVar central_const("central_const","",mean_pt(1)-result_fitpt->ParError(1),mean_pt(1)+result_fitpt->ParError(1));
     	RooRealVar central_index("central_index","",mean_pt(2)-result_fitpt->ParError(2),mean_pt(2)+result_fitpt->ParError(2));
     
     	RooMultiVarGaussian mvg_pt("mvg_pt","mvg_pt",RooArgList(central_slope,central_const,central_index),mean_pt,cov_pt);
     	RooDataSet* toymcdata_pt = mvg_pt.generate(RooArgSet(central_slope,central_const,central_index),NTOY);
     	std::ostringstream toymchistname;
     	TF1 *h_toymc_pt[NTOY];
     
     	TMatrixDSym cov_vtx = result_fitvtx->GetCovarianceMatrix(); 
     	result_fitvtx->Print("V");     
     	TVectorD mean_vtx(2) ;
     	mean_vtx(0) = result_fitvtx->Parameter(0); 
     	mean_vtx(1) = result_fitvtx->Parameter(1);
     	RooRealVar central_p0("central_p0","", mean_vtx(0)-result_fitvtx->ParError(0), mean_vtx(0)+result_fitvtx->ParError(0));
     	RooRealVar central_p1("central_p1","",mean_vtx(1)-result_fitvtx->ParError(1), mean_vtx(1)+result_fitvtx->ParError(1));
     
     	RooMultiVarGaussian mvg_vtx("mvg_vtx","mvg_vtx",RooArgList(central_p0,central_p1),mean_vtx,cov_vtx);
     	RooDataSet* toymcdata_vtx = mvg_vtx.generate(RooArgSet(central_p0,central_p1),NTOY);
     	TF1 *h_toymc_vtx[NTOY];
     
     	ofstream myfile;
     	myfile.open("ToyFakeRate_FullEcal.txt", std::ios_base::app | std::ios_base::out);
     	TH1D *p_scalefactor = new TH1D("p_scalefactor","fake rate scale factor; scale factor;",100,1000,2000);
     
     	float toyptvalue[graphPtBins][NTOY];
     	float lowtoyptvalue[graphPtBins];
     	float hightoyptvalue[graphPtBins];
     	for(unsigned ii(0); ii < graphPtBins; ii++){
     		lowtoyptvalue[ii] = 1;
     		hightoyptvalue[ii] = 0;
     	}
     	float toyvtxvalue[graphVtxBins][NTOY];
     	float lowtoyvtxvalue[graphVtxBins];
     	float hightoyvtxvalue[graphVtxBins];
     	for(unsigned ii(0); ii < graphVtxBins; ii++){
     		lowtoyvtxvalue[ii] = 1;
     		hightoyvtxvalue[ii] = 0;
     	}
     
     	for(int i(0); i<NTOY; i++){
     		double data1 = toymcdata_pt->get(i)->getRealValue("central_slope");
     		double data2 = toymcdata_pt->get(i)->getRealValue("central_const");
     		double data3 = toymcdata_pt->get(i)->getRealValue("central_index");
     		if(data1 > 1e6 || data2 > 1e6 || data3 > 1e6)continue;
     		toymchistname.str("");
     		toymchistname << "h_toymc_pt_" << i;
     		h_toymc_pt[i] = new TF1(toymchistname.str().c_str(), fakerate_ptDependence, 30,1000, 3);
     		h_toymc_pt[i]->SetParameter(0, data1);
     		h_toymc_pt[i]->SetParameter(1, data2);
     		h_toymc_pt[i]->SetParameter(2, data3);
     
     		double data4 = toymcdata_vtx->get(i)->getRealValue("central_p0");
     		double data5 = toymcdata_vtx->get(i)->getRealValue("central_p1");
     		toymchistname.str("");
     		toymchistname << "h_toymc_vtx_" << i;
     		h_toymc_vtx[i] = new TF1(toymchistname.str().c_str(), "pol1", 0,50);
     		h_toymc_vtx[i]->SetParameters(data4, data5);
     
     		invmass_prednum->Reset();
     		for(unsigned iEvt(0); iEvt < etreeEt.size(); iEvt++){
     			float w_ele = h_toymc_pt[i]->Eval(etreeEt[iEvt])*h_toymc_vtx[i]->Eval(etreeVtx[iEvt])*fit_fakerate_eta->Eval(etreeEta[iEvt]);
					if(w_ele > 0.5 || w_ele < 0)std::cout << "toy " << etreeEt[iEvt] << " " << etreeEta[iEvt] << " " << w_ele << std::endl;
     			invmass_prednum->Fill(etreeInvmass[iEvt], w_ele); 
     		}
     		invmass_prednum->Sumw2();
     
     		if(invmass_prednum->Integral(lowBinNumber,highBinNumber) > 0 && invmass_prednum->Integral(lowBinNumber,highBinNumber) < 1e20){
     			myfile << random_totalnum[i]/invmass_prednum->Integral(lowBinNumber,highBinNumber) << " " << data1 << " " << data2 << " " << data3 <<  " " << data4 << " " << data5 << std::endl;
     			p_scalefactor->Fill(random_totalnum[i]/invmass_prednum->Integral(lowBinNumber,highBinNumber));
     			for(unsigned ibin(1); ibin <= graphPtBins; ibin++){
     				double estimated = h_toymc_pt[i]->Eval(ptlowedge+0.5*ibin);
     				toyptvalue[ibin-1][i] = estimated;
     				if(estimated < lowtoyptvalue[ibin-1])lowtoyptvalue[ibin-1] = estimated;
     				if(estimated > hightoyptvalue[ibin-1])hightoyptvalue[ibin-1] = estimated;
     			}
     			for(unsigned ibin(1); ibin <= graphVtxBins; ibin++){
     				double estimated = h_toymc_vtx[i]->Eval(0.5*ibin);
     				toyvtxvalue[ibin-1][i] = estimated;
     				if(estimated < lowtoyvtxvalue[ibin-1])lowtoyvtxvalue[ibin-1] = estimated;
     				if(estimated > hightoyvtxvalue[ibin-1])hightoyvtxvalue[ibin-1] = estimated;
     			}
     		}
     	}
     
     // *************************   Calculated errors  ***************************************************************************************************************//
     	for(unsigned ibin(0); ibin < graphPtBins; ibin++){
     		TH1D *h_toyptdis = new TH1D("h_toyptdis","",50,lowtoyptvalue[ibin],hightoyptvalue[ibin]);
     		for(unsigned i(0); i < NTOY; i++)h_toyptdis->Fill(toyptvalue[ibin][i]);
     		h_toyptdis->Fit("gaus");
     		float fiterror = h_toyptdis->GetFunction("gaus")->GetParameter(2);
     		fr_pt_sigmaband->SetPoint(ibin, ptlowedge+0.5*ibin,fit_fakerate_pt->Eval(ptlowedge+0.5*ibin)); 
     		fr_pt_sigmaband->SetPointError(ibin, 0.25, fiterror);
     		fr_pt_ratioError->SetPoint(ibin, ptlowedge+0.5*ibin, 1);
     		fr_pt_ratioError->SetPointError(ibin, 0.25, fiterror/fit_fakerate_pt->Eval(ptlowedge+0.5*ibin));
     		delete h_toyptdis; 
     	}
     	for(unsigned ibin(0); ibin < graphVtxBins; ibin++){
     		TH1D *h_toyvtxdis = new TH1D("h_toyvtxdis","",50,lowtoyvtxvalue[ibin],hightoyvtxvalue[ibin]);
     		for(unsigned i(0); i < NTOY; i++)h_toyvtxdis->Fill(toyvtxvalue[ibin][i]);
     		h_toyvtxdis->Fit("gaus");
     		float fiterror = h_toyvtxdis->GetFunction("gaus")->GetParameter(2);
     		fr_vtx_sigmaband->SetPoint(ibin, vtxlowedge+0.5*ibin,fit_fakerate_vtx->Eval(vtxlowedge+0.5*ibin)); 
     		fr_vtx_sigmaband->SetPointError(ibin, 0.25, fiterror);
     		fr_vtx_ratioError->SetPoint(ibin, vtxlowedge+0.5*ibin, 1);
     		fr_vtx_ratioError->SetPointError(ibin, 0.25, fiterror/fit_fakerate_vtx->Eval(vtxlowedge+0.5*ibin));
     		delete h_toyvtxdis; 
     	}
     
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
     	fr_pt_ratio->SetMarkerStyle(20);
     	fr_pt_ratio->SetMarkerColor(kBlue);
     	fr_pt_ratio->SetLineColor(kBlue);
     	fr_pt_ratio->SetLineWidth(2);
     	fr_pt_ratio->SetFillColor(0);
     	fr_vtx_ratio->SetMarkerStyle(20);
     	fr_vtx_ratio->SetMarkerColor(kBlue);
     	fr_vtx_ratio->SetLineColor(kBlue);
     	fr_vtx_ratio->SetLineWidth(2);
     	fr_vtx_ratio->SetFillColor(0);
     
     	TLegend *leg = new TLegend(0.5,0.7,0.85,0.85);
     	leg->AddEntry(fr_bothcount_pt, "data");
     	leg->AddEntry(fit_fakerate_pt,"fit result");
     	leg->AddEntry(fr_pt_sigmaband,"total uncertainty");
     
     	TCanvas *canpt = new TCanvas("canpt","",600,600);
     	canpt->cd();
     	TPad *canpt_pad1 = new TPad("canpt_pad1", "pad1", 0, 0.3, 1, 1.0);
     	canpt_pad1->SetBottomMargin(0);
     	canpt_pad1->Draw();          
     	canpt_pad1->cd();          
     	TH1D *dummy_pt = new TH1D("",";Pt(GeV);fake rate",170,34,204);
     	dummy_pt->SetMaximum(0.06);
     	dummy_pt->Draw();
     	fr_bothcount_pt->Draw("EP same");
     	fit_fakerate_pt->Draw("same");
     	fr_pt_sigmaband->SetFillStyle(3005);
     	fr_pt_sigmaband->Draw("E2 same");
     	leg->Draw("same");
     
     	canpt->cd();   
     	TPad *canpt_pad2 = new TPad("canpt_pad2", "pad2", 0, 0.05, 1, 0.3);
     	canpt_pad2->SetTopMargin(0);
     	canpt_pad2->SetBottomMargin(0.3);
     	canpt_pad2->Draw();
     	canpt_pad2->cd(); 	
     	TH1D *dummy_ptratio = new TH1D("dummy_ptratio",";Pt(GeV);data/fit",17,30,200);
     	dummy_ptratio->SetMaximum(2);
     	dummy_ptratio->SetMinimum(0);
     	dummy_ptratio->Draw();
     	fr_pt_ratio->Draw("EP same");
     	fr_pt_ratioError->SetFillStyle(3005);
     	fr_pt_ratioError->Draw("E2 same");
     	canpt->SaveAs("elefake_pt_systematic_new.pdf");	
     	  
     	TCanvas *canvtx = new TCanvas("canvtx","",600,600);
     	canvtx->cd();
     	TPad *canvtx_pad1 = new TPad("canvtx_pad1", "pad1", 0, 0.3, 1, 1.0);
     	canvtx_pad1->SetBottomMargin(0);
     	canvtx_pad1->Draw();          
     	canvtx_pad1->cd();          
     	TH1D *dummy_vtx = new TH1D("",";nVtx;fake rate",44,0,44);
     	dummy_vtx->SetMaximum(0.1);
     	dummy_vtx->Draw();
     	fr_bothcount_vtx->Draw("EP same");
     	fit_fakerate_vtx->Draw("same");
     	fr_vtx_sigmaband->SetFillStyle(3005);
     	fr_vtx_sigmaband->Draw("E2 same");
     	leg->Draw("same");
     
     	canvtx->cd();   
     	TPad *canvtx_pad2 = new TPad("canvtx_pad2", "pad2", 0, 0.05, 1, 0.3);
     	canvtx_pad2->SetTopMargin(0);
     	canvtx_pad2->SetBottomMargin(0.3);
     	canvtx_pad2->Draw();
     	canvtx_pad2->cd(); 	
     	TH1D *dummy_vtxratio = new TH1D("dummy_vtxratio",";nVtx;data/fit",44,0,44);
     	dummy_vtxratio->SetMaximum(2);
     	dummy_vtxratio->SetMinimum(0);
     	dummy_vtxratio->Draw();
     	fr_vtx_ratio->Draw("EP same");
     	fr_vtx_ratioError->SetFillStyle(3005);
     	fr_vtx_ratioError->Draw("E2 same");
     	canvtx->SaveAs("elefake_vtx_systematic_new.pdf");	
     
     
     	TCanvas *caneta = new TCanvas("caneta","",600,600);
     	caneta->cd();
     	TPad *caneta_pad1 = new TPad("caneta_pad1", "pad1", 0, 0.05, 1, 1.0);
     	caneta_pad1->SetBottomMargin(0.1);
     	caneta_pad1->Draw();          
     	caneta_pad1->cd();          
     	TH1D *dummy_eta = new TH1D("",";eta;fake rate",50,0,2.5);
     	dummy_eta->SetMaximum(0.08);
     	dummy_eta->Draw();
     	fr_bothcount_eta->Draw("EPL same");
     	fr_eta_sigmaband->SetFillStyle(3005);
     	fr_eta_sigmaband->Draw("E2 same");
     	caneta->SaveAs("elefake_eta_systematic_new.pdf");
     
     	TCanvas *canscale = new TCanvas("scale","scale",600,600);
     	p_scalefactor->Draw("hist");
     	canscale->SaveAs("elefake_scalefactor_new.pdf");		
}


