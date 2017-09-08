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
#include "TMatrixDSym.h"
#include "../../../include/tdrstyle.C"
#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_tools.h"
#include "../../../include/analysis_mcData.h"

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
#define NTOY 1

bool isElectron(int PID, int momID){
   bool isEle;
   if(fabs(PID) == 11){ 
	   switch(momID){
	     case 1: isEle = true; break;
	     case 2: isEle = true; break;
	     case 3: isEle = true; break;
	     case 4: isEle = true; break;
	     case 5: isEle = true; break;
	     case 6: isEle = true; break;
	     case 21: isEle = true; break;
	     case 23: isEle = true; break;
	     default: isEle = false; break;
	   }
  }
  else isEle = false;

  return isEle;
}


double etaRates[]={0.0176057,0.0141562,0.0148432,0.0138896,0.0123009,0.0121965,0.0112304,0.0109691,0.0109968,0.0110605,0.0101498,0.010515,0.0121492,0.0112751,0.0100622,0.0107934,0.00985726,0.00986542,0.0114425,0.00914913,0.0098849,0.0107726,0.00980432,0.00995957,0.011031,0.0121443,0.0118273,0.0150705,0.0150705};

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
	double eta = x[0];
	for(int ieta(0); ieta < 29; ieta++){
		if(eta > ieta*0.05 && eta < (ieta+1)*0.05)return etaRates[ieta];
		else return 0.02;
	}
}
	

void fitFakeFunc(){//main 

	bool DoPlotPt(true);
	bool DoPlotPolEta(false);
	bool DoPlotEta(true);
	bool DoPlotVtx(true);
 
	gStyle->SetOptStat(0);
	int PtBins[] ={35,40,45,50,55,60,65,70,75,80,90,100,120,150,180};
	unsigned nPtBins = sizeof(PtBins)/sizeof(int);
	float EtaBins[30];
	for(int i(0); i<30; i++)EtaBins[i]=0.05*i;
	unsigned nEtaBins = sizeof(EtaBins)/sizeof(float);
	float VtxBins[]={0,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,48,50};
	unsigned nVtxBins = sizeof(VtxBins)/sizeof(float);
	TGraphErrors *fr_bothcount_pt = new TGraphErrors(nPtBins-1);
	TGraphErrors *fr_bothcount_eta= new TGraphErrors(nEtaBins-1);
	TGraphErrors *fr_bothcount_vtx= new TGraphErrors(nVtxBins-1);
	TGraph *fr_pt_upper = new TGraph(290);
	TGraph *fr_pt_lower = new TGraph(290);
	TGraph *fr_vtx_upper = new TGraph(92);
	TGraph *fr_vtx_lower = new TGraph(92);
	TGraphErrors *fr_true_eta= new TGraphErrors(nEtaBins-1);

	std::string bintype;
	std::string numtype;
	float lowcut;
	float signal1, error1;
	float fitmean, fitrms;
	float DYmean, DYrms;
	float Polmean, Polrms;
	float tmpsig, tmperror;
	float truevalue;
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

	double trueeta_den[nEtaBins-1];
	double trueeta_num[nEtaBins-1];

	std::ifstream Pt_file("MCEleFakeRate-2016ReReco-Bw-ker-pt-60-120.txt");
	std::ifstream Pt_DYfile("MCEleFakeRate-2016ReReco-Bw-ker-pt-60-120.txt");
	std::ifstream Pt_Polfile("MCEleFakeRate-2016ReReco-Bw-ker-pt-60-120.txt");
	std::ifstream Eta_file("MCEleFakeRate-2016ReReco-Bw-ker-eta-60-120.txt");
	std::ifstream Eta_DYfile("MCEleFakeRate-2016ReReco-Bw-ker-eta-60-120.txt");
	std::ifstream Eta_Polfile("MCEleFakeRate-2016ReReco-Bw-ker-eta-60-120.txt");
	std::ifstream Vtx_file("MCEleFakeRate-2016ReReco-Bw-ker-vtx-60-120.txt");
	std::ifstream Vtx_DYfile("MCEleFakeRate-2016ReReco-Bw-ker-vtx-60-120.txt");
	std::ifstream Vtx_Polfile("MCEleFakeRate-2016ReReco-Bw-ker-vtx-60-120.txt");

 
	if(DoPlotPt){
		if(Pt_file.is_open() && Pt_DYfile.is_open() && Pt_Polfile.is_open()){
			for(unsigned i(0); i<nPtBins-1; i++){ 
				Pt_file >> bintype >> numtype >> lowcut >> signal1 >> error1 >> fitmean >> fitrms >> truevalue;
				pt_den[i] = signal1;
				Pt_DYfile >> bintype >> numtype >> lowcut >> DYmean >> DYrms >> signal1 >> error1 >> truevalue; 
				Pt_Polfile >> bintype >> numtype >> lowcut >> Polmean >> Polrms >> signal1 >> error1 >> truevalue;
				double sysdiff = fabs(DYmean - fitmean) > fabs(truevalue - fitmean)? fabs(DYmean - fitmean):fabs(truevalue - fitmean);
				pt_denerror[i] = sqrt(fitrms*fitrms + sysdiff*sysdiff);
			}    
			for(unsigned i(0); i<nPtBins-1; i++){ 
				Pt_file >> bintype >> numtype >> lowcut >> fitmean >> fitrms >> signal1 >> error1 >> truevalue; 
				pt_num[i] = fitmean; 
				Pt_DYfile >> bintype >> numtype >> lowcut >> DYmean >> DYrms >> signal1 >> error1 >> truevalue; 
				Pt_Polfile >> bintype >> numtype >> lowcut >> Polmean >> Polrms >> signal1 >> error1 >> truevalue;
				double sysdiff = fabs(DYmean - fitmean) > fabs(truevalue - fitmean)? fabs(DYmean - fitmean):fabs(truevalue - fitmean);
				pt_numerror[i] = sqrt(fitrms*fitrms + sysdiff*sysdiff);
			}    
		}
		for(unsigned i(0); i<nPtBins-1; i++){
			double fakerate = pt_num[i]/pt_den[i];
			double error = sqrt(fakerate*fakerate*pt_denerror[i]*pt_denerror[i]/(pt_den[i]*pt_den[i])+ pt_numerror[i]*pt_numerror[i]/pt_den[i]/pt_den[i]);
			fr_bothcount_pt->SetPoint(i,(PtBins[i]+PtBins[i+1])/2.0, fakerate );
			fr_bothcount_pt->SetPointError(i, (PtBins[i+1]-PtBins[i])/2.0, error);
		}
	}

	if(DoPlotEta){
		if(Eta_file.is_open() && Eta_DYfile.is_open() && Eta_Polfile.is_open()){
			for(unsigned i(0); i<nEtaBins-1; i++){ 
				Eta_file >> bintype >> numtype >> lowcut >> fitmean >> fitrms >> signal1 >> error1 >> truevalue;
				eta_den[i] = fitmean;
				trueeta_den[i] = truevalue;
				Eta_DYfile >> bintype >> numtype >> lowcut >> DYmean >> DYrms >> signal1 >> error1 >> truevalue;
				Eta_Polfile >> bintype >> numtype >> lowcut >> Polmean >> Polrms >> signal1 >> error1 >> truevalue;
				double sysdiff = fabs(DYmean - fitmean) > fabs(truevalue - fitmean)? fabs(DYmean - fitmean):fabs(truevalue - fitmean);
				eta_denerror[i] = sqrt(fitrms*fitrms + sysdiff*sysdiff);
			}    
			for(unsigned i(0); i<nEtaBins-1; i++){ 
				Eta_file >> bintype >> numtype >> lowcut >> fitmean >> fitrms >> signal1 >> error1 >> truevalue;
				eta_num[i] = fitmean;
				trueeta_num[i] = truevalue;
				Eta_DYfile >> bintype >> numtype >> lowcut >>  DYmean >> DYrms >> signal1 >> error1 >> truevalue;
				Eta_Polfile >> bintype >> numtype >> lowcut >> Polmean >> Polrms >> signal1 >> error1 >> truevalue;
				double sysdiff = fabs(DYmean - fitmean) > fabs(truevalue - fitmean)? fabs(DYmean - fitmean):fabs(truevalue - fitmean);
				eta_numerror[i] = sqrt(fitrms*fitrms + sysdiff*sysdiff);
			}    
		}
		for(unsigned i(0); i<nEtaBins-1; i++){
			double fakerate = eta_num[i]/eta_den[i];
			double error = sqrt(fakerate*fakerate*eta_denerror[i]*eta_denerror[i]/(eta_den[i]*eta_den[i])+ eta_numerror[i]*eta_numerror[i]/eta_den[i]/eta_den[i]);
			std::cout << "eta " << EtaBins[i] << " " << fakerate << std::endl;
			fr_bothcount_eta->SetPoint(i,(EtaBins[i]+EtaBins[i+1])/2.0, fakerate );
			fr_bothcount_eta->SetPointError(i, (EtaBins[i+1]-EtaBins[i])/2.0, error);
			fr_true_eta->SetPoint(i,(EtaBins[i]+EtaBins[i+1])/2.0, trueeta_num[i]/trueeta_den[i]);
		}
	}

	if(DoPlotVtx){
		if(Vtx_file.is_open() && Vtx_DYfile.is_open() && Vtx_Polfile.is_open()){
			for(unsigned i(0); i<nVtxBins-1; i++){ 
				Vtx_file >> bintype >> numtype >> lowcut >> fitmean >> fitrms >> signal1 >> error1 >> truevalue;
				vtx_den[i] = fitmean;
				Vtx_DYfile >> bintype >> numtype >> lowcut >>  DYmean >> DYrms >> signal1 >> error1 >> truevalue;
				Vtx_Polfile >> bintype >> numtype >> lowcut >> Polmean >> Polrms >> signal1 >> error1 >> truevalue;
				double sysdiff = fabs(DYmean - fitmean) > fabs(truevalue - fitmean)? fabs(DYmean - fitmean):fabs(truevalue - fitmean);
				vtx_denerror[i] = sqrt(fitrms*fitrms + sysdiff*sysdiff);
			}    
			for(unsigned i(0); i<nVtxBins-1; i++){ 
				Vtx_file >> bintype >> numtype >> lowcut >> fitmean >> fitrms >> signal1 >> error1 >> truevalue;
				vtx_num[i] = fitmean;
				Vtx_DYfile >> bintype >> numtype >> lowcut >> DYmean >> DYrms >> signal1 >> error1 >> truevalue;
				Vtx_Polfile >> bintype >> numtype >> lowcut >> Polmean >> Polrms >> signal1 >> error1 >> truevalue;
				double sysdiff = fabs(DYmean - fitmean) > fabs(truevalue - fitmean)? fabs(DYmean - fitmean):fabs(truevalue - fitmean);
				vtx_numerror[i] = sqrt(fitrms*fitrms + sysdiff*sysdiff);
			}    
		}
		for(unsigned i(0); i<nVtxBins-1; i++){
			double fakerate = vtx_num[i]/vtx_den[i];
			double error = sqrt(fakerate*fakerate*vtx_denerror[i]*vtx_denerror[i]/(vtx_den[i]*vtx_den[i])+ vtx_numerror[i]*vtx_numerror[i]/vtx_den[i]/vtx_den[i]);
			fr_bothcount_vtx->SetPoint(i,(VtxBins[i]+VtxBins[i+1])/2.0, fakerate );
			fr_bothcount_vtx->SetPointError(i, (VtxBins[i+1]-VtxBins[i])/2.0, error);
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

	TF1 *fakerate_pt;
	TF1 *fakerate_vtx = new TF1("vtxfake","pol1",0,100);
  TF1 *fakerate_eta = new TF1("fakerate_eta", fakerate_etaDependence, 0, 1.5, 0); 
	TFitResultPtr result_fitpt;
	TFitResultPtr result_fitvtx;

  TH1F *dummy_pt = new TH1F("",";Pt(GeV);fake rate",140,34,174);
  TCanvas *canvas_pt = new TCanvas("fake-rate vs. p_{T}","",600,600);
	TCanvas *canvas_vtx = new TCanvas("fake-rate vs. nVtx"," fake rate",600,600);
  canvas_pt->cd();
  dummy_pt->SetMaximum(0.06);
  dummy_pt->Draw();
	TF1 *f1 = new TF1("f1", fakerate_ptDependence,30,1000,3);

	TVirtualFitter::SetMaxIterations(1000000);
	f1->SetParameter(0, 50);
	f1->SetParLimits(0, 0, 10000);
	f1->SetParameter(1, 88);
	f1->SetParameter(2, -0.5);
	//f1->SetParameter(3, 0.00252886);
	f1->SetParNames("slope","constant","index");
	fr_bothcount_pt->Draw("P same");
	TLegend *leg=new TLegend(0.6,0.7,0.9,0.9);
	result_fitpt = fr_bothcount_pt->Fit("f1","R S");
	fakerate_pt = fr_bothcount_pt->GetFunction("f1");
  canvas_pt->SaveAs("elefakepho_rate_pt.pdf");

	TH1F *dummy_eta = new TH1F("",";#eta;fake rate",14,0,1.4);
	TCanvas *canvas_eta = new TCanvas("fake-rate vs. #eta","fake rate",600,600);
	canvas_eta->cd();
	dummy_eta->SetMaximum(0.08);
	dummy_eta->Draw();
	fr_bothcount_eta->Draw("PL same");
	fr_true_eta->SetLineColor(kRed);
	fr_true_eta->Draw("same");
	canvas_eta->SaveAs("elefakepho_rate_eta.pdf");

	canvas_vtx->cd();
	TH1F *dummy_vtx = new TH1F("",";nVtx;fake rate",44,0,44);
	dummy_vtx->SetMaximum(0.08);
	dummy_vtx->Draw();
	fr_bothcount_vtx->SetTitle("CMS Preliminary");
	fr_bothcount_vtx->GetXaxis()->SetTitle("nVertex");
	fr_bothcount_vtx->GetYaxis()->SetTitle("Fake rate");
	fr_bothcount_vtx->Draw("P same");
	fakerate_vtx->SetParameter(0,0.02);
	fakerate_vtx->SetParameter(1,1);
	result_fitvtx = fr_bothcount_vtx->Fit("vtxfake","S");
	fakerate_vtx->Draw("same");
	canvas_vtx->SaveAs("elefakepho_rate_vtx.pdf");

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libRooFitClasses.so");

  TChain *etree = new TChain("FakeRateTree");
  etree->Add("/uscms_data/d3/mengleis/plot_elefakepho_DYToLL.root");

  float invmass=0; 
  float tagPt=0; 
  float probePt=0; 
  float probeEta=0;
  float probePhi=0;
  bool  vetovalue=0;
	bool  FSRveto=0;
  int   nVertex=0;
	std::vector<int>   *mcPID_tnp=0;
	std::vector<float> *mcEta_tnp=0;
	std::vector<float> *mcPhi_tnp=0;
	std::vector<float> *mcPt_tnp=0;
	std::vector<int> *mcMomPID_tnp=0;
	std::vector<int> *mcGMomPID_tnp=0;
  etree->SetBranchAddress("invmass",   &invmass); 
  etree->SetBranchAddress("tagPt",     &tagPt);
  etree->SetBranchAddress("probePt",   &probePt);
  etree->SetBranchAddress("probeEta",  &probeEta);
  etree->SetBranchAddress("probePhi",  &probePhi);
  etree->SetBranchAddress("vetovalue", &vetovalue);
  etree->SetBranchAddress("FSRveto",   &FSRveto);
  etree->SetBranchAddress("nVertex",   &nVertex);
	etree->SetBranchAddress("mcPID",			&mcPID_tnp);
	etree->SetBranchAddress("mcEta",			&mcEta_tnp);
	etree->SetBranchAddress("mcPhi",			&mcPhi_tnp);
	etree->SetBranchAddress("mcPt",			&mcPt_tnp);
	etree->SetBranchAddress("mcMomPID",	&mcMomPID_tnp);
	etree->SetBranchAddress("mcGMomPID",	&mcGMomPID_tnp);

  TH1F* invmass_den = new TH1F("invmass_den", "invmass_den",60,60,120);
  TH1F* invmass_num = new TH1F("invmass_num", "invmass_num",60,60,120);
	TH1F* invmass_true= new TH1F("invmass_true","invmass_true",60,60,120);
	TH1F* h_invmass_bg  = new TH1F("invmass_bg",  "invmass_bg", 60,60,120);
  TH1F* invmass_prednum = new TH1F("invmass_prednum", "invmass_prednum",60,60,120);
  
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
	  if(vetovalue== false)invmass_den->Fill(invmass);
	  else if(vetovalue==true && FSRveto == true)invmass_num->Fill(invmass);

		bool isZee(false);
		double mindRprobe(0.5);
		unsigned probeIndex(0);
		for(unsigned iMC(0); iMC<mcPID_tnp->size(); iMC++){
			double dR2 = DeltaR((*mcEta_tnp)[iMC], (*mcPhi_tnp)[iMC], probeEta,probePhi);
			double dE2 = fabs((*mcPt_tnp)[iMC] - probePt)/probePt;
			if(dR2 < mindRprobe && dE2 < 0.1){mindRprobe=dR2; probeIndex=iMC;}
		}
		if(mindRprobe < 0.1){
			isZee = isElectron(fabs((*mcPID_tnp)[probeIndex]), fabs((*mcMomPID_tnp)[probeIndex]));
		}
		if(isZee && vetovalue== true && FSRveto == true)invmass_true->Fill(invmass);


		float w_ele = fakerate_pt->Eval(probePt)*fakerate_vtx->Eval(nVertex)*fakerate_eta->Eval(fabs(probeEta));
		if(vetovalue==false){
			invmass_prednum->Fill(invmass, w_ele);
			etreeEt.push_back(probePt);
			etreeVtx.push_back(nVertex);
			etreeEta.push_back(fabs(probeEta));
			etreeInvmass.push_back(invmass);
		} 
  }
	invmass_prednum->Sumw2();
	TCanvas *canpred = new TCanvas("canpred","",600,600);
	canpred->cd();
	invmass_prednum->Draw();

  TChain *bgtree = new TChain("BGTree");
  //bgtree->Add("/uscms_data/d3/mengleis/plot_bgtemplate_DYJetsToLL.root");
  bgtree->Add("/uscms_data/d3/mengleis/ReMiniAOD/plot_bgtemplate-SingleMu_ReMiniAOD.root");
  float invmass_bg=0; 
  bool  vetovalue_bg=0;
  bgtree->SetBranchAddress("invmass",   &invmass_bg); 
  bgtree->SetBranchAddress("vetovalue", &vetovalue_bg);

  for(unsigned iEvt(0); iEvt < bgtree->GetEntries(); iEvt++){
		bgtree->GetEntry(iEvt);
	  if(vetovalue_bg== true){h_invmass_bg->Fill(invmass_bg);}
  }

	TH1F  *h_DYinvmass = new TH1F("h_DYinvmass","h_DYinvmass",60,60,120);
	TChain *DYtree = new TChain("FakeRateTree");
	DYtree->Add("/uscms_data/d3/mengleis/plot_elefakepho_DYToLL.root");

	float DY_invmass=0; 
	float DY_tagPt=0; 
	float DY_tagEta=0; 
	float DY_tagPhi=0; 
	float DY_probePt=0; 
	float DY_probeEta=0;
	float DY_probePhi=0;
	bool  DY_vetovalue=0;
	int   DY_nVertex=0;
	std::vector<int>   *mcPID=0;
	std::vector<float> *mcEta=0;
	std::vector<float> *mcPhi=0;
	std::vector<float> *mcPt=0;
	std::vector<int> *mcMomPID=0;
	std::vector<int> *mcGMomPID=0;
	DYtree->SetBranchAddress("invmassUncalib",   &DY_invmass); 
	DYtree->SetBranchAddress("tagPt",     &DY_tagPt);
	DYtree->SetBranchAddress("tagEta",    &DY_tagEta);
	DYtree->SetBranchAddress("tagPhi",    &DY_tagPhi);
	DYtree->SetBranchAddress("probePt",   &DY_probePt);
	DYtree->SetBranchAddress("probeEta",  &DY_probeEta);
	DYtree->SetBranchAddress("probePhi",  &DY_probePhi);
	DYtree->SetBranchAddress("vetovalue", &DY_vetovalue);
	DYtree->SetBranchAddress("nVertex",   &DY_nVertex);
	DYtree->SetBranchAddress("mcPID",			&mcPID);
	DYtree->SetBranchAddress("mcEta",			&mcEta);
	DYtree->SetBranchAddress("mcPhi",			&mcPhi);
	DYtree->SetBranchAddress("mcPt",			&mcPt);
	DYtree->SetBranchAddress("mcMomPID",	&mcMomPID);
	DYtree->SetBranchAddress("mcGMomPID",	&mcGMomPID);
	for(unsigned iEvt(0); iEvt < DYtree->GetEntries(); iEvt++){
		DYtree->GetEntry(iEvt);
		if(DY_probePt < 35)continue;
		bool isZee(false);
		double mindRtag(0.3), mindRprobe(0.3);
		unsigned tagIndex(0), probeIndex(0);
		for(unsigned iMC(0); iMC<mcPID->size(); iMC++){
			double dR1 = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], DY_tagEta, DY_tagPhi);
			double dR2 = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], DY_probeEta,DY_probePhi);
			double dE1 = fabs((*mcPt)[iMC] - DY_tagPt)/DY_tagPt;
			double dE2 = fabs((*mcPt)[iMC] - DY_probePt)/DY_probePt;
			if(dR1 < mindRtag && dE1 < 0.1){mindRtag=dR1; tagIndex=iMC;}
			if(dR2 < mindRprobe && dE2 < 0.1){mindRprobe=dR2; probeIndex=iMC;}
		}
		if(mindRtag < 0.1 && mindRprobe < 0.1){
			bool isZe(false),isZg(false);
			isZe = isElectron(fabs((*mcPID)[tagIndex]), fabs((*mcMomPID)[tagIndex]));
			isZg = isElectron(fabs((*mcPID)[probeIndex]), fabs((*mcMomPID)[probeIndex]));
			if(isZe && isZg)isZee=true; 
		}
		if(isZee)h_DYinvmass->Fill(DY_invmass);
	}
	h_DYinvmass->Sumw2();

  RooRealVar mass_axis("invmass","invmass",70,110);
  TCanvas *c_fitMass = new TCanvas("c_fitMass", "", 600, 600);
  c_fitMass->cd();

  RooDataHist datahist_data("both", "", mass_axis, invmass_num);
  RooDataHist datahist_true("both", "", mass_axis, invmass_true);
	RooDataHist datahist_bg("bg","",mass_axis, h_invmass_bg);
	RooHistPdf  pdf_bg("pdf_bg","pdf_bg",mass_axis, datahist_bg);
	RooDataHist datahist_DY("DYDataSet","DYDataSet", mass_axis, h_DYinvmass);
	RooHistPdf  DYpdf("DYpdf","DYpdf", mass_axis, datahist_DY);
	DYpdf.setInterpolationOrder(1);

  RooRealVar m0( "m0", "m0", 91.188, 80,100);
  RooRealVar width( "width", "width", 2.495, 0, 15);
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
  //RooFFTConvPdf signalRes("pdf", "pdf",mass_axis, bw, *cb);
	RooFFTConvPdf signalRes("pdf", "pdf",mass_axis, DYpdf, gauss);

  double iniSig = 0.8*invmass_num->Integral(1,60);
  double iniBkg = 100*(invmass_num->GetBinContent(1)+invmass_num->GetBinContent(60));
  RooRealVar nSig("nSig", "", iniSig, 0, invmass_num->GetEntries()*1.2);
  RooRealVar nBkg("nBkg", "", iniBkg, 0, invmass_num->GetEntries());
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
	datahist_true.plotOn(mass_Frame, RooFit::LineColor(kGreen));
  mass_Frame->Draw();

  c_fitMass->SaveAs("fit_totalNum.pdf");

  mass_axis.setRange("signal",70,110);
  RooAbsReal* igx_sig = signalRes.createIntegral(mass_axis,RooFit::NormSet(mass_axis),RooFit::Range("signal"));
  std::cout << "num = " <<  igx_sig->getVal()*(nSig.getVal()) << " " << igx_sig->getVal()*(nSig.getError()) <<  " true= "<< invmass_true->Integral(10,50) << std::endl;
	std::cout << "predict = " << invmass_prednum->Integral(10,50); 
	std::cout << "scale factor = " << igx_sig->getVal()*(nSig.getVal())/invmass_prednum->Integral(10,50) << std::endl;

	TCanvas *cancompare = new TCanvas("cancompare","",600,600);
  RooPlot* compare_Frame = mass_axis.frame(RooFit::Title("compare"),RooFit::Bins(40));
  model->plotOn(compare_Frame,
			 RooFit::Components(signalRes),
			 RooFit::LineStyle(kDotted),
			 RooFit::Normalization(1.0, RooAbsReal::RelativeExpected));
	invmass_prednum->Scale(igx_sig->getVal()*(nSig.getVal())*1.0/invmass_prednum->Integral(11,50));
  RooDataHist datahist_prednum("both", "", mass_axis, invmass_prednum);
	datahist_prednum.plotOn(compare_Frame, RooFit::MarkerColor(kRed));
	datahist_true.plotOn(compare_Frame, RooFit::MarkerColor(kGreen));
	compare_Frame->Draw();
	cancompare->SaveAs("compare_predvsnum.pdf");
	
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
	myfile.open("ToyFakeRate.txt", std::ios_base::app | std::ios_base::out);
	TH1F *p_scalefactor = new TH1F("p_scalefactor","p_scalefactor",5000,0,5000);

	TProfile *pro_pt = new TProfile("pro_pt","",290,35,180);
	pro_pt->SetErrorOption("s");
	TProfile *pro_vtx = new TProfile("pro_vtx","",92,0,46);
	pro_vtx->SetErrorOption("s");
	for(int i(0); i<NTOY; i++){
		double data1 = toymcdata_pt->get(i)->getRealValue("central_slope");
		double data2 = toymcdata_pt->get(i)->getRealValue("central_const");
		double data3 = toymcdata_pt->get(i)->getRealValue("central_index");
		if(data1 > 1e6 || data2 > 1e6 || data3 > 1e6)continue;
		toymchistname.str("");
		toymchistname << "h_toymc_pt_" << i;
		//h_toymc_pt[i] = new TF1(toymchistname.str().c_str(), fakerate_ptDependence, 35,10000, 4);
		h_toymc_pt[i] = new TF1(toymchistname.str().c_str(), fakerate_ptDependence, 30,1000, 3);
		//h_toymc_pt[i]->SetParameters(data1, data2, data3, data4);
		h_toymc_pt[i]->SetParameter(0, data1);
		h_toymc_pt[i]->SetParameter(1, data2);
		h_toymc_pt[i]->SetParameter(2, data3);
		//Int_t status = fr_bothcount_pt->Fit(toymchistname.str().c_str(),"R");
		//if(h_toymc_pt[i]->GetParError(0)/h_toymc_pt[i]->GetParameter(0) <= 1.0 && h_toymc_pt[i]->GetParError(1)/h_toymc_pt[i]->GetParameter(1) <= 1.0 && h_toymc_pt[i]->GetParError(2)/h_toymc_pt[i]->GetParameter(2) <= 1) myfile << "good point " << data1 << " " << data2 << " " << data3 << " status = " << status  << std::endl;

		double data5 = toymcdata_vtx->get(i)->getRealValue("central_p0");
		double data6 = toymcdata_vtx->get(i)->getRealValue("central_p1");
		toymchistname.str("");
		toymchistname << "h_toymc_vtx_" << i;
		h_toymc_vtx[i] = new TF1(toymchistname.str().c_str(), "pol1", 0,50);
		h_toymc_vtx[i]->SetParameters(data5, data6);
	//	h_toymc_vtx[i]->SetLineColorAlpha(kBlue, 0.35);
	//	h_toymc_vtx[i]->Draw("same");


		invmass_prednum->Reset();
		for(unsigned iEvt(0); iEvt < etreeEt.size(); iEvt++){
			float w_ele = h_toymc_pt[i]->Eval(etreeEt[iEvt])*h_toymc_vtx[i]->Eval(etreeVtx[iEvt])*fakerate_eta->Eval(etreeEta[iEvt]);
			invmass_prednum->Fill(etreeInvmass[iEvt], w_ele); 
		}
		invmass_prednum->Sumw2();
		//std::cout << "predict = " << invmass_prednum->Integral(11,50); 
		//std::cout << "scale factor = " << igx_sig->getVal()*(nSig.getVal())/invmass_prednum->Integral(11,50) << std::endl;
		if(invmass_prednum->Integral(11,50) > 0 && invmass_prednum->Integral(11,50) < 1e20){
			myfile << igx_sig->getVal()*(nSig.getVal())/invmass_prednum->Integral(11,50) << " " << data1 << " " << data2 << " " << data3 <<  " " << data5 << " " << data6 << std::endl;
			p_scalefactor->Fill(igx_sig->getVal()*(nSig.getVal())/invmass_prednum->Integral(11,50));
			for(unsigned ibin(1); ibin <= 290; ibin++){
				double estimated = h_toymc_pt[i]->Eval(pro_pt->GetBinCenter(ibin));
				pro_pt->Fill( 35+(180.0-35.0)/290.0*ibin, estimated);
			}
			for(unsigned ibin(1); ibin <= 92; ibin++){
				double estimated = h_toymc_vtx[i]->Eval(pro_vtx->GetBinCenter(ibin));
				pro_vtx->Fill( 46.0/92.0*ibin, estimated);
			}
		}
	}
	for(int i(0); i<290; i++)fr_pt_upper->SetPoint(i, (180.0-35.0)/290*i + 35.0, pro_pt->GetBinContent(i+1)+ pro_pt->GetBinError(i+1));
	for(int i(0); i<290; i++)fr_pt_lower->SetPoint(i, (180.0-35.0)/290*i + 35.0, pro_pt->GetBinContent(i+1)- pro_pt->GetBinError(i+1));
	for(int i(0); i<92; i++)fr_vtx_upper->SetPoint(i, 46.0/92.0*i, pro_vtx->GetBinContent(i+1)+ pro_vtx->GetBinError(i+1));
	for(int i(0); i<92; i++)fr_vtx_lower->SetPoint(i, 46.0/92.0*i, pro_vtx->GetBinContent(i+1)- pro_vtx->GetBinError(i+1));

	TCanvas *ptsyscan = new TCanvas("ptsyscan","",600,600);
	ptsyscan->cd();
	dummy_pt->Draw();
	f1->SetParameter(0, 50);
	f1->SetParLimits(0, 0, 10000);
	f1->SetParameter(1, 88);
	f1->SetParameter(2, -0.5);
	f1->SetParNames("slope","constant","index");
	fr_bothcount_pt->SetTitle("CMS Preliminary");
	fr_bothcount_pt->GetXaxis()->SetTitle("p_{T} (GeV)");
	fr_bothcount_pt->GetYaxis()->SetTitle("Fake rate");
	fr_bothcount_pt->Draw("P same");
	fr_bothcount_pt->Fit("f1","R S");
	fr_pt_upper->SetLineColor(kBlack);
	fr_pt_upper->SetLineStyle(2);
	fr_pt_upper->Draw("L same");
	fr_pt_lower->SetLineColor(kBlack);
	fr_pt_lower->SetLineStyle(2);
	fr_pt_lower->Draw("L same");
	for(unsigned it(0); it < NTOY; it++)h_toymc_pt[it]->Draw("same");
	ptsyscan->SaveAs("elefake_pt_systematic.pdf");	
	  
	TCanvas *vtxsyscan = new TCanvas("vtxsyscan","",600,600);
	vtxsyscan->cd();
	canvas_vtx->cd();
	dummy_vtx->SetMaximum(0.08);
	dummy_vtx->Draw();
	fr_bothcount_vtx->SetTitle("CMS Preliminary");
	fr_bothcount_vtx->GetXaxis()->SetTitle("nVertex");
	fr_bothcount_vtx->GetYaxis()->SetTitle("Fake rate");
	fr_bothcount_vtx->Draw("P same");
	fr_bothcount_vtx->Fit("pol1","R S");
	fr_vtx_upper->SetLineColor(kBlack);
	fr_vtx_upper->SetLineStyle(2);
	fr_vtx_upper->Draw("L same");
	fr_vtx_lower->SetLineColor(kBlack);
	fr_vtx_lower->SetLineStyle(2);
	fr_vtx_lower->Draw("L same");
	vtxsyscan->SaveAs("elefake_vtx_systematic.pdf");	
}


