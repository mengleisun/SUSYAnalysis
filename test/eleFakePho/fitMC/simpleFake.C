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
#include "RooMultiVarGaussian.h"
#include "../../../include/RooCMSShape.h"
#include "../../../include/RooDCBShape.h"
#include "../../../include/RooUserPoly.h"
#include "../../../include/analysis_fakes.h"
#define NTOY 100

Double_t fakerate_ptDependence(Double_t *x, Double_t *par)
{
  double slope = par[0];
  double constant = par[1];
  double index = par[2];
  double coeff = par[3]; 

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
	

void simpleFake(){//main 

	gStyle->SetOptStat(0);
	int PtBins[] ={35,40,45,50,55,60,65,70,75,80,90,100,120,150,180};
	unsigned nPtBins = sizeof(PtBins)/sizeof(int);
	TGraphErrors *fr_bothcount_pt = new TGraphErrors(nPtBins-1);
	TGraphErrors *fr_bothcount_true = new TGraphErrors(nPtBins-1);
	TGraphErrors *fr_bothcount_ratio = new TGraphErrors(nPtBins-1);

	std::string bintype;
	std::string numtype;
	float lowcut;
	float signal1, error1;
	float truemean;
	double pt_den[nPtBins-1];
	double pt_denerror[nPtBins-1];
	double pt_num[nPtBins-1];
	double pt_numerror[nPtBins-1];
	double true_den[nPtBins-1];
	double true_num[nPtBins-1];

	std::ifstream Pt_file("MCEleFakeRate-2016ReReco-Bw-ker-pt-60-120.txt");
 
		if(Pt_file.is_open()){ 
			for(unsigned i(0); i<nPtBins-1; i++){ 
				Pt_file >> bintype >> numtype >> lowcut >> signal1 >> error1 >> truemean; 
				pt_den[i] = signal1;
				pt_denerror[i] = error1;
				true_den[i] = truemean;
			}    
			for(unsigned i(0); i<nPtBins-1; i++){ 
				Pt_file >> bintype >> numtype >> lowcut >> signal1 >> error1 >> truemean; 
				pt_num[i] = signal1;
				pt_numerror[i] = error1;
				true_num[i] = truemean; 
			}    
		}
		for(unsigned i(0); i<nPtBins-1; i++){
			double fakerate = pt_num[i]/pt_den[i];
			double error = sqrt(fakerate*fakerate*pt_denerror[i]*pt_denerror[i]/(pt_den[i]*pt_den[i])+ pt_numerror[i]*pt_numerror[i]/pt_den[i]/pt_den[i]);
			fr_bothcount_pt->SetPoint(i,(PtBins[i]+PtBins[i+1])/2.0, fakerate );
			fr_bothcount_pt->SetPointError(i, (PtBins[i+1]-PtBins[i])/2.0, error);
			fr_bothcount_true->SetPoint(i, (PtBins[i]+PtBins[i+1])/2.0, true_num[i]/true_den[i]);
			fr_bothcount_ratio->SetPoint(i, (PtBins[i]+PtBins[i+1])/2.0, fakerate/(true_num[i]/true_den[i]));
		}
	
	TH1F *dummy = new TH1F("dummy","",20,0,200);
	dummy->SetMaximum(2);
	dummy->Draw();
	fr_bothcount_pt->Draw("P same");
	fr_bothcount_true->SetLineColor(kRed);
	fr_bothcount_true->Draw("same");
	fr_bothcount_ratio->Draw("same");

}


