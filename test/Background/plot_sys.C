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
void plot_sys(){//main  

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");
	
	TFile *file_VG  = TFile::Open("bkgTree_eg_VGBkg_ReMiniAOD.root");

	TH1F *p_PhoEt = (TH1F*)file_VG->Get("p_PhoEt");
	TH1F *p_LepPt = (TH1F*)file_VG->Get("p_LepPt");
	TH1F *p_MET   = (TH1F*)file_VG->Get("p_MET");
	TH1F *p_Mt    = (TH1F*)file_VG->Get("p_Mt");
	TH1F *p_dPhi  = (TH1F*)file_VG->Get("p_dPhiEleMET");
	TH1F *p_HT    = (TH1F*)file_VG->Get("p_HT");

	TH1F *jesup_MET = (TH1F*)file_VG->Get("jesup_MET");
	TH1F *jesup_Mt = (TH1F*)file_VG->Get("jesup_Mt");
	TH1F *jesup_HT = (TH1F*)file_VG->Get("jesup_HT");
	TH1F *jesup_dPhiEleMET = (TH1F*)file_VG->Get("jesup_dPhiEleMET");

	TH1F *jesdo_MET = (TH1F*)file_VG->Get("jesdo_MET");
	TH1F *jesdo_Mt = (TH1F*)file_VG->Get("jesdo_Mt");
	TH1F *jesdo_HT = (TH1F*)file_VG->Get("jesdo_HT");
	TH1F *jesdo_dPhiEleMET = (TH1F*)file_VG->Get("jesdo_dPhiEleMET");

	TH1F *jerup_MET = (TH1F*)file_VG->Get("jerup_MET");
	TH1F *jerup_Mt = (TH1F*)file_VG->Get("jerup_Mt");
	TH1F *jerup_dPhiEleMET = (TH1F*)file_VG->Get("jerup_dPhiEleMET");

	TH1F *jerdo_MET = (TH1F*)file_VG->Get("jerdo_MET");
	TH1F *jerdo_Mt = (TH1F*)file_VG->Get("jerdo_Mt");
	TH1F *jerdo_dPhiEleMET = (TH1F*)file_VG->Get("jerdo_dPhiEleMET");

	TH1F *scaleup_PhoEt = (TH1F*)file_VG->Get("scaleup_PhoEt");
	TH1F *scaleup_PhoEta = (TH1F*)file_VG->Get("scaleup_PhoEta");
	TH1F *scaleup_LepPt = (TH1F*)file_VG->Get("scaleup_LepPt");
	TH1F *scaleup_LepEta = (TH1F*)file_VG->Get("scaleup_LepEta");
	TH1F *scaleup_MET = (TH1F*)file_VG->Get("scaleup_MET");
	TH1F *scaleup_Mt = (TH1F*)file_VG->Get("scaleup_Mt");
	TH1F *scaleup_HT = (TH1F*)file_VG->Get("scaleup_HT");
	TH1F *scaleup_dPhiEleMET = (TH1F*)file_VG->Get("scaleup_dPhiEleMET");

	TH1F *scaledo_PhoEt = (TH1F*)file_VG->Get("scaledo_PhoEt");
	TH1F *scaledo_PhoEta = (TH1F*)file_VG->Get("scaledo_PhoEta");
	TH1F *scaledo_LepPt = (TH1F*)file_VG->Get("scaledo_LepPt");
	TH1F *scaledo_LepEta = (TH1F*)file_VG->Get("scaledo_LepEta");
	TH1F *scaledo_MET = (TH1F*)file_VG->Get("scaledo_MET");
	TH1F *scaledo_Mt = (TH1F*)file_VG->Get("scaledo_Mt");
	TH1F *scaledo_HT = (TH1F*)file_VG->Get("scaledo_HT");
	TH1F *scaledo_dPhiEleMET = (TH1F*)file_VG->Get("scaledo_dPhiEleMET");

	//for(int ibin(1); ibin < p_PhoEt->GetSize(); ibin++){
	//	double syserror(0);
	//	syserror += pow((scaleup_PhoEt->GetBinContent(ibin)-p_PhoEt->GetBinContent(ibin)),2);
	//	syserror += pow((scaledo_PhoEt->GetBinContent(ibin)-p_PhoEt->GetBinContent(ibin)),2);
	//	std::cout << "Et " << sqrt(syserror)/p_PhoEt->GetBinContent(ibin) << std::endl;
	//}	
	//for(int ibin(1); ibin < p_LepPt->GetSize(); ibin++){
	//	double syserror(0);
	//	syserror += pow((scaleup_LepPt->GetBinContent(ibin)-p_LepPt->GetBinContent(ibin)),2);
	//	syserror += pow((scaledo_LepPt->GetBinContent(ibin)-p_LepPt->GetBinContent(ibin)),2);
	//	std::cout << "Pt " << sqrt(syserror)/p_LepPt->GetBinContent(ibin) << std::endl;
	//}	
	for(int ibin(1); ibin < p_MET->GetSize(); ibin++){
		double syserror(0);
		syserror += pow((jesup_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)),2);
		syserror += pow((jesdo_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)),2);
		std::cout << "MET " << sqrt(syserror)/p_MET->GetBinContent(ibin) << std::endl;
	}	
	for(int ibin(1); ibin < p_Mt->GetSize(); ibin++){
		double syserror(0);
		syserror += pow((jesup_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)),2);
		syserror += pow((jesdo_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)),2);
		std::cout << "Mt " << sqrt(syserror)/p_Mt->GetBinContent(ibin) << std::endl;
	}	
	for(int ibin(1); ibin < p_HT->GetSize(); ibin++){
		double syserror(0);
		syserror += pow((jesup_HT->GetBinContent(ibin)-p_HT->GetBinContent(ibin)),2);
		syserror += pow((jesdo_HT->GetBinContent(ibin)-p_HT->GetBinContent(ibin)),2);
		std::cout << "HT " << syserror << std::endl;
	}	
}


