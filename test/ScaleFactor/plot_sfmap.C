#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>

#include "TPaletteAxis.h"
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

void plot_sfmap(){//main  

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");
	gStyle->SetOptStat(0);
	//setTDRStyle();

//	TFile *photonIDFile  = TFile::Open("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/sf/egammaEffi.txt_EGM2D.root");
//	TH2F  *photonIDESF   = (TH2F*)photonIDFile->Get("EGamma_SF2D");
//	photonIDESF->SetTitle("#gamma scale factor");
//	gStyle->SetOptStat(0);
//	TCanvas *c_mt = new TCanvas("Mt", "Mt",600,600);
//	c_mt->SetTopMargin(0.2);
//	c_mt->cd();
//	TPad *pho_pad = new TPad("pho_pad", "pho_pad", 0, 0, 1, 1.0);
//	pho_pad->SetTopMargin(0.1);
//	pho_pad->SetBottomMargin(0.12);
//	pho_pad->SetLeftMargin(0.18);
//	pho_pad->SetRightMargin(0.15);
//	pho_pad->Draw();  
//	pho_pad->cd();  
//	gPad->SetLogy();
//	gStyle->SetPaintTextFormat("4.2f");
//	photonIDESF->SetTitle("#gamma scale factor");
//	photonIDESF->Draw("colz text E");
//	c_mt->SaveAs("SF_Photon.pdf");	
//



//	TFile *electronIDFile = TFile::Open("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/sf/ele_scaleFactors.root");
//	TH2F *electronIDESF = (TH2F*)electronIDFile->Get("GsfElectronToCutBasedSpring15M");
//	TH2F *electronISOESF= (TH2F*)electronIDFile->Get("MVAVLooseElectronToMini");
//	TCanvas *c_eleid = new TCanvas("Mt", "Mt",1200,600);
//	c_eleid->SetTopMargin(0.2);
//	c_eleid->Divide(2);
//	c_eleid->cd(1);
//	TPad *eleid_pad = new TPad("eleid_pad", "eleid_pad", 0, 0, 1, 1.0);
//	eleid_pad->SetTopMargin(0.1);
//	eleid_pad->SetBottomMargin(0.12);
//	eleid_pad->SetLeftMargin(0.18);
//	eleid_pad->SetRightMargin(0.15);
//	eleid_pad->Draw();  
//	eleid_pad->cd();  
//	gPad->SetLogx();
//	gStyle->SetPaintTextFormat("4.2f");
//	electronIDESF->SetTitle("e ID scale factor");
//	electronIDESF->Draw("colz text E");
//	c_eleid->cd(2);
//	TPad *eleiso_pad = new TPad("eleiso_pad", "eleiso_pad", 0, 0, 1, 1.0);
//	eleiso_pad->SetTopMargin(0.1);
//	eleiso_pad->SetBottomMargin(0.12);
//	eleiso_pad->SetLeftMargin(0.18);
//	eleiso_pad->SetRightMargin(0.15);
//	eleiso_pad->Draw();  
//	eleiso_pad->cd();  
//	gPad->SetLogx();
//	gStyle->SetPaintTextFormat("4.2f");
//	electronISOESF->SetTitle("e mini-Isolation scale factor");
//	electronISOESF->Draw("colz text E");
//	c_eleid->SaveAs("SF_Ele.pdf");	
	
	TFile *muonIDFile     = TFile::Open("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/sf/TnP_NUM_MediumID_DENOM_generalTracks_VAR_map_pt_eta.root");
	TFile *muonIsoFile    = TFile::Open("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/sf/TnP_NUM_MiniIsoTight_DENOM_MediumID_VAR_map_pt_eta.root");
	TH2F *muonIDESF    = (TH2F*)muonIDFile->Get("SF");
	TH2F *muonISOESF   = (TH2F*)muonIsoFile->Get("SF");
	TCanvas *c_muid = new TCanvas("Mt", "Mt",1200,600);
	c_muid->SetTopMargin(0.2);
	c_muid->SetRightMargin(0.1);
	c_muid->Divide(2);
	c_muid->cd(1);
	TPad *muid_pad = new TPad("muid_pad", "muid_pad", 0, 0, 1, 1.0);
	muid_pad->SetTopMargin(0.1);
	muid_pad->SetBottomMargin(0.12);
	muid_pad->SetLeftMargin(0.18);
	muid_pad->SetRightMargin(0.1);
	muid_pad->Draw();  
	muid_pad->cd();  
	gPad->SetLogx();
	gStyle->SetPaintTextFormat("4.3f");
	muonIDESF->SetTitle("muon ID scale factor");
	muonIDESF->GetZaxis()->SetTitle("");
  TPaletteAxis *paletteID = (TPaletteAxis*)muonIDESF->GetListOfFunctions()->FindObject("palette");
  paletteID->SetX1NDC(0.9);
  paletteID->SetX2NDC(0.95);
  c_muid->Modified();
  c_muid->Update();
	muonIDESF->Draw("col text E");
	c_muid->cd(2);
	TPad *muiso_pad = new TPad("muiso_pad", "muiso_pad", 0, 0, 1, 1.0);
	muiso_pad->SetTopMargin(0.1);
	muiso_pad->SetBottomMargin(0.12);
	muiso_pad->SetLeftMargin(0.18);
	muiso_pad->SetRightMargin(0.2);
	muiso_pad->Draw();  
	muiso_pad->cd();  
	gPad->SetLogx();
	gStyle->SetPaintTextFormat("4.3f");
	muonISOESF->SetTitle("muon mini-Isolation scale factor");
	muonISOESF->GetZaxis()->SetTitle("");
  TPaletteAxis *palette = (TPaletteAxis*)muonISOESF->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.82);
  palette->SetX2NDC(0.88);
  c_muid->Modified();
  c_muid->Update();
	muonISOESF->Draw("colz text E");
	c_muid->SaveAs("SF_Muon.pdf");	
}


