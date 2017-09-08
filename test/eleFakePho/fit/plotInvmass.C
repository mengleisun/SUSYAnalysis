#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include "TH1F.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TF1.h"
#include "TMath.h"
#include "TSystem.h"
#include "TChain.h"
#include "TTree.h"

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
#include "../../../include/RooCMSShape.h"
#include "../../../include/RooDCBShape.h"
#include "../../../include/RooUserPoly.h"
#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_tools.h"
#include "../../../include/analysis_mcData.h"
//#include "../include/RooCBExGaussShape.h"

enum BinType{
  byPt = 0,
  byEta = 1,
  byVtx = 2,
  nonType = 3
};

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

void plotInvmass(int inputbintype, int inputfittype, float lowercut, float uppercut){//main  

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libRooFitClasses.so");
	TFile phoscale("../egammaEffi.txt_EGM2D.root");
	TFile elescale("../scaleFactors.root");
	TH2D  *h_elescale = (TH2D*)elescale.Get("GsfElectronToCutBasedSpring15M");
	TH2D  *h_eleiso   = (TH2D*)elescale.Get("MVAVLooseElectronToMini");
	TH2F  *h_phoscale = (TH2F*)phoscale.Get("EGamma_SF2D");
  
  BinType bintype = (BinType)inputbintype;

//************** Process DYJet Tree *******************************************//
	TH1F  *h_DYinvmass = new TH1F("h_DYinvmass","h_DYinvmass",60,60,120);
	TH1F  *h_DYinvmasscorr = new TH1F("h_DYinvmasscorr","h_DYinvmass",60,60,120);
	h_DYinvmasscorr->SetLineColor(kRed);

  TChain *DYtree = new TChain("FakeRateTree");
  DYtree->Add("/uscms_data/d3/mengleis/plot_elefakepho_DYJetsToLL_2.root");

  float DY_invmass=0; 
  float DY_invmassUncalib=0; 
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
  DYtree->SetBranchAddress("invmass",   &DY_invmass); 
  DYtree->SetBranchAddress("invmassUncalib",   &DY_invmassUncalib); 
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
	  bool isZee(false);
	  double mindRtag(0.3), mindRprobe(0.3);
	  unsigned tagIndex(0), probeIndex(0);
	  for(unsigned iMC(0); iMC<mcPID->size(); iMC++){
			double dR1 = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], DY_tagEta, DY_tagPhi);
			double dR2 = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], DY_probeEta,DY_probePhi);
			if(dR1 < mindRtag){mindRtag=dR1; tagIndex=iMC;}
			if(dR2 < mindRprobe){mindRprobe=dR2; probeIndex=iMC;}
	  }
	  if(mindRtag < 0.1 && mindRprobe < 0.1){
			bool isZe(false),isZg(false);
			isZe = isElectron(fabs((*mcPID)[tagIndex]), fabs((*mcMomPID)[tagIndex]));
			isZg = isElectron(fabs((*mcPID)[probeIndex]), fabs((*mcMomPID)[probeIndex]));
  		if(isZe && isZg)isZee=true; 
    }
    float keyvariable(0);
    if(bintype == byPt)keyvariable = DY_probePt;
    else if(bintype == byEta)keyvariable = fabs(DY_probeEta);
    else if(bintype == byVtx)keyvariable = DY_nVertex;
		int tagbin_x = h_elescale->GetXaxis()->FindBin(DY_tagPt);
		int tagbin_y = h_elescale->GetYaxis()->FindBin(fabs(DY_tagEta));
		int probebin_x=h_phoscale->GetXaxis()->FindBin(DY_probeEta);
		int probebin_y=h_phoscale->GetYaxis()->FindBin(fabs(DY_probePt));
		double eventweight = h_elescale->GetBinContent(h_elescale->GetBin(tagbin_x, tagbin_y))*h_eleiso->GetBinContent(h_eleiso->GetBin(tagbin_x, tagbin_y))*h_phoscale->GetBinContent(h_phoscale->GetBin(probebin_x, probebin_y));
    if(keyvariable >= lowercut && keyvariable < uppercut){
			h_DYinvmass->Fill(DY_invmassUncalib); 
      h_DYinvmasscorr->Fill(DY_invmass, eventweight);
      //if(isZee && inputfittype == 0 && DY_vetovalue== false){ h_DYinvmasscorr->Fill(DY_invmass, eventweight); }
      //else if(isZee && inputfittype == 1 && DY_vetovalue== true){h_DYinvmasscorr->Fill(DY_invmass, eventweight); }
    }
  }
	std::cout << h_DYinvmass->GetEntries() << std::endl;
	TCanvas *can  = new TCanvas("can","can",600,600);
	can->cd();
	h_DYinvmass->Sumw2();
	h_DYinvmasscorr->Sumw2();
	h_DYinvmasscorr->Scale(h_DYinvmass->GetBinContent(31)/h_DYinvmasscorr->GetBinContent(31));
	h_DYinvmass->Draw();
	h_DYinvmasscorr->Draw("same");
	can->SaveAs("Invmass.pdf");
}
