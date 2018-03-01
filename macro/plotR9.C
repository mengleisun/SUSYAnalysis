#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TFrame.h"
#include "TF1.h"
#include "TMath.h"
#include "TProfile2D.h"
#include "TPad.h"
#include "TPaveText.h"
#include "../include/tdrstyle.C"

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

float DeltaR(float eta1,float phi1,float eta2,float phi2)
{
	float deltaPhi = TMath::Abs(phi1-phi2);
	float deltaEta = eta1-eta2;
	if(deltaPhi > TMath::Pi())
	deltaPhi = TMath::TwoPi() - deltaPhi;
		return TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
}

void plotR9(){//main  

	std::ostringstream treename;
	bool plotLeading(false);

	treename.str("");
	if(plotLeading)treename << "egTree";
	else treename << "eeTree";

	gStyle->SetOptStat(0);
	setTDRStyle();    
	gStyle->SetErrorX(0.5);
	TFile *file = TFile::Open("/uscms_data/d3/mengleis/FullStatusOct/plot_R9_DY.root");

  TCanvas *can[3];
  std::ostringstream canvas; 
  for(unsigned iC(0); iC<3; iC++){
	 canvas.str("");
	 canvas << "canvas" << iC;
	 can[iC] = new TCanvas(canvas.str().c_str(),canvas.str().c_str(),600,600);
	 setCanvas(can[iC]);
	 can[iC]->SetLeftMargin(0.15);
	 can[iC]->SetBottomMargin(0.15);
  }
	TLatex cms;
	cms.SetTextSize(0.04);

	TTree *egtree = (TTree*)file->Get(treename.str().c_str());
	float tagPt(0);
	float tagEta(0);
	float tagPhi(0);
	float tagR9(0);
	float probeEt(0);
	float probeEta(0);
	float probePhi(0);
	float probeR9(0);
	bool  probeMatchLeading;
	bool  probeMatchTrailing;
	float invmass;
	std::vector<int>   *mcPID=0;
	std::vector<float> *mcEta=0;
	std::vector<float> *mcPhi=0;
	std::vector<float> *mcPt=0;
	std::vector<int>   *mcMomPID=0;
	std::vector<int>   *mcGMomPID=0;

	egtree->SetBranchAddress("tagPt",                 &tagPt);
	egtree->SetBranchAddress("tagEta",                &tagEta);
	egtree->SetBranchAddress("tagPhi",                &tagPhi);
	egtree->SetBranchAddress("tagR9",                 &tagR9);
	egtree->SetBranchAddress("probeEt",            		&probeEt);
	egtree->SetBranchAddress("probeEta",           		&probeEta);
	egtree->SetBranchAddress("probePhi",           		&probePhi);
	egtree->SetBranchAddress("probeR9",            		&probeR9);
	egtree->SetBranchAddress("probeMatchLeading",  		&probeMatchLeading);
	egtree->SetBranchAddress("probeMatchTrailing", 		&probeMatchTrailing);
	egtree->SetBranchAddress("invmass",               &invmass);
	egtree->SetBranchAddress("mcPID",			   					&mcPID);
	egtree->SetBranchAddress("mcEta",			   					&mcEta);
	egtree->SetBranchAddress("mcPhi",			   					&mcPhi);
	egtree->SetBranchAddress("mcPt",				   				&mcPt);
	egtree->SetBranchAddress("mcMomPID",			   			&mcMomPID);
	egtree->SetBranchAddress("mcGMomPID",		   				&mcGMomPID);

  Double_t plotPtBins[]={10,11,12,13,14,14,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,34,36,38,40,44,48,52,56,60,70,100};
  Double_t plotPt2DBins[]={10,20,30,40,50,200};
  Double_t plotLeadEtaBins[]={-1.444, -0.8, 0, 0.8, 1.444};
  Double_t plotTrailEtaBins[]={-2.5,-2.0,-1.566,-1.444, -0.8, 0, 0.8, 1.444, 1.566, 2.0, 2.5};
  TProfile *p_effEB  = new TProfile("p_effEB", ";p_{T}(GeV);Efficiency",  33, plotPtBins); 
  TProfile *p_effEE  = new TProfile("p_effEE", ";p_{T}(GeV);Efficiency",  33, plotPtBins); 

	TProfile2D *p_effEB2D;
  p_effEB2D  = new TProfile2D("Traileff_data", "",  5, plotPt2DBins, 10, plotTrailEtaBins); 
	TH2D *p_SF = new TH2D("SF","",5, plotPt2DBins, 10, plotTrailEtaBins);

  for(unsigned ievt(0); ievt < egtree->GetEntries(); ievt++){
    egtree->GetEntry(ievt);
		if(invmass < 80 || invmass > 101)continue;
		//if(invmass < 60 || invmass > 120)continue;
	  bool isZee(false);
	  double mindRtag(0.3), mindRprobe(0.3);
	  unsigned tagIndex(0), probeIndex(0);
	  for(unsigned iMC(0); iMC<mcPID->size(); iMC++){
			double dR1 = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], tagEta, tagPhi);
			double dR2 = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], probeEta,probePhi);
			double dE1 = fabs((*mcPt)[iMC] - tagPt)/tagPt;
			double dE2 = fabs((*mcPt)[iMC] - probeEt)/probeEt;
			if(dR1 < mindRtag && dE1 < 0.1){mindRtag=dR1; tagIndex=iMC;}
			if(dR2 < mindRprobe && dE2 < 0.1){mindRprobe=dR2; probeIndex=iMC;}
	  }
	  if(mindRtag < 0.1 && mindRprobe < 0.1){
			bool isZe(false),isZg(false);
			isZe = isElectron(fabs((*mcPID)[tagIndex]), fabs((*mcMomPID)[tagIndex]));
			isZg = isElectron(fabs((*mcPID)[probeIndex]), fabs((*mcMomPID)[probeIndex]));
  		if(isZe && isZg)isZee=true; 
    }
		bool mcTrue(false);
		if(isZee)mcTrue = true;
		else mcTrue = false;
		//if(!mcTrue)continue;

		if(!isZee){
			if(mindRtag < 0.1 && mindRprobe < 0.1){
				std::cout << "tag " << fabs((*mcPID)[tagIndex]) << " " << fabs((*mcMomPID)[tagIndex]) << std::endl;
				std::cout << "pro " << fabs((*mcPID)[probeIndex]) << " " << fabs((*mcMomPID)[probeIndex]) << std::endl;
			}
		}

		float photonEt(0);
		if(probeEt < 500)photonEt = probeEt;
		else photonEt = 499.5;

		if(fabs(probeEta) < 1.444){
			if(probeR9 > 0.5)p_effEB->Fill(photonEt, 1);
			else p_effEB->Fill(photonEt, 0);
		}
		else if(fabs(probeEta) > 1.56){
			if(probeR9 > 0.8)p_effEE->Fill(photonEt, 1);
			else p_effEE->Fill(photonEt, 0);
		}

		if(fabs(probeEta) < 1.444){
			if(probeR9 > 0.5)p_effEB2D->Fill(photonEt,probeEta, 1);
			else p_effEB2D->Fill(photonEt,probeEta, 0);
		}
		else if(fabs(probeEta) > 1.56){
			if(probeR9 > 0.8)p_effEB2D->Fill(photonEt,probeEta, 1);
			else p_effEB2D->Fill(photonEt,probeEta, 0);
		}
  }

	for(unsigned i(1); i<6; i++)
		for(unsigned j(1); j<=10; j++)
			std::cout << "DY " << p_effEB2D->GetBinContent(i,j) << std::endl;

	gStyle->SetPaintTextFormat("4.4f");
	Int_t PaletteColors[] = {9, kBlue, kBlue-4,kCyan, kTeal, kGreen,kSpring, 5, 2};
	gStyle->SetPalette(9, PaletteColors);
	can[0]->cd();
	p_effEB->Draw();

	can[1]->cd();
	p_effEE->Draw();

	std::ostringstream outputname;
	outputname.str("");
	if(plotLeading)outputname << "diphoton_phoR9_DY.root";
	else outputname << "diphoton_eR9_DY.root";
	TFile *outputfile = TFile::Open(outputname.str().c_str(),"RECREATE");
	outputfile->cd();
	p_effEB->Write();
	p_effEE->Write();
	outputfile->Write();
	outputfile->Close();

/*********************  Fitting result ************************/
	//std::ifstream FittingResult("eeR960-120-FullIntegral.txt");
	std::ifstream FittingResult("eeR960-120.txt");
	int nEtaBins(10);
	double pt_den[5][nEtaBins];
	double pt_num[5][nEtaBins];
	double fit_eff[5][nEtaBins];
	std::string bintype;
	double ptvalue(0), etavalue(0);
	double fitmean(0), fitrms(0);
	if(FittingResult.is_open()){
		for(int i(1); i<=5; i++){
			for(int j(1); j <= nEtaBins; j++){
				FittingResult >> bintype >> ptvalue >> etavalue >> fitmean >> fitrms;
				pt_den[i-1][j-1] = fitmean;
			}
		}
		for(int i(1); i<=5; i++){
			for(int j(1); j <= nEtaBins; j++){
				FittingResult >> bintype >> ptvalue >> etavalue >> fitmean >> fitrms;
				pt_num[i-1][j-1] = fitmean;
				fit_eff[i-1][j-1]= pt_num[i-1][j-1]/pt_den[i-1][j-1];
				std::cout << fit_eff[i-1][j-1] << std::endl;
			}
		}
	}
		for(int i(1); i<=5; i++){
			for(int j(1); j <= nEtaBins; j++){
				std::cout << "ratio " << p_effEB2D->GetBinContent(i,j)/fit_eff[i-1][j-1] << std::endl;
				p_SF->SetBinContent( i,j, fit_eff[i-1][j-1]/p_effEB2D->GetBinContent(i,j));
			}
		}
		p_SF->Draw("colz text");
}


