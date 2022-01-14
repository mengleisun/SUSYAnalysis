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

void plotTrigger(){//main  

	std::ostringstream treename;
	bool plotLeading(false);

	treename.str("");
	if(plotLeading)treename << "egTree";
	else treename << "eeTree";

	gStyle->SetOptStat(0);
	setTDRStyle();    
	gStyle->SetErrorX(0.5);
	TFile *file = TFile::Open("/uscms_data/d3/mengleis/FullStatusOct/plot_egTrigger_ReMiniAOD.root");

  TCanvas *can[8];
  std::ostringstream canvas; 
  for(unsigned iC(0); iC<8; iC++){
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

  Double_t plotPtBins[]={10,11,12,13,14,14,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,34,36,38,40,44,48,52,56,60,70,100};
  Double_t plotPt2DBins[]={10,20,35,50,90,150,500};
  Double_t plotLeadEtaBins[]={-1.444, -0.8, 0, 0.8, 1.444};
  Double_t plotTrailEtaBins[]={-2.5,-2.0,-1.566,-1.444, -0.8, 0, 0.8, 1.444, 1.566, 2.0, 2.5};
  TProfile *p_leadeffEB  = new TProfile("p_leadeffEB", ";p_{T}(GeV);Efficiency",  33, plotPtBins); 
  TProfile *p_leadeffEE  = new TProfile("p_leadeffEE", ";p_{T}(GeV);Efficiency",  33, plotPtBins); 
  TProfile *p_traileffEB = new TProfile("p_traileffEB",";p_{T}(GeV);Efficiency",  33, plotPtBins); 
  TProfile *p_traileffEE = new TProfile("p_traileffEE",";p_{T}(GeV);Efficiency",  33, plotPtBins); 

	TProfile2D *p_effEB2D;
  p_effEB2D  = new TProfile2D("Traileff_data", "",  6, plotPt2DBins, 10, plotTrailEtaBins); 

  for(unsigned ievt(0); ievt < egtree->GetEntries(); ievt++){
    egtree->GetEntry(ievt);
		if(invmass < 80 || invmass > 101)continue;
		float photonEt(0);
		if(probeEt < 500)photonEt = probeEt;
		else photonEt = 499.5;

		if(fabs(probeEta) < 1.444 && probeR9 > 0.5 ){
			if(probeMatchLeading && probeMatchTrailing)p_leadeffEB->Fill(photonEt, 1);
			else p_leadeffEB->Fill(photonEt, 0);
			if(probeMatchTrailing)p_traileffEB->Fill(photonEt, 1);
			else p_traileffEB->Fill(photonEt, 0);
		}
		else if(fabs(probeEta) > 1.56 && probeR9 > 0.8){
			if(probeMatchLeading && probeMatchTrailing)p_leadeffEE->Fill(photonEt, 1);
			else p_leadeffEE->Fill(photonEt, 0);
			if(probeMatchTrailing)p_traileffEE->Fill(photonEt, 1);
			else p_traileffEE->Fill(photonEt, 0);
		}
	
		if(plotLeading){	
			if(probeMatchLeading && probeMatchTrailing)p_effEB2D->Fill(photonEt, probeEta, 1);
			else p_effEB2D->Fill(photonEt, probeEta, 0);
		}
		else{
			if(probeMatchTrailing)p_effEB2D->Fill(photonEt, probeEta, 1);
			else p_effEB2D->Fill(photonEt, probeEta, 0);
		}
  }


/*********************  Fitting result ************************/
	std::ifstream FittingResult("TriggerResult_eg/egTrigger-Trailing-Fit.txt");
	int nEtaBins(10);
	double pt_den[6][nEtaBins];
	double pt_num[6][nEtaBins];
	double fit_eff[6][nEtaBins];
	std::string bintype;
	double ptvalue(0), etavalue(0);
	double fitmean(0), fitrms(0);
	if(FittingResult.is_open()){
		for(int i(1); i<=6; i++){
			for(int j(1); j <= nEtaBins; j++){
				FittingResult >> bintype >> ptvalue >> etavalue >> fitmean >> fitrms;
				pt_den[i-1][j-1] = fitmean;
			}
		}
		for(int i(1); i<=6; i++){
			for(int j(1); j <= nEtaBins; j++){
				FittingResult >> bintype >> ptvalue >> etavalue >> fitmean >> fitrms;
				pt_num[i-1][j-1] = fitmean;
				fit_eff[i-1][j-1]= pt_num[i-1][j-1]/pt_den[i-1][j-1];
			}
		}
	}
/***************************************************************************  Start MC ***************************************************************************/
	TFile *mcfile = TFile::Open("/eos/uscms/store/user/tmishra/Trigger/plot_egTrigger_DY.root");

	TTree *DYtree = (TTree*)mcfile->Get(treename.str().c_str());
	float DY_tagPt(0);
	float DY_tagEta(0);
	float DY_tagPhi(0);
	float DY_tagR9(0);
	float DY_probeEt(0);
	float DY_probeEta(0);
	float DY_probePhi(0);
	float DY_probeR9(0);
	bool  DY_probeMatchLeading;
	bool  DY_probeMatchTrailing;
	float DY_invmass;
	std::vector<int>   *mcPID=0;
	std::vector<float> *mcEta=0;
	std::vector<float> *mcPhi=0;
	std::vector<float> *mcPt=0;
	std::vector<int>   *mcMomPID=0;
	std::vector<int>   *mcGMomPID=0;

	DYtree->SetBranchAddress("tagPt",                 &DY_tagPt);
	DYtree->SetBranchAddress("tagEta",                &DY_tagEta);
	DYtree->SetBranchAddress("tagPhi",                &DY_tagPhi);
	DYtree->SetBranchAddress("tagR9",                 &DY_tagR9);
	DYtree->SetBranchAddress("probeEt",            		&DY_probeEt);
	DYtree->SetBranchAddress("probeEta",           		&DY_probeEta);
	DYtree->SetBranchAddress("probePhi",           		&DY_probePhi);
	DYtree->SetBranchAddress("probeR9",            		&DY_probeR9);
	DYtree->SetBranchAddress("probeMatchLeading",  		&DY_probeMatchLeading);
	DYtree->SetBranchAddress("probeMatchTrailing", 		&DY_probeMatchTrailing);
	DYtree->SetBranchAddress("invmass",               &DY_invmass);
	DYtree->SetBranchAddress("mcPID",			   					&mcPID);
	DYtree->SetBranchAddress("mcEta",			   					&mcEta);
	DYtree->SetBranchAddress("mcPhi",			   					&mcPhi);
	DYtree->SetBranchAddress("mcPt",				   				&mcPt);
	DYtree->SetBranchAddress("mcMomPID",			   			&mcMomPID);
	DYtree->SetBranchAddress("mcGMomPID",		   				&mcGMomPID);

	TProfile2D *p_effEBMC;
  p_effEBMC  = new TProfile2D("Traileff_MC", "",  6, plotPt2DBins, 10, plotTrailEtaBins); 

  for(unsigned ievt(0); ievt < DYtree->GetEntries(); ievt++){
    DYtree->GetEntry(ievt);
		float photonEt(0);
		if(probeEt < 500)photonEt = DY_probeEt;
		else photonEt = 499.5;
		
	  bool isZee(false);
	  double mindRtag(0.3), mindRprobe(0.3);
	  unsigned tagIndex(0), probeIndex(0);
	  for(unsigned iMC(0); iMC<mcPID->size(); iMC++){
			double dR1 = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], DY_tagEta, DY_tagPhi);
			double dR2 = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], DY_probeEta,DY_probePhi);
			double dE1 = fabs((*mcPt)[iMC] - DY_tagPt)/DY_tagPt;
			double dE2 = fabs((*mcPt)[iMC] - photonEt)/photonEt;
			if(dR1 < mindRtag && dE1 < 0.1){mindRtag=dR1; tagIndex=iMC;}
			if(dR2 < mindRprobe && dE2 < 0.1){mindRprobe=dR2; probeIndex=iMC;}
	  }
	  if(mindRtag < 0.1 && mindRprobe < 0.1){
			bool isZe(false),isZg(false);
			isZe = isElectron(fabs((*mcPID)[tagIndex]), fabs((*mcMomPID)[tagIndex]));
			isZg = isElectron(fabs((*mcPID)[probeIndex]), fabs((*mcMomPID)[probeIndex]));
  		if(isZe && isZg)isZee=true; 
    }

		if(plotLeading){
			if(isZee){
				if(DY_probeMatchLeading && DY_probeMatchTrailing)p_effEBMC->Fill(photonEt, DY_probeEta, 1);
				else p_effEBMC->Fill(photonEt, DY_probeEta, 0);
			}
		}
		else{
			if(isZee){
				if(DY_probeMatchTrailing)p_effEBMC->Fill(photonEt, DY_probeEta, 1);
				else p_effEBMC->Fill(photonEt, DY_probeEta, 0);
			}
		}
  }

	TH2F *p_TriggerScale;
	if(plotLeading)p_TriggerScale = new TH2F("ESF_Lead","photon trigger ESF;p_{T} (GeV); #eta",6, plotPt2DBins, 10, plotTrailEtaBins); 
	else p_TriggerScale = new TH2F("p_TriggerScale","photon trigger ESF;p_{T} (GeV); #eta",6, plotPt2DBins, 10, plotTrailEtaBins);
	for(int i(1); i<=6; i++){
		for(int j(1); j <= 10; j++){
			p_TriggerScale->SetBinContent(i, j,p_effEB2D->GetBinContent(i,j)/p_effEBMC->GetBinContent(i,j));
			double staterror = p_effEB2D->GetBinError(i,j);
			double syserror = i<6? fabs(fit_eff[i-1][j-1]-p_effEB2D->GetBinContent(i,j)):0 ;
			p_TriggerScale->SetBinError(i, j, sqrt(staterror*staterror + syserror*syserror));
		}
	}

	std::ostringstream outputname;
	outputname.str("");
	if(plotLeading)outputname << "diphoton_pholeg.root";
	else outputname << "diphoton_eleg.root";
	TFile *outputfile = TFile::Open(outputname.str().c_str(),"RECREATE");
	outputfile->cd();
	p_effEB2D->Write();
	p_effEBMC->Write();
	p_TriggerScale->Write();
	outputfile->Write();
	outputfile->Close();

	gStyle->SetPaintTextFormat("4.4f");
	Int_t PaletteColors[] = {9, kBlue, kBlue-4,kCyan, kTeal, kGreen,kSpring, 5, 2};
	gStyle->SetPalette(9, PaletteColors);
	can[0]->cd();
	p_leadeffEB->Draw();
	if(plotLeading)can[0]->SaveAs("egTrigger_LeadingEB_re.pdf");

	can[1]->cd();
	p_traileffEB->Draw();
	if(!plotLeading)can[1]->SaveAs("egTrigger_TrailingEB_re.pdf");
	
	can[2]->cd();
	p_traileffEE->Draw();
	if(!plotLeading)can[2]->SaveAs("egTrigger_TrailingEE_re.pdf");

	can[3]->cd();
	p_traileffEB->SetLineColor(kBlack);
	p_traileffEB->Draw();
	p_traileffEE->SetLineColor(kBlue);
	p_traileffEE->SetMarkerColor(kBlue);
	p_traileffEE->Draw("same");
	TLegend *leg = new TLegend(0.5,0.2,0.8,0.4);
	leg->AddEntry(p_traileffEB, "EB");
	leg->AddEntry(p_traileffEE, "EE");
	leg->Draw("same");
	if(!plotLeading)can[3]->SaveAs("egTrigger_TrailingAll_re.pdf");

  can[4]->cd();
	p_effEB2D->Draw("colz text");

  can[5]->cd();
	p_effEBMC->Draw("colz text");

  can[7]->cd();
	gPad->SetLogx();
	if(!plotLeading)p_TriggerScale->SetTitle("electron trigger ESF");
	p_TriggerScale->GetXaxis()->SetRangeUser(20,500);
	p_TriggerScale->Draw("colz E text");
	if(plotLeading)can[7]->SaveAs("egTrigger_LeadingESF_re.pdf");
	else can[7]->SaveAs("egTrigger_TrailingESF_re.pdf");

	for(unsigned ibin(1); ibin < p_leadeffEB->GetSize(); ibin++){
		std::cout << p_leadeffEB->GetBinCenter(ibin) << " EB " << p_leadeffEB->GetBinContent(ibin) << " EE " << p_traileffEE->GetBinContent(ibin) << std::endl;
	}
}


