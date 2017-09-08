#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<ctime>

#include "TFile.h"
#include "TTree.h"
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
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TProfile2D.h"

float DeltaR(float eta1,float phi1,float eta2,float phi2)
{
	float deltaPhi = TMath::Abs(phi1-phi2);
	float deltaEta = eta1-eta2;
	if(deltaPhi > TMath::Pi())
	deltaPhi = TMath::TwoPi() - deltaPhi;
		return TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
}
void test_WG(){//main 

	Double_t plotEtBins[]={130,150,200,250,300,350,400,600,800};
	Double_t plotPtBins[]={25,30,35,40,45,50,55,60,65,70,75,80, 85,90,95,100,105,110,115,120,125,130, 135,140,146,152,158,164,170,177,184,192, 200,208,216,224,232,240,250,260,275,290, 305,325,345,370,400};
	TH1F *p_phoEt_WG     = new TH1F("p_phoEt_WG","",8,plotEtBins);
	TH1F *p_mcEt = new TH1F("p_mcEt","",100,0,400);
	TH1F *p_phoEta_WG    = new TH1F("p_phoEta_WG","; |eta|;",30,0,3);
	TH1F *p_phoEta_WG_neg    = new TH1F("p_phoEta_WG_neg","",30,0,3);
	TH1F *p_phoPhi_WG    = new TH1F("p_phoPhi_WG","",64,-3.2, 3.2);
	TH1F *p_lepPt_WG     = new TH1F("p_lepPt_WG","",46,plotPtBins);
	TH1F *p_lepEta_WG    = new TH1F("p_lepEta_WG","",60,-3,3);
	TH1F *p_lepPhi_WG    = new TH1F("p_lepPhi_WG","",64,-3.2, 3.2);
	TH1F *p_sigMET_WG    = new TH1F("p_sigMET_WG","",100,0,400);
	TH1F *p_sigMT_WG     = new TH1F("p_sigMT_WG","",100,0,400);
	TH1F *p_sigMETPhi_WG = new TH1F("p_sigMETPhi_WG","",64,-3.2, 3.2);
	TH1F *p_dPhiLepMET_WG= new TH1F("p_dPhiLepMET_WG","",64,-3.2, 3.2);
	TH1F *p_nVertex_WG   = new TH1F("p_nVertex_WG","",100,0,100);
	TH1F *p_dRPhoLep_WG  = new TH1F("p_dRPhoLep_WG","",60,0,6);
	TH1F *p_HT_WG        = new TH1F("p_HT_WG","",100,0,1000);
	TH1F *p_nJet_WG      = new TH1F("p_nJet_WG","",10,0,10);
	TH1F *p_invmass_WG   = new TH1F("p_invmass_WG","",100,0,200);

	TH1F *p_phoEt_new     = new TH1F("p_phoEt_new","",8,plotEtBins);
	TH1F *p_phoEta_new    = new TH1F("p_phoEta_new","",30,0,3);
	TH1F *p_phoPhi_new    = new TH1F("p_phoPhi_new","",64,-3.2, 3.2);
	TH1F *p_lepPt_new     = new TH1F("p_lepPt_new","",46,plotPtBins);
	TH1F *p_lepEta_new    = new TH1F("p_lepEta_new","",60,-3,3);
	TH1F *p_lepPhi_new    = new TH1F("p_lepPhi_new","",64,-3.2, 3.2);
	TH1F *p_sigMET_new    = new TH1F("p_sigMET_new","",100,0,400);
	TH1F *p_sigMT_new     = new TH1F("p_sigMT_new","",100,0,400);
	TH1F *p_sigMETPhi_new = new TH1F("p_sigMETPhi_new","",64,-3.2, 3.2);
	TH1F *p_dPhiLepMET_new= new TH1F("p_dPhiLepMET_new","",64,-3.2, 3.2);
	TH1F *p_nVertex_new   = new TH1F("p_nVertex_new","",100,0,100);
	TH1F *p_dRPhoLep_new  = new TH1F("p_dRPhoLep_new","",60,0,6);
	TH1F *p_HT_new        = new TH1F("p_HT_new","",100,0,1000);
	TH1F *p_nJet_new      = new TH1F("p_nJet_new","",10,0,10);
	TH1F *p_invmass_new   = new TH1F("p_invmass_new","",100,0,200);

//************ Signal Tree **********************//
  TChain *tree = new TChain("MCTree");
	tree->Add("/uscms_data/d3/mengleis/test/resTree_MC_ZG.root");
	float MCweight;
  std::vector<int>   *mcPID=0;
  std::vector<float> *mcEta=0;
  std::vector<float> *mcPhi=0;
  std::vector<float> *mcPt=0;
  std::vector<int>   *mcMomPID=0;
  std::vector<int>   *mcGMomPID=0;
 
	tree->SetBranchAddress("MCweight", &MCweight); 
  tree->SetBranchAddress("mcPID",    &mcPID);
  tree->SetBranchAddress("mcEta",    &mcEta);
  tree->SetBranchAddress("mcPhi",    &mcPhi);
  tree->SetBranchAddress("mcPt",     &mcPt);
  tree->SetBranchAddress("mcMomPID", &mcMomPID);
  tree->SetBranchAddress("mcGMomPID",&mcGMomPID);

  for(unsigned ievt(0); ievt<tree->GetEntries(); ++ievt){//loop on entries
		tree->GetEntry(ievt);
		double weight = MCweight; 

		for(unsigned it(0); it < mcPID->size(); it++){
		//	if(fabs((*mcMomPID)[it]) < 25){
		  {
				if((*mcPt)[it] < 135)continue;
				p_phoEt_WG->Fill( (*mcPt)[it], weight); 
				p_phoEta_WG->Fill( (*mcEta)[it], weight);
				if((*mcEta)[it] < 0)p_phoEta_WG_neg->Fill(fabs((*mcEta)[it]), weight);
				p_phoPhi_WG->Fill((*mcPhi)[it], weight);
			}
		}
	}//loop on  events




//************ Signal Tree **********************//
  TChain *newtree = new TChain("MCTree");
	newtree->Add("/uscms_data/d3/mengleis/test/resTree_MC_LO130.root");
	float new_MCweight;
  std::vector<int>   *new_mcPID=0;
  std::vector<float> *new_mcEta=0;
  std::vector<float> *new_mcPhi=0;
  std::vector<float> *new_mcPt=0;
  std::vector<int>   *new_mcMomPID=0;
  std::vector<int>   *new_mcGMomPID=0;
 
	newtree->SetBranchAddress("MCweight", &new_MCweight); 
  newtree->SetBranchAddress("mcPID",    &new_mcPID);
  newtree->SetBranchAddress("mcEta",    &new_mcEta);
  newtree->SetBranchAddress("mcPhi",    &new_mcPhi);
  newtree->SetBranchAddress("mcPt",     &new_mcPt);
  newtree->SetBranchAddress("mcMomPID", &new_mcMomPID);
  newtree->SetBranchAddress("mcGMomPID",&new_mcGMomPID);

  for(unsigned ievt(0); ievt<newtree->GetEntries(); ++ievt){//loop on entries
		newtree->GetEntry(ievt);
		double weight = new_MCweight; 

		for(unsigned it(0); it < new_mcPID->size(); it++){
			//if(fabs((*new_mcMomPID)[it]) < 25){
			{
				if( (*new_mcPt)[it] < 135)continue;
				p_phoEt_new->Fill( (*new_mcPt)[it], weight); 
				p_phoEta_new->Fill( (*new_mcEta)[it], weight);
				p_phoPhi_new->Fill((*new_mcPhi)[it], weight);
			}
		}
	}//loop on  events

	
	TCanvas *can_phoEt     = new TCanvas("can_phoEt",       "can_phoEt", 600,600); 
	TCanvas *can_phoEta    = new TCanvas("can_phoEta",      "can_phoEta",600,600); 
	TCanvas *can_phoPhi    = new TCanvas("can_phoPhi",      "can_phoPhi",600,600); 

	//float scalefactor = p_phoEt_WG->Integral(21,24)/p_phoEt_new->Integral(21,24);
	//std::cout << "scale = " << scalefactor<<std::endl;
	float scalefactor = p_phoEt_WG->Integral(1,2)/p_phoEt_new->Integral(1,2);
	//float scalefactor = p_phoEt_WG->GetEntries()/p_phoEt_new->GetEntries(); 
	//float scalefactor = 1;
	can_phoEt->cd();
	gPad->SetLogy();
	p_phoEt_WG->Draw();
	p_phoEt_new->SetLineColor(kRed);
	p_phoEt_new->Scale(scalefactor);
	p_phoEt_new->Draw("same");


	can_phoEta->cd();
	p_phoEta_WG->Draw();
	p_phoEta_WG_neg->SetLineColor(kGreen);
	p_phoEta_WG_neg->Draw("same");
	p_phoEta_new->SetLineColor(kRed);
	p_phoEta_new->Scale(scalefactor);
	p_phoEta_new->Draw("same");

	can_phoPhi->cd();
	p_phoPhi_WG->Draw();
	p_phoPhi_new->SetLineColor(kRed);
	p_phoPhi_new->Scale(scalefactor);
	p_phoPhi_new->Draw("same");

}


