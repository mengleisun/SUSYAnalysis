#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<ctime>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
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

#include "../../include/analysis_rawData.h"
#include "../../include/analysis_photon.h"
#include "../../include/analysis_muon.h"
#include "../../include/analysis_ele.h"
#include "../../include/analysis_mcData.h"
#include "../../include/analysis_tools.h"
#include "../../include/analysis_fakes.h"
#include "../../include/analysis_scalefactor.h"

void test_ZG(){//main 

	int channelType = 2;
  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

	esfScaleFactor  objectESF;


	TH1D *p_iso_data = new TH1D("p_iso_data","p_iso_data",130,0,1.3);		
	TH1D *p_iso_ZG = new TH1D("p_iso_ZG","p_iso_ZG",130,0,1.3);		
	TH1D *p_status = new TH1D("p_status","p_status",7,0,7);

	Double_t plotEtBins[]={35,50,100,150,200,250,300,400,600,800};
	Double_t plotPtBins[]={25,30,35,40,45,50,55,60,65,70,75,80, 85,90,95,100,105,110,115,120,125,130, 135,140,146,152,158,164,170,177,184,192, 200,208,216,224,232,240,250,260,275,290, 305,325,345,370,400};
	TH1D *p_phoEt_data     = new TH1D("p_phoEt_data","",9,plotEtBins);
	TH1D *p_mcEt = new TH1D("p_mcEt","",100,0,400);
	TH1D *p_phoEta_data    = new TH1D("p_phoEta_data","",60,-3,3);
	TH1D *p_phoPhi_data    = new TH1D("p_phoPhi_data","",64,-3.2, 3.2);
	TH1D *p_lepPt_data     = new TH1D("p_lepPt_data","",9,plotEtBins);
	TH1D *p_trailPt_data     = new TH1D("p_trail_data","",9,plotEtBins);
	TH1D *p_lepEta_data    = new TH1D("p_lepEta_data","",60,-3,3);
	TH1D *p_lepPhi_data    = new TH1D("p_lepPhi_data","",64,-3.2, 3.2);
	TH1D *p_sigMET_data    = new TH1D("p_sigMET_data","",100,0,400);
	TH1D *p_sigMT_data     = new TH1D("p_sigMT_data","",100,0,400);
	TH1D *p_sigMETPhi_data = new TH1D("p_sigMETPhi_data","",64,-3.2, 3.2);
	TH1D *p_dPhiLepMET_data= new TH1D("p_dPhiLepMET_data","",64,-3.2, 3.2);
	TH1D *p_nVertex_data   = new TH1D("p_nVertex_data","",100,0,100);
	TH1D *p_dRPhoLep_data  = new TH1D("p_dRPhoLep_data","",60,0.3,3.3);
	TH1D *p_HT_data        = new TH1D("p_HT_data","",100,0,1000);
	TH1D *p_nJet_data      = new TH1D("p_nJet_data","",10,-0.5,9.5);
	TH1D *p_invmass_data   = new TH1D("p_invmass_data","",100,30,130);
	TH1D *p_llmass_data   = new TH1D("p_llmass_data","",100,30,130);

	TH1D *p_phoEt_ZG     = new TH1D("p_phoEt_ZG","",9,plotEtBins);
	TH1D *p_phoEta_ZG    = new TH1D("p_phoEta_ZG","",60,-3,3);
	TH1D *p_phoPhi_ZG    = new TH1D("p_phoPhi_ZG","",64,-3.2, 3.2);
	TH1D *p_lepPt_ZG     = new TH1D("p_lepPt_ZG","",9,plotEtBins);
	TH1D *p_trailPt_ZG     = new TH1D("p_trail_ZG","",9,plotEtBins);
	TH1D *p_lepEta_ZG    = new TH1D("p_lepEta_ZG","",60,-3,3);
	TH1D *p_lepPhi_ZG    = new TH1D("p_lepPhi_ZG","",64,-3.2, 3.2);
	TH1D *p_sigMET_ZG    = new TH1D("p_sigMET_ZG","",100,0,400);
	TH1D *p_sigMT_ZG     = new TH1D("p_sigMT_ZG","",100,0,400);
	TH1D *p_sigMETPhi_ZG = new TH1D("p_sigMETPhi_ZG","",64,-3.2, 3.2);
	TH1D *p_dPhiLepMET_ZG= new TH1D("p_dPhiLepMET_ZG","",64,-3.2, 3.2);
	TH1D *p_nVertex_ZG   = new TH1D("p_nVertex_ZG","",100,0,100);
	TH1D *p_dRPhoLep_ZG  = new TH1D("p_dRPhoLep_ZG","",60,0.3,3.3);
	TH1D *p_HT_ZG        = new TH1D("p_HT_ZG","",100,0,1000);
	TH1D *p_nJet_ZG      = new TH1D("p_nJet_ZG","",10,-0.5,9.5);
	TH1D *p_invmass_ZG   = new TH1D("p_invmass_ZG","",20,80,100);
	TH1D *p_llmass_ZG   = new TH1D("p_llmass_ZG","",100,30,130);

	TProfile *p_scalefactor  = new TProfile("p_scalefactor","p_scalefactor",50,-2.5,2.5);
//************ Signal Tree **********************//
  TChain *tree = new TChain("ZTree");
  tree->Add("/uscms_data/d3/mengleis/test/resTree_ZISR_DY.root");
  float phoEt(0);
  float phoEta(0);
  float phoPhi(0);
	float phoChIso(0);
	float phoNeuIso(0);
  float lepPt(0);
  float lepEta(0);
  float lepPhi(0);
	float trailPt(0);
  float sigMT(0);
  float sigMET(0);
  float sigMETPhi(0);
  float dPhiLepMET(0);
  int   nVertex(0);
  float dRPhoLep(0);
  float HT(0);
  float nJet(0);
	float threeMass(0);
	float dilepMass(0);

//  tree->SetBranchAddress("phoEt",     &phoEt);
//  tree->SetBranchAddress("phoEta",    &phoEta);
//  tree->SetBranchAddress("phoPhi",    &phoPhi);
//	tree->SetBranchAddress("phoChIso",  &phoChIso);
//	tree->SetBranchAddress("phoNeuIso", &phoNeuIso);
  tree->SetBranchAddress("lepPt",     &lepPt);
  tree->SetBranchAddress("lepEta",    &lepEta);
  tree->SetBranchAddress("lepPhi",    &lepPhi);
  tree->SetBranchAddress("trailPt",   &trailPt);
  tree->SetBranchAddress("sigMT",     &sigMT);
  tree->SetBranchAddress("sigMET",    &sigMET);
  tree->SetBranchAddress("sigMETPhi", &sigMETPhi);
  tree->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
  tree->SetBranchAddress("nVertex",   &nVertex);
//  tree->SetBranchAddress("dRPhoLep",  &dRPhoLep);
  tree->SetBranchAddress("HT",        &HT);
  tree->SetBranchAddress("nJet",      &nJet);
//	tree->SetBranchAddress("threeMass", &threeMass);
	tree->SetBranchAddress("dilepMass", &dilepMass);
	
  for(unsigned ievt(0); ievt<tree->GetEntries(); ++ievt){//loop on entries
		tree->GetEntry(ievt);

		if(fabs(nJet-0) >= 0.1)continue;
		if(dilepMass < 80 || dilepMass > 100)continue;
		//double weight = 35.8*1000*4895.0/122053259; 
		double weight = 1; 

		if(fabs(nJet-0) < 0.1)sum_data_1+=weight;
		if( dilepMass >= 80 && dilepMass <= 100)sum_data_2 +=weight;
	
		p_iso_data->Fill(phoNeuIso, weight);			
		p_phoEt_data->Fill(phoEt, weight); 
		p_phoEta_data->Fill(phoEta, weight);
		p_phoPhi_data->Fill(phoPhi, weight);
		p_lepPt_data->Fill(lepPt, weight);
		p_lepEta_data->Fill(lepEta, weight);
		p_lepPhi_data->Fill(lepPhi, weight);
		p_trailPt_data->Fill(trailPt, weight);
		p_sigMET_data->Fill(sigMET, weight);
		p_sigMT_data->Fill(sigMT, weight);
		p_sigMETPhi_data->Fill(sigMETPhi, weight);
		p_dPhiLepMET_data->Fill(dPhiLepMET, weight);
		p_nVertex_data->Fill(nVertex, weight);
		p_dRPhoLep_data->Fill(dRPhoLep, weight);
		p_HT_data->Fill(HT, weight);
		p_nJet_data->Fill(nJet, weight);
		p_invmass_data->Fill(91.5, 1.44);
		p_llmass_data->Fill(dilepMass, weight);
	}//loop on  events




//************ Signal Tree **********************//
  TChain *ZGtree = new TChain("ZTree");
  ZGtree->Add("/uscms_data/d3/mengleis/test/resTree_ZISR_DYLO.root");
  float ZG_phoEt(0);
  float ZG_phoEta(0);
  float ZG_phoPhi(0);
	float ZG_phoChIso(0);
	float ZG_phoNeuIso(0);
  float ZG_lepPt(0);
  float ZG_lepEta(0);
  float ZG_lepPhi(0);
	float ZG_trailPt(0);
	float ZG_trailEta(0);
	float ZG_trailPhi(0);
  float ZG_sigMT(0);
  float ZG_sigMET(0);
  float ZG_sigMETPhi(0);
  float ZG_dPhiLepMET(0);
  int   ZG_nVertex(0);
  float ZG_dRPhoLep(0);
  float ZG_dRPhoTrail(0);
  float ZG_HT(0);
  float ZG_nJet(0);
	float ZG_threeMass(0);
	float ZG_dilepMass(0);
  std::vector<int>   *ZG_mcPID=0;
  std::vector<float> *ZG_mcEta=0;
  std::vector<float> *ZG_mcPhi=0;
  std::vector<float> *ZG_mcPt=0;
  std::vector<int>   *ZG_mcMomPID=0;
  std::vector<int>   *ZG_mcGMomPID=0;
	std::vector<int>   *ZG_mcStatus=0;

 // ZGtree->SetBranchAddress("phoEt",     &ZG_phoEt);
 // ZGtree->SetBranchAddress("phoEta",    &ZG_phoEta);
 // ZGtree->SetBranchAddress("phoPhi",    &ZG_phoPhi);
 // ZGtree->SetBranchAddress("phoChIso",  &ZG_phoChIso);
 // ZGtree->SetBranchAddress("phoNeuIso", &ZG_phoNeuIso);
  ZGtree->SetBranchAddress("lepPt",     &ZG_lepPt);
  ZGtree->SetBranchAddress("lepEta",    &ZG_lepEta);
  ZGtree->SetBranchAddress("lepPhi",    &ZG_lepPhi);
	ZGtree->SetBranchAddress("trailPt",   &ZG_trailPt);
//	ZGtree->SetBranchAddress("trailEta",  &ZG_trailEta);
//	ZGtree->SetBranchAddress("trailPhi",  &ZG_trailPhi);
  ZGtree->SetBranchAddress("sigMT",     &ZG_sigMT);
  ZGtree->SetBranchAddress("sigMET",    &ZG_sigMET);
  ZGtree->SetBranchAddress("sigMETPhi", &ZG_sigMETPhi);
  ZGtree->SetBranchAddress("dPhiLepMET",&ZG_dPhiLepMET);
  ZGtree->SetBranchAddress("nVertex",   &ZG_nVertex);
//  ZGtree->SetBranchAddress("dRPhoLep",  &ZG_dRPhoLep);
//  ZGtree->SetBranchAddress("dRPhoTrail",&ZG_dRPhoTrail);
  ZGtree->SetBranchAddress("HT",        &ZG_HT);
  ZGtree->SetBranchAddress("nJet",      &ZG_nJet);
//	ZGtree->SetBranchAddress("threeMass", &ZG_threeMass);
	ZGtree->SetBranchAddress("dilepMass", &ZG_dilepMass);
  ZGtree->SetBranchAddress("mcPID",    &ZG_mcPID);
  ZGtree->SetBranchAddress("mcEta",    &ZG_mcEta);
  ZGtree->SetBranchAddress("mcPhi",    &ZG_mcPhi);
  ZGtree->SetBranchAddress("mcPt",     &ZG_mcPt);
  ZGtree->SetBranchAddress("mcMomPID", &ZG_mcMomPID);
  ZGtree->SetBranchAddress("mcGMomPID",&ZG_mcGMomPID);
	ZGtree->SetBranchAddress("mcStatus", &ZG_mcStatus);

  for(unsigned ievt(0); ievt<ZGtree->GetEntries(); ++ievt){//loop on entries
		ZGtree->GetEntry(ievt);

		if( fabs(ZG_nJet-0) > 0.25 )continue;
		if(ZG_dilepMass < 80 || ZG_dilepMass > 100)continue;
		double weight = 35.8*1000*4895.0/48166771; 

	//	bool   isTruePho(false);
	//	double mindR(0.3);
	//	unsigned phoIndex(0);
	//	for(unsigned iMC(0); iMC< ZG_mcPID->size(); iMC++){
	//		double dR1 = DeltaR((*ZG_mcEta)[iMC], (*ZG_mcPhi)[iMC], ZG_phoEta, ZG_phoPhi);
	//		double dE1 = fabs((*ZG_mcPt)[iMC] - ZG_phoEt)/ZG_phoEt;
	//		if(dR1 < mindR && dE1 < 0.2){mindR=dR1; phoIndex=iMC;}
	//	}
	//	if(mindR < 0.1){
	//		if((*ZG_mcPID)[phoIndex] == 22 && (fabs((*ZG_mcMomPID)[phoIndex]) == 23 || fabs((*ZG_mcMomPID)[phoIndex]) == 13 || fabs((*ZG_mcMomPID)[phoIndex])==999)){
	//			isTruePho=true;
	//			for(int i(0); i < 7; i++)if(( ((*ZG_mcStatus)[i] >> i)&1) == 1)p_status->Fill(i);
	//		}
	//	}

		p_iso_ZG->Fill(ZG_phoNeuIso, weight);	
	
		p_phoEt_ZG->Fill(ZG_phoEt, weight); 
		p_phoEta_ZG->Fill(ZG_phoEta, weight);
		p_phoPhi_ZG->Fill(ZG_phoPhi, weight);
		p_lepPt_ZG->Fill(ZG_lepPt, weight);
		p_lepEta_ZG->Fill(ZG_lepEta, weight);
		p_lepPhi_ZG->Fill(ZG_lepPhi, weight);
		p_trailPt_ZG->Fill(ZG_trailPt, weight);
		p_sigMET_ZG->Fill(ZG_sigMET, weight);
		p_sigMT_ZG->Fill(ZG_sigMT, weight);
		p_sigMETPhi_ZG->Fill(ZG_sigMETPhi, weight);
		p_dPhiLepMET_ZG->Fill(ZG_dPhiLepMET, weight);
		p_nVertex_ZG->Fill(ZG_nVertex, weight);
		p_dRPhoLep_ZG->Fill(ZG_dRPhoLep, weight);
		p_HT_ZG->Fill(ZG_HT, weight);
		p_nJet_ZG->Fill(ZG_nJet, weight);
		p_invmass_ZG->Fill(ZG_threeMass, weight);
		p_llmass_ZG->Fill(ZG_dilepMass, weight);
	}//loop on  events

	TCanvas *can_phoEt     = new TCanvas("can_phoEt",       "can_phoEt", 600,600); 
	TCanvas *can_phoEta    = new TCanvas("can_phoEta",      "can_phoEta",600,600); 
	TCanvas *can_phoPhi    = new TCanvas("can_phoPhi",      "can_phoPhi",600,600); 
	TCanvas *can_lepPt     = new TCanvas("can_lepPt",       "can_lepPt",600,600); 
	TCanvas *can_lepEta    = new TCanvas("can_lepEta",      "can_lepEta",600,600); 
	TCanvas *can_lepPhi    = new TCanvas("can_lepPhi",      "can_lepPhi",600,600); 
	TCanvas *can_sigMET    = new TCanvas("can_sigMET",      "can_sigMET",600,600); 
	TCanvas *can_sigMT     = new TCanvas("can_sigMT",       "can_sigMT",600,600); 
	TCanvas *can_sigMETPhi = new TCanvas("can_sigMETPhi",   "can_sigMETPhi",600,600); 
	TCanvas *can_dPhiLepMET= new TCanvas("can_dPhiLepMET",  "can_dPhiLepMET",600,600); 
	TCanvas *can_nVertex   = new TCanvas("can_nVertex",     "can_nVertex",600,600); 
	TCanvas *can_dRPhoLep  = new TCanvas("can_dRPhoLep",    "can_dRPhoLep",600,600); 
	TCanvas *can_HT        = new TCanvas("can_HT",          "can_HT",600,600); 
	TCanvas *can_nJet      = new TCanvas("can_nJet",        "can_nJet",600,600); 
	TCanvas *can_invmass   = new TCanvas("can_invmass",     "can_invmass", 600,600);    
	TCanvas *can_llmass   = new TCanvas("can_llmass",     "can_llmass", 600,600);   
	TCanvas *can_scalefactor = new TCanvas("can_scale",   "can_scale", 600,600); 

	//float scalefactor = p_llmass_data->Integral(50,70)/p_llmass_ZG->Integral(50,70);
	float scalefactor = p_nJet_data->Integral(1,2)/p_nJet_ZG->Integral(1,2);
	can_phoEt->cd();
	TPad *phoEt_pad1 = new TPad("phoEt_pad1", "phoEt_pad1", 0, 0.3, 1, 1.0);
	phoEt_pad1->SetBottomMargin(0.1);
	phoEt_pad1->Draw();  
	phoEt_pad1->cd();  
	phoEt_pad1->SetLogy();
	p_phoEt_data->Draw();
	p_phoEt_ZG->SetLineColor(kRed);
	p_phoEt_ZG->Scale(scalefactor);
	p_phoEt_ZG->Draw("same");
	TLegend *leg=new TLegend(0.6,0.9,0.75,0.9);
	leg->AddEntry(p_phoEt_data,"observed");
	leg->AddEntry(p_phoEt_ZG,"ZG");
	leg->Draw("same");

	can_phoEt->cd();
	TPad *phoEt_pad2 = new TPad("phoEt_pad2", "phoEt_pad2", 0, 0.05, 1, 0.25);
	phoEt_pad2->Draw();
	phoEt_pad2->cd();
	TLine *flatratio_eventcount = new TLine(0,1,600,1);
	TH1D *ratio_phoEt = (TH1D*)p_phoEt_ZG->Clone("ratio_phoEt");
	ratio_phoEt->Divide(p_phoEt_data);
	ratio_phoEt->Draw();
	flatratio_eventcount->Draw("same");
	for(int ibin(1); ibin < ratio_phoEt->GetSize(); ibin++){
		std::cout << ibin << " " << ratio_phoEt->GetBinContent(ibin) << std::endl;
	}
	

	can_phoEta->cd();
	p_phoEta_data->Draw();
	p_phoEta_ZG->SetLineColor(kRed);
	p_phoEta_ZG->Scale(scalefactor);
	p_phoEta_ZG->Draw("same");

	can_phoPhi->cd();
	p_phoPhi_data->Draw();
	p_phoPhi_ZG->SetLineColor(kRed);
	p_phoPhi_ZG->Scale(scalefactor);
	p_phoPhi_ZG->Draw("same");

	can_lepPt->cd();
	gPad->SetLogy();
	p_lepPt_data->Draw();
	p_lepPt_ZG->SetLineColor(kRed);
	p_lepPt_ZG->Scale(scalefactor);
	p_lepPt_ZG->Draw("same");

	can_lepEta->cd();
	p_lepEta_data->Draw();
	p_lepEta_ZG->SetLineColor(kRed);
	p_lepEta_ZG->Scale(scalefactor);
	p_lepEta_ZG->Draw("same");
	for(unsigned ibin(1); ibin < 60; ibin++){
		if(p_lepEta_data->GetBinContent(ibin) > 0)std::cout << p_lepEta_ZG->GetBinCenter(ibin) << " " << p_lepEta_ZG->GetBinContent(ibin)/p_lepEta_data->GetBinContent(ibin) << std::endl;
	} 

	can_lepPhi->cd();
	p_lepPhi_data->Draw();
	p_lepPhi_ZG->SetLineColor(kRed);
	p_lepPhi_ZG->Scale(scalefactor);
	p_lepPhi_ZG->Draw("same");

	can_sigMET->cd();
	p_sigMET_data->Draw();
	p_sigMET_ZG->SetLineColor(kRed);
	p_sigMET_ZG->Scale(scalefactor);
	p_sigMET_ZG->Draw("same");

	can_sigMT->cd();
	p_sigMT_data->Draw();
	p_sigMT_ZG->SetLineColor(kRed);
	p_sigMT_ZG->Scale(scalefactor);
	p_sigMT_ZG->Draw("same");

	can_sigMETPhi->cd();
	p_sigMETPhi_data->Draw();
	p_sigMETPhi_ZG->SetLineColor(kRed);
	p_sigMETPhi_ZG->Scale(scalefactor);
	p_sigMETPhi_ZG->Draw("same");

	can_dPhiLepMET->cd();
	p_dPhiLepMET_data->Draw();
	p_dPhiLepMET_ZG->SetLineColor(kRed);
	p_dPhiLepMET_ZG->Scale(scalefactor);
	p_dPhiLepMET_ZG->Draw("same");

	can_nVertex->cd();
	p_nVertex_data->Draw();
	p_nVertex_ZG->SetLineColor(kRed);
	p_nVertex_ZG->Scale(scalefactor);
	p_nVertex_ZG->Draw("same");

	can_dRPhoLep->cd();
	p_dRPhoLep_data->Draw();
	p_dRPhoLep_ZG->SetLineColor(kRed);
	p_dRPhoLep_ZG->Scale(scalefactor);
	p_dRPhoLep_ZG->Draw("same");

	can_HT->cd();
	gPad->SetLogy();
	p_HT_data->Draw();
	p_HT_ZG->SetLineColor(kRed);
	p_HT_ZG->Scale(scalefactor);
	p_HT_ZG->Draw("same");

	can_nJet->cd();
	p_nJet_data->Draw();
	p_nJet_ZG->SetLineColor(kRed);
	p_nJet_ZG->Scale(scalefactor);
	p_nJet_ZG->Draw("same");

	can_invmass->cd();
	p_invmass_data->Draw();
	p_invmass_ZG->SetLineColor(kRed);
	p_invmass_ZG->Scale(scalefactor);
	p_invmass_ZG->Draw("same");

	can_llmass->cd();
	p_llmass_data->Draw();
	p_llmass_ZG->SetLineColor(kRed);
	p_llmass_ZG->Scale(scalefactor);
	p_llmass_ZG->Draw("same");

	can_scalefactor->cd();
	p_status->Draw();
}


