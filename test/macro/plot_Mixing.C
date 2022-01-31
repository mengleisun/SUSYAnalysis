// input is VGamma/src/analysis_Mixing.C
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
#include "TLatex.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TProfile2D.h"
int RunYear = 2017;

float DeltaR(float eta1,float phi1,float eta2,float phi2)
{
	float deltaPhi = TMath::Abs(phi1-phi2);
	float deltaEta = eta1-eta2;
	if(deltaPhi > TMath::Pi())
	deltaPhi = TMath::TwoPi() - deltaPhi;
		return TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
}

void plot_Mixing(){//main 
// for pt<50, pt>50 pt > 130 samples
  double scalefactor1,scalefactor2,scalefactor3;
  if(RunYear==2016){
  	scalefactor1 = 36.47*1000*412.7/18193471;
	scalefactor2 = 36.47*1000*19.75/2989420;
	scalefactor3 = 36.47*1000*0.8099/2834237;
   }
  if(RunYear==2017){
  	scalefactor1 = 27.13*1000*412.7/10283062;
	scalefactor2 = 27.13*1000*19.75/3598774;
	scalefactor3 = 27.13*1000*0.8099/3639621;
   }
   if(RunYear==2018){
  	scalefactor1 = 59.96*1000*412.7/9850083;
	scalefactor2 = 59.96*1000*19.75/4764595;
	scalefactor3 = 59.96*1000*0.8099/4708995;
   }
    gROOT->SetBatch(kTRUE);
  TH1D *mugamma_phoEt_1 = new TH1D("mugamma_phoEt_1","",165,35,200); 
  TH1D *mugamma_phoEt_2 = new TH1D("mugamma_phoEt_2","",165,35,200); 
  TH1D *mugamma_phoEt_3 = new TH1D("mugamma_phoEt_3","",165,35,200); 
	TH1D *egamma_phoEt_1 =  new TH1D("egamma_phoEt_1","", 165,35,200); 
	TH1D *egamma_phoEt_2 =  new TH1D("egamma_phoEt_2","", 165,35,200); 
	TH1D *egamma_phoEt_3 =  new TH1D("egamma_phoEt_3","", 165,35,200); 
  TH1D *mugamma_phoEt_total = new TH1D("mugamma_phoEt_total",";p_{T} (GeV);Events",165,35,200); 
	TH1D *egamma_phoEt_total  = new TH1D("egamma_phoEt_total",";p_{T} (GeV);Events", 165,35,200); 

  mugamma_phoEt_1->Sumw2();
  mugamma_phoEt_2->Sumw2();
  mugamma_phoEt_3->Sumw2();
	egamma_phoEt_1->Sumw2();
	egamma_phoEt_2->Sumw2();
	egamma_phoEt_3->Sumw2();
  mugamma_phoEt_total->Sumw2();
	egamma_phoEt_total->Sumw2();
  // WGToLNuG tree
  TChain *mgtree = new TChain("mgTree");
	mgtree->Add(Form("/eos/uscms/store/user/tmishra/egMC/mixing_WGToLNuG_%d_TH1D.root",RunYear));
  float mg_phoEt=0;
  float mg_phoEta=0;
  float mg_phoPhi=0;

  mgtree->SetBranchAddress("phoEt",    &mg_phoEt);
  mgtree->SetBranchAddress("phoEta",   &mg_phoEta);
  mgtree->SetBranchAddress("phoPhi",   &mg_phoPhi);

	for(unsigned ievt(0); ievt < mgtree->GetEntries(); ievt++){
		mgtree->GetEntry(ievt);
		if(mg_phoEt > 50)continue;
		//if(fabs(mg_phoEta) > 1.4442)continue;
		mugamma_phoEt_1->Fill(mg_phoEt, scalefactor1);	
		mugamma_phoEt_total->Fill(mg_phoEt, scalefactor1);	
	}

  // WGJets_PtG-40-130 tree
  TChain *mg40tree = new TChain("mgTree");
	mg40tree->Add(Form("/eos/uscms/store/user/tmishra/egMC/mixing_WGJet40_%d_TH1D.root",RunYear));
  float mg40_phoEt=0;
  float mg40_phoEta=0;
  float mg40_phoPhi=0;

  mg40tree->SetBranchAddress("phoEt",    &mg40_phoEt);
  mg40tree->SetBranchAddress("phoEta",   &mg40_phoEta);
  mg40tree->SetBranchAddress("phoPhi",   &mg40_phoPhi);

	for(unsigned ievt(0); ievt < mg40tree->GetEntries(); ievt++){
		mg40tree->GetEntry(ievt);
		if(mg40_phoEt <= 50)continue;
		//if(fabs(mg40_phoEta) > 1.4442)continue;
		mugamma_phoEt_2->Fill(mg40_phoEt, scalefactor2);			
		mugamma_phoEt_total->Fill(mg40_phoEt, scalefactor2);
	}

  // WGJets_PtG-130 tree
  TChain *mg130tree = new TChain("mgTree");
	mg130tree->Add(Form("/eos/uscms/store/user/tmishra/egMC/mixing_WGJet130_%d_TH1D.root",RunYear));
  float mg130_phoEt=0;
  float mg130_phoEta=0;
  float mg130_phoPhi=0;

  mg130tree->SetBranchAddress("phoEt",    &mg130_phoEt);
  mg130tree->SetBranchAddress("phoEta",   &mg130_phoEta);
  mg130tree->SetBranchAddress("phoPhi",   &mg130_phoPhi);

	for(unsigned ievt(0); ievt < mg130tree->GetEntries(); ievt++){
		mg130tree->GetEntry(ievt);
		//if(fabs(mg130_phoEta) > 1.4442)continue;
		mugamma_phoEt_3->Fill(mg130_phoEt, scalefactor3);			
		mugamma_phoEt_total->Fill(mg130_phoEt, scalefactor3);
	}

  TChain *egtree = new TChain("egTree");
	egtree->Add(Form("/eos/uscms/store/user/tmishra/egMC/mixing_WGToLNuG_%d_TH1D.root",RunYear));
  float eg_phoEt=0;
  float eg_phoEta=0;
  float eg_phoPhi=0;

  egtree->SetBranchAddress("phoEt",    &eg_phoEt);
  egtree->SetBranchAddress("phoEta",   &eg_phoEta);
  egtree->SetBranchAddress("phoPhi",   &eg_phoPhi);
  cout<<egtree->GetEntries()<<endl;

	for(unsigned ievt(0); ievt < egtree->GetEntries(); ievt++){
		egtree->GetEntry(ievt);
		if(eg_phoEt > 50)continue;
		//if(fabs(eg_phoEta) > 1.4442)continue;
		egamma_phoEt_1->Fill(eg_phoEt, scalefactor1);		
		egamma_phoEt_total->Fill(eg_phoEt, scalefactor1);	
	}


  TChain *eg40tree = new TChain("egTree");
	eg40tree->Add(Form("/eos/uscms/store/user/tmishra/egMC/mixing_WGJet40_%d_TH1D.root",RunYear));
  float eg40_phoEt=0;
  float eg40_phoEta=0;
  float eg40_phoPhi=0;

  eg40tree->SetBranchAddress("phoEt",    &eg40_phoEt);
  eg40tree->SetBranchAddress("phoEta",   &eg40_phoEta);
  eg40tree->SetBranchAddress("phoPhi",   &eg40_phoPhi);
  cout<<eg40tree->GetEntries()<<endl;
	for(unsigned ievt(0); ievt < eg40tree->GetEntries(); ievt++){
		eg40tree->GetEntry(ievt);
		if(eg40_phoEt <= 50)continue;
		//if(fabs(eg40_phoEta) > 1.4442)continue;
		egamma_phoEt_2->Fill(eg40_phoEt, scalefactor2);			
		egamma_phoEt_total->Fill(eg40_phoEt, scalefactor2);
	}

  TChain *eg130tree = new TChain("egTree");
	eg130tree->Add(Form("/eos/uscms/store/user/tmishra/egMC/mixing_WGJet130_%d_TH1D.root",RunYear));
  float eg130_phoEt=0;
  float eg130_phoEta=0;
  float eg130_phoPhi=0;

  eg130tree->SetBranchAddress("phoEt",    &eg130_phoEt);
  eg130tree->SetBranchAddress("phoEta",   &eg130_phoEta);
  eg130tree->SetBranchAddress("phoPhi",   &eg130_phoPhi);
  cout<<eg130tree->GetEntries()<<endl;

	for(unsigned ievt(0); ievt < eg130tree->GetEntries(); ievt++){
		eg130tree->GetEntry(ievt);
		//if(fabs(eg130_phoEta) > 1.4442)continue;
		egamma_phoEt_3->Fill(eg130_phoEt, scalefactor3);			
		egamma_phoEt_total->Fill(eg130_phoEt, scalefactor3);
	}

  double ratio1 = egamma_phoEt_1->Integral(6,8)/egamma_phoEt_2->Integral(6,8);
  double ratio2 = mugamma_phoEt_1->Integral(6,8)/mugamma_phoEt_2->Integral(6,8);

  double ratio1_130 = egamma_phoEt_1->Integral(15,30)/egamma_phoEt_3->Integral(15,30);
  double ratio2_130 = mugamma_phoEt_1->Integral(15,30)/mugamma_phoEt_3->Integral(15,30);

	gStyle->SetOptStat(0);
	TCanvas *can1 = new TCanvas("can1","",600,600);
  can1->cd();
	gPad->SetLogy();
	egamma_phoEt_total->SetFillColor(17);
	egamma_phoEt_total->SetLineColor(17);
	egamma_phoEt_total->SetFillStyle(1001);
	egamma_phoEt_total->Draw("hist");
	egamma_phoEt_1->SetMarkerStyle(7);
	egamma_phoEt_1->SetMarkerColor(kBlue);
	egamma_phoEt_1->SetLineColor(kBlue);
	egamma_phoEt_1->Draw("EP same");
	egamma_phoEt_2->SetMarkerStyle(7);
	egamma_phoEt_2->SetMarkerColor(kRed);
	egamma_phoEt_2->SetLineColor(kRed);
	egamma_phoEt_2->Draw("EP same");
	egamma_phoEt_3->SetMarkerStyle(7);
	egamma_phoEt_3->SetMarkerColor(kGreen);
	egamma_phoEt_3->SetLineColor(kGreen);
	egamma_phoEt_3->Draw("EP same");
	TLegend *leg =  new TLegend(0.4,0.6,0.87,0.9);
	leg->SetTextSize(0.025);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->AddEntry(egamma_phoEt_total, "mixed W#gamma");
	leg->AddEntry(egamma_phoEt_1,"WGToLNuG");
	leg->AddEntry(egamma_phoEt_2,"WGJets_MonoPhoton_PtG-40to130");
	leg->AddEntry(egamma_phoEt_3,"WGJets_MonoPhoton_PtG-130");
	leg->Draw("same");
	can1->SaveAs(Form("/eos/uscms/store/user/tmishra/VGamma/WGMixing_%d.pdf",RunYear));
}
