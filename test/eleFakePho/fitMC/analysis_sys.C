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
#include "TPad.h"

#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_mcData.h"
#include "../../../include/analysis_tools.h"
#include "../../../include/analysis_fakes.h"

#define NTOY 150

void analysis_sys(){//main  

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");


  int channelType = 1; // eg = 1; mg =2;
	double met_lowcut(0);
	double met_upcut(1000000000);

	TF3 *h_toymc_fakerate[NTOY];
	std::ostringstream funcname;
	std::ifstream elefake_file("ToyFakeRate2.txt");
	double scalefactor(0);
	double ptslope(0);
	double ptconstant(0);
	double ptindex(0);
	double ptcoeff(0);
	double vtxconst(0);
	double vtxslope(0);
	if(elefake_file.is_open()){
  	for(int i(0); i<NTOY; i++){ 
			elefake_file >> scalefactor >> ptslope >> ptconstant >> ptindex >> ptcoeff >> vtxconst >> vtxslope;
			funcname.str("");
			funcname << "h_toymc_fakerate_" << i;
			h_toymc_fakerate[i] = new TF3(funcname.str().c_str(), fakerate_func,10,1000,0,100,0,1.5,7);
			h_toymc_fakerate[i]->SetParameters(scalefactor, ptslope, ptconstant, ptindex, ptcoeff, vtxconst, vtxslope); 
	  }
	}
	elefake_file.close();
  /**********************************/
	/*	double normfactor = par[0]; 	*/  
  /*	double slope = par[1];				*/
  /*	double constant = par[2];			*/
  /*	double index = par[3];				*/
  /*	double coeff = par[4]; 				*/
	/*	double vtx_constant = par[5];	*/
	/*	double vtx_slope = par[6];		*/
  /**********************************/

	//*********** histo list **********************//
		std::ostringstream histname;

		TH1F *p_elebkgPhoEt[1000];
		TH1F *p_elebkgLepPt[1000];
		TH1F *p_elebkgMET[1000];
		TH1F *p_elebkgMt[1000];
		TH1F *p_eledPhiEleMET[1000];
		TProfile *pro_elebkgPhoEt = new TProfile("pro_elebkgPhoEt", "pro_elebkgPhoEt" , 75,0,1500);
		TProfile *pro_elebkgLepPt = new TProfile("pro_elebkgLepPt", "pro_elebkgLepPt" , 1500,0,1500);
		TProfile *pro_elebkgMET   = new TProfile("pro_elebkgMET", "pro_elebkgMET", 50,0,500);
		TProfile *pro_elebkgMt    = new TProfile("pro_elebkgMt", "pro_elebkgMt", 50,0,500);
		TProfile *pro_eledPhiEleMET = new TProfile("pro_eledPhiEleMET","pro_eledPhiEleMET", 64,-3.2, 3.2);
		TH1F *errorscale_elebkgPhoEt	 = new TH1F("errorscale_elebkgPhoEt", "pro_elebkgPhoEt" , 75,0,1500);
		TH1F *errorscale_elebkgLepPt	 = new TH1F("errorscale_elebkgLepPt", "pro_elebkgLepPt" , 1500,0,1500);
		TH1F *errorscale_elebkgMET		 = new TH1F("errorscale_elebkgMET", "pro_elebkgMET", 50,0,500);
		TH1F *errorscale_elebkgMt			 = new TH1F("errorscale_elebkgMt", "pro_elebkgMt", 50,0,500);
		TH1F *errorscale_eledPhiEleMET = new TH1F("errorscale_eledPhiEleMET","pro_eledPhiEleMET", 64,-3.2, 3.2);
		for(unsigned ih(0); ih < 1000; ih++){
			histname.str("");
			histname << "p_eledPhiEleMET_" << ih;
			p_eledPhiEleMET[ih] = new TH1F(histname.str().c_str(), histname.str().c_str(), 128,-3.2, 3.2);
			histname.str("");
      histname << "p_elebkgPhoEt_" << ih;
			p_elebkgPhoEt[ih] = new TH1F(histname.str().c_str(), histname.str().c_str(), 1500,0,1500);
			histname.str("");
      histname << "p_elebkgLepPt_" << ih;
			p_elebkgLepPt[ih] = new TH1F(histname.str().c_str(), histname.str().c_str(), 1500,0,1500);
			histname.str("");
      histname << "p_p_elebkgMET_" << ih;
			p_elebkgMET[ih] = new TH1F(histname.str().c_str(), histname.str().c_str(), 500,0,500);
			histname.str("");
      histname << "p_p_elebkgMt_" << ih;
			p_elebkgMt[ih] = new TH1F(histname.str().c_str(), histname.str().c_str(), 500,0,500);
		}
		pro_elebkgPhoEt->SetErrorOption("s");	
		pro_elebkgLepPt->SetErrorOption("s");	
		pro_elebkgMET->SetErrorOption("s");	
		pro_elebkgMt->SetErrorOption("s");	
		pro_eledPhiEleMET->SetErrorOption("s");	
		pro_elebkgPhoEt->SetMarkerStyle(20);
		pro_elebkgLepPt->SetMarkerStyle(20);
		pro_elebkgMET->SetMarkerStyle(20);
		pro_elebkgMt->SetMarkerStyle(20);
		pro_eledPhiEleMET->SetMarkerStyle(20);
	//************ Proxy Tree **********************//
		TChain *proxytree = new TChain("proxyTree");
		if(channelType==1)proxytree->Add("/uscms_data/d3/mengleis/ReMiniAOD/resTree_egsignal_DoubleEG_ReMiniAOD.root");
		if(channelType==2)proxytree->Add("/uscms_data/d3/mengleis/Rereco/resTree_mgsignal_2016Rereco_MuonEG.root");

		float proxyphoEt(0);
		float proxyphoEta(0);
		float proxyphoPhi(0);
		float proxylepPt(0);
		float proxylepEta(0);
		float proxylepPhi(0);
		float proxysigMT(0);
		float proxysigMET(0);
		float proxysigMETPhi(0);
		float proxydPhiLepMET(0);
		int   proxynVertex(0);
		float proxydRPhoLep(0);
		float proxyHT(0);
		float proxynJet(0);
		
		proxytree->SetBranchAddress("phoEt",     &proxyphoEt);
		proxytree->SetBranchAddress("phoEta",    &proxyphoEta);
		proxytree->SetBranchAddress("phoPhi",    &proxyphoPhi);
		proxytree->SetBranchAddress("lepPt",     &proxylepPt);
		proxytree->SetBranchAddress("lepEta",    &proxylepEta);
		proxytree->SetBranchAddress("lepPhi",    &proxylepPhi);
		proxytree->SetBranchAddress("sigMT",     &proxysigMT);
		proxytree->SetBranchAddress("sigMET",    &proxysigMET);
		proxytree->SetBranchAddress("sigMETPhi", &proxysigMETPhi);
		proxytree->SetBranchAddress("dPhiLepMET",&proxydPhiLepMET);
		proxytree->SetBranchAddress("nVertex",   &proxynVertex);
		proxytree->SetBranchAddress("dRPhoLep",  &proxydRPhoLep);
		proxytree->SetBranchAddress("HT",        &proxyHT);
		proxytree->SetBranchAddress("nJet",      &proxynJet);
	 
		for (unsigned ievt(0); ievt<proxytree->GetEntries(); ++ievt){//loop on entries
			proxytree->GetEntry(ievt);
			if(proxyphoEt <= 35)continue;
			if(proxysigMET > met_upcut || proxysigMET < met_lowcut)continue;
			for(unsigned it(0); it < NTOY; it++){
				double w_ele = h_toymc_fakerate[it]->Eval(proxyphoEt, proxynVertex, fabs(proxyphoEta));
				p_elebkgPhoEt[it]->Fill(proxyphoEt,w_ele);
				p_elebkgMET[it]->Fill(proxysigMET, w_ele);
				p_elebkgMt[it]->Fill(proxysigMT, w_ele);
				p_elebkgLepPt[it]->Fill(proxylepPt, w_ele);

				float deltaPhi = proxylepPhi - proxysigMETPhi;
				if(fabs(deltaPhi) > TMath::Pi()){
					if(deltaPhi > 0)deltaPhi = -1.0*(TMath::TwoPi() - fabs(deltaPhi));
					else deltaPhi = TMath::TwoPi() - fabs(deltaPhi);
				}
				p_eledPhiEleMET[it]->Fill(deltaPhi, w_ele);
			}
		}

	gStyle->SetOptStat(0);
	Double_t plotPtBins[]={35,40,45,50,55,60,65,70,75,80,85,90,95,100,108,116,124,132,140,148,156,164,172,180,190,200,210,220,240,260,300};
	TCanvas *c_pt = new TCanvas("Photon_Pt", "Photon P_{T}",1200,600);
	c_pt->Divide(2);
	c_pt->cd(1);
	gPad->SetLogy();
  for(unsigned i(0); i<NTOY; i++){
		p_elebkgPhoEt[i]->Rebin(20);	
		p_elebkgPhoEt[i]->SetLineColor(kRed);
		if(i==0)p_elebkgPhoEt[i]->Draw("hist");
		else p_elebkgPhoEt[i]->Draw("hist same");
	}
	for(unsigned it(0); it < NTOY; it++){
	  for(int ibin(1); ibin < p_elebkgPhoEt[0]->GetSize(); ibin++)	
			pro_elebkgPhoEt->Fill(p_elebkgPhoEt[it]->GetBinCenter(ibin), p_elebkgPhoEt[it]->GetBinContent(ibin));
	}
	pro_elebkgPhoEt->SetLineColor(kBlue);
	pro_elebkgPhoEt->SetMarkerColor(kBlue);
	pro_elebkgPhoEt->Draw("same");
	for(int ibin(1); ibin < pro_elebkgPhoEt->GetSize(); ibin++)	
		if(pro_elebkgPhoEt->GetBinContent(ibin)> 0){
			errorscale_elebkgPhoEt->Fill(pro_elebkgPhoEt->GetBinCenter(ibin), pro_elebkgPhoEt->GetBinError(ibin)/pro_elebkgPhoEt->GetBinContent(ibin));
			errorscale_elebkgPhoEt->SetBinError(ibin, 0);
		}
	c_pt->cd(2);
	errorscale_elebkgPhoEt->Draw();

	TCanvas *c_met = new TCanvas("MET", "MET",1200,600);
	c_met->Divide(2);
	c_met->cd(1);
	gPad->SetLogy();
  for(unsigned i(0); i<NTOY; i++){
		p_elebkgMET[i]->Rebin(10);	
		p_elebkgMET[i]->SetLineColor(kRed);
		if(i==0)p_elebkgMET[i]->Draw("hist");
		else p_elebkgMET[i]->Draw("hist same");
	}
	for(unsigned it(0); it < NTOY; it++){
	  for(int ibin(1); ibin < p_elebkgMET[0]->GetSize(); ibin++)	
		  pro_elebkgMET->Fill(p_elebkgMET[it]->GetBinCenter(ibin), p_elebkgMET[it]->GetBinContent(ibin));
	}		
	pro_elebkgMET->SetLineColor(kBlue);
	pro_elebkgMET->SetMarkerColor(kBlue);
	pro_elebkgMET->Draw("same");
	for(int ibin(1); ibin < pro_elebkgMET->GetSize(); ibin++)	
		if(pro_elebkgMET->GetBinContent(ibin)> 0){
			errorscale_elebkgMET->Fill(pro_elebkgMET->GetBinCenter(ibin), pro_elebkgMET->GetBinError(ibin)/pro_elebkgMET->GetBinContent(ibin));
			errorscale_elebkgMET->SetBinError(ibin, 0);
		}
	c_met->cd(2);
	errorscale_elebkgMET->Draw();


	TCanvas *c_mt = new TCanvas("Photon_Mt", "MT",600,600);
	c_mt->cd();
	gPad->SetLogy();
  for(unsigned i(0); i<NTOY; i++){
		p_elebkgMt[i]->Rebin(10);	
		p_elebkgMt[i]->SetLineColor(kRed);
		if(i==0)p_elebkgMt[i]->Draw("hist");
		else p_elebkgMt[i]->Draw("hist same");
	}
	for(unsigned it(0); it < NTOY; it++){
	  for(int ibin(1); ibin < p_elebkgMt[0]->GetSize(); ibin++)	
			pro_elebkgMt->Fill(p_elebkgMt[it]->GetBinCenter(ibin), p_elebkgMt[it]->GetBinContent(ibin));
	}		


	TCanvas *c_dphi = new TCanvas("Photon_dPhiEleMET", "MT",1200,600);
	c_dphi->Divide(2);
	c_dphi->cd(1);
  for(unsigned i(0); i<NTOY; i++){
		p_eledPhiEleMET[i]->Rebin(2);	
		p_eledPhiEleMET[i]->SetLineColor(kRed);
		if(i==0)p_eledPhiEleMET[i]->Draw("hist");
		else p_eledPhiEleMET[i]->Draw("hist same");
	}
	for(unsigned it(0); it < NTOY; it++){
	  for(int ibin(1); ibin < p_eledPhiEleMET[0]->GetSize(); ibin++)	
			pro_eledPhiEleMET->Fill(p_eledPhiEleMET[it]->GetBinCenter(ibin), p_eledPhiEleMET[it]->GetBinContent(ibin));
	}		
	pro_eledPhiEleMET->SetLineColor(kBlue);
	pro_eledPhiEleMET->SetMarkerColor(kBlue);
	pro_eledPhiEleMET->Draw("same");
	for(int ibin(1); ibin < pro_eledPhiEleMET->GetSize(); ibin++)	
		if(pro_eledPhiEleMET->GetBinContent(ibin)> 0){
			errorscale_eledPhiEleMET->Fill(pro_eledPhiEleMET->GetBinCenter(ibin), pro_eledPhiEleMET->GetBinError(ibin)/pro_eledPhiEleMET->GetBinContent(ibin));
			errorscale_eledPhiEleMET->SetBinError(ibin, 0);
		}
	c_dphi->cd(2);
	errorscale_eledPhiEleMET->Draw();

}


