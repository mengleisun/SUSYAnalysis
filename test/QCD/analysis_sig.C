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

#include "../../include/analysis_rawData.h"
#include "../../include/analysis_photon.h"
#include "../../include/analysis_muon.h"
#include "../../include/analysis_ele.h"
#include "../../include/analysis_mcData.h"
#include "../../include/analysis_tools.h"
#include "../../include/analysis_fakes.h"

void analysis_sig(int ichannel, double inputmetlow, double inputmetup, int leplow, int lephigh){//main  

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");
  int channelType = ichannel; // eg = 1; mg =2;
	double met_lowcut = inputmetlow;
	double met_upcut  = inputmetup;
	//*********** histo list **********************//
	std::ostringstream histname;
	Double_t plotEtBins[]={35,40,45,50,55,60,65,70,75,80, 85,90,95,100,105,110,115,120,125,130, 135,140,146,152,158,164,170,177,184,192, 200,208,216,224,232,240,250,260,275,290, 305,325,345,370,400};
	TH1F *p_PhoEt = new TH1F("p_PhoEt","#gamma E_{T}; E_{T} (GeV)",44,plotEtBins);
	TH1F *p_PhoEta = new TH1F("p_PhoEta","#gamma #eta; #eta;",60,-3,3);
	Double_t plotPtBins[]={25,30,35,40,45,50,55,60,65,70,75,80, 85,90,95,100,105,110,115,120,125,130, 135,140,146,152,158,164,170,177,184,192, 200,208,216,224,232,240,250,260,275,290, 305,325,345,370,400};
	TH1F *p_LepPt = new TH1F("p_LepPt","p_LepPt",46,plotPtBins);
	TH1F *p_LepEta = new TH1F("p_LepEta","p_LepEta",60,-3,3);
	TH1F *p_MET = new TH1F("p_MET","MET; MET (GeV);",100,0,1000);
	TH1F *p_Mt = new TH1F("p_Mt","M_{T}; M_{T} (GeV);",200,0,1000);
	TH1F *p_HT = new TH1F("p_HT","HT; HT (GeV);",100,0,1200); 
	TH1F *p_dPhiEleMET = new TH1F("p_dPhiEleMET","dPhiEleMET",32,0,3.2); 
	TH1F *p_PU = new TH1F("p_PU","",100,0,100);
	//************ Signal Tree **********************//
	TChain *sigtree = new TChain("signalTree");
	if(channelType==1)sigtree->Add("/uscms_data/d3/mengleis/test_WG.root");
	if(channelType==2)sigtree->Add("/uscms_data/d3/mengleis/FullStatusData/resTree_mgsignal_MuonEG_FebReminiAOD.root");

	float phoEt(0);
	float phoEta(0);
	float phoPhi(0);
	float lepPt(0);
	float lepEta(0);
	float lepPhi(0);
	float sigMT(0);
	float sigMET(0);
	float dPhiLepMET(0);
	float HT(0);
	int   nVertex(0);
	float dRPhoLep(0);	
	sigtree->SetBranchAddress("phoEt",     &phoEt);
	sigtree->SetBranchAddress("phoEta",    &phoEta);
	sigtree->SetBranchAddress("phoPhi",    &phoPhi);
	sigtree->SetBranchAddress("lepPt",     &lepPt);
	sigtree->SetBranchAddress("lepEta",    &lepEta);
	sigtree->SetBranchAddress("lepPhi",    &lepPhi);
	sigtree->SetBranchAddress("sigMT",     &sigMT);
	sigtree->SetBranchAddress("sigMET",    &sigMET);
  sigtree->SetBranchAddress("HT",        &HT);
	sigtree->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
	sigtree->SetBranchAddress("nVertex",   &nVertex);
	sigtree->SetBranchAddress("dRPhoLep",  &dRPhoLep);

	for (unsigned ievt(0); ievt<sigtree->GetEntries(); ++ievt){//loop on entries
		sigtree->GetEntry(ievt);
		p_PU->Fill(nVertex);
		if(sigMET > met_upcut || sigMET < met_lowcut)continue;
		if(phoEt < 35 || lepPt < 25)continue;
		if(fabs(phoEta) > 1.4442 || fabs(lepEta) > 2.5)continue;

		if(lepPt < leplow || lepPt > lephigh)continue;
		//if(fabs(lepEta) > 1.4442)continue;
		//if(phoEt < 80)continue;
		//if(sigMT > 100)continue;

		p_PhoEt->Fill(phoEt);
		p_PhoEta->Fill(phoEta);
		p_LepPt->Fill(lepPt);
		p_LepEta->Fill(lepEta);
		p_MET->Fill(sigMET);
		p_Mt->Fill(sigMT);
		p_HT->Fill(HT);
		p_dPhiEleMET->Fill(fabs(dPhiLepMET),1);
	}        

	std::ostringstream outputname;
	
	if(channelType==1 && inputmetup < 100)outputname << "WGcontrolTree_eg_signal_ReMiniAOD" << leplow << "_" << lephigh <<".root";
	else if(channelType==1 && inputmetup >= 100)outputname << "bkgTree_eg_signal_ReMiniAOD.root";
	else if(channelType==2 && inputmetup < 100)outputname << "WGcontrolTree_mg_signal_ReMiniAOD"<< leplow << "_" << lephigh <<".root";
	else if(channelType==2 && inputmetup >= 100)outputname << "bkgTree_mg_signal_ReMiniAOD.root";
	TFile *outputfile = TFile::Open(outputname.str().c_str(),"RECREATE");
	outputfile->cd();
	p_PhoEt->Write();
	p_PhoEta->Write();
	p_LepPt->Write();
	p_LepEta->Write();
	p_MET->Write();
	p_Mt->Write();
	p_HT->Write();
	p_dPhiEleMET->Write();
	p_PU->Write();
	outputfile->Write();
	outputfile->Close();
}


