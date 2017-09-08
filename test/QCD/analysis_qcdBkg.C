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

void analysis_qcdBkg(int ichannel, int inputmetlow, int inputmetup,int leplow, int lephigh, int isocut ){//main  
	bool toDeriveScale(true);

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

	double factor_egFakeLepton(0.046);
	double factorerror_egFakeLepton(0.0039);
	double factor_mgFakeLepton(0.1);
	double factorerror_mgFakeLepton(0.05);

  int channelType = ichannel; // eg = 1; mg =2;
	double met_lowcut = inputmetlow;
	double met_upcut  = inputmetup;
	double factorQCD(1);
	double factorQCDUP = factorQCD*(1+0);
	double factorQCDDO = factorQCD*(1-0);
	if(toDeriveScale){
		factorQCD = 1;
		factorQCDUP = 1;
		factorQCDDO = 1;
	}
	else{
		if(channelType == 1){
			factorQCD = factor_egFakeLepton;
			factorQCDUP = factor_egFakeLepton+factorerror_egFakeLepton;
			factorQCDDO = factor_egFakeLepton-factorerror_egFakeLepton;
		}
		else if(channelType == 2){
			factorQCD = factor_mgFakeLepton;
			factorQCDUP = factor_mgFakeLepton + factorerror_mgFakeLepton;
			factorQCDDO = factor_mgFakeLepton - factorerror_mgFakeLepton;
		}
	}

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
	TH1F *p_dPhiEleMET = new TH1F("p_dPhiEleMET","dPhiEleMET",32,0,3.2); 
	TH1F *p_PU = new TH1F("p_PU","",100,0,100);

// ********** fake lepton tree ************** //
  TChain *fakeEtree = new TChain("fakeETree","fakeETree");
	if(channelType==1)fakeEtree->Add("/uscms_data/d3/mengleis/test_WG.root");
	if(channelType==2)fakeEtree->Add("/uscms_data/d3/mengleis/FullStatusData/resTree_mgsignal_MuonEG_FebReminiAOD.root");
  float phoEt(0);
  float phoEta(0);
  float phoPhi(0);
  float lepPt(0);
  float lepEta(0);
  float lepPhi(0);
	float fakeElepMiniIso(0);
  float sigMT(0);
  float sigMET(0);
  float sigMETPhi(0);
  float dPhiLepMET(0);
  int   nVertex(0);
  float dRPhoLep(0);
  float HT(0);
  float nJet(0);
  
  fakeEtree->SetBranchAddress("phoEt",     &phoEt);
  fakeEtree->SetBranchAddress("phoEta",    &phoEta);
  fakeEtree->SetBranchAddress("phoPhi",    &phoPhi);
  fakeEtree->SetBranchAddress("lepPt",     &lepPt);
  fakeEtree->SetBranchAddress("lepEta",    &lepEta);
  fakeEtree->SetBranchAddress("lepPhi",    &lepPhi);
  fakeEtree->SetBranchAddress("fakeElepMiniIso", &fakeElepMiniIso);
  fakeEtree->SetBranchAddress("sigMT",     &sigMT);
  fakeEtree->SetBranchAddress("sigMET",    &sigMET);
  fakeEtree->SetBranchAddress("sigMETPhi", &sigMETPhi);
  fakeEtree->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
  fakeEtree->SetBranchAddress("nVertex",   &nVertex);
  fakeEtree->SetBranchAddress("dRPhoLep",  &dRPhoLep);
  fakeEtree->SetBranchAddress("HT",        &HT);
  fakeEtree->SetBranchAddress("nJet",      &nJet);

	for(unsigned ievt(0); ievt < fakeEtree->GetEntries(); ievt++){
		fakeEtree->GetEntry(ievt);
		p_PU->Fill(nVertex);
		if(sigMET > met_upcut || sigMET < met_lowcut)continue;
		if(phoEt < 35 || lepPt < 25)continue;
		if(fabs(phoEta) > 1.4442 || fabs(lepEta) > 2.5)continue;

		if(lepPt < leplow || lepPt > lephigh)continue;
    if(fakeElepMiniIso > isocut*0.1)continue;
//		if(fabs(lepEta) > 1.4442)continue;
		//if(phoEt < 80)continue;
//		if(sigMT > 100)continue;
	
		p_PhoEt->Fill(phoEt, factorQCD);
		p_PhoEta->Fill(phoEta, factorQCD);
		p_LepPt->Fill(lepPt, factorQCD);
		p_LepEta->Fill(lepEta, factorQCD);
		p_MET->Fill(sigMET, factorQCD);
		p_Mt->Fill(sigMT, factorQCD);
		p_dPhiEleMET->Fill(fabs(dPhiLepMET), 1);
	}

	std::ostringstream outputname;
	
	if(channelType==1 && inputmetup < 100)outputname << "WGcontrolTree_eg_qcd_ReMiniAOD_" << leplow << "_" << lephigh << "_met" << inputmetlow << "_iso" << isocut << ".root";
	else if(channelType==1 && inputmetup >= 100)outputname << "bkgTree_eg_qcd_ReMiniAOD.root";
	else if(channelType==2 && inputmetup < 100)outputname << "WGcontrolTree_mg_qcd_ReMiniAOD_"<< leplow << "_" << lephigh << "_met" << inputmetlow << "_iso" << isocut << ".root";
	else if(channelType==2 && inputmetup >= 100)outputname << "bkgTree_mg_qcd_ReMiniAOD.root";
	TFile *outputfile = TFile::Open(outputname.str().c_str(),"RECREATE");
	outputfile->cd();
	p_PhoEt->Write();
	p_PhoEta->Write();
	p_LepPt->Write();
	p_LepEta->Write();
	p_MET->Write();
	p_Mt->Write();
	p_dPhiEleMET->Write();
	p_PU->Write();
	outputfile->Write();
	outputfile->Close();
}


