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
#include "../../include/analysis_scalefactor.h"

void analysis_VGBkg(int ichannel, double inputmetlow, double inputmetup,int leplow, int lephigh ){//main  
	
	bool toDeriveScale(true);

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");
	
	esfScaleFactor  objectESF;
	double ZGammaFactor(1.90);
	double WGammaFactor(1.36);
	double factor_egVGamma(1.14);
	double factorerror_egVGamma(0.07);
	double factor_mgVGamma(0.83);
	double factorerror_mgVGamma(0.11);


  int channelType = ichannel; // eg = 1; mg =2;
	double met_lowcut = inputmetlow;
	double met_upcut  = inputmetup;
	double factorMCE(1);
	double factorMCEUP = factorMCE*(1+0);
	double factorMCEDO = factorMCE*(1-0);
	if(toDeriveScale){
		factorMCE = 1;
		factorMCEUP = 1;
		factorMCEDO = 1;
	}
	else{
		if(channelType == 1){
			factorMCE = factor_egVGamma;
			factorMCEUP = factor_egVGamma+factorerror_egVGamma;
			factorMCEDO = factor_egVGamma-factorerror_egVGamma;
		}
		else if(channelType == 2){
			factorMCE = factor_mgVGamma;
			factorMCEUP = factor_mgVGamma + factorerror_mgVGamma;
			factorMCEDO = factor_mgVGamma - factorerror_mgVGamma;
		}
	}
	//*********** histo list **********************//
	std::ostringstream outputname;
	if(channelType==1 && inputmetup < 100)outputname << "controlTree_eg_VGBkg_ReMiniAOD" << leplow << "_" << lephigh <<".root";
	else if(channelType==1 && inputmetup >= 100)outputname << "bkgTree_eg_VGBkg_ReMiniAOD.root";
	else if(channelType==2 && inputmetup < 100)outputname << "controlTree_mg_VGBkg_ReMiniAOD"<< leplow << "_" << lephigh <<".root";
	else if(channelType==2 && inputmetup >= 100)outputname << "bkgTree_mg_VGBkg_ReMiniAOD.root";
	TFile *outputfile = TFile::Open(outputname.str().c_str(),"RECREATE");
	outputfile->cd();
	std::ostringstream histname;

	Double_t plotEtBins[]={35,40,45,50,55,60,65,70,75,80, 85,90,95,100,105,110,115,120,125,130, 135,140,146,152,158,164,170,177,184,192, 200,208,216,224,232,240,250,260,275,290, 305,325,345,370,400};
	TH1F *p_PhoEt = new TH1F("p_PhoEt","#gamma E_{T}; E_{T} (GeV)",44,plotEtBins);
	TH1F *p_PhoEt_low_WG = new TH1F("p_PhoEt_low_WG","#gamma E_{T}; E_{T} (GeV)",1000,0,1000);
	TH1F *p_PhoEt_high_WG = new TH1F("p_PhoEt_high_WG","#gamma E_{T}; E_{T} (GeV)",1000,0,1000);
	TH1F *p_PhoEt_low_ZG = new TH1F("p_PhoEt_low_ZG","#gamma E_{T}; E_{T} (GeV)",1000,0,1000);
	TH1F *p_PhoEt_high_ZG = new TH1F("p_PhoEt_high_ZG","#gamma E_{T}; E_{T} (GeV)",1000,0,1000);
	TH1F *p_PhoEta = new TH1F("p_PhoEta","#gamma #eta; #eta;",60,-3,3);
	Double_t plotPtBins[]={25,30,35,40,45,50,55,60,65,70,75,80, 85,90,95,100,105,110,115,120,125,130, 135,140,146,152,158,164,170,177,184,192, 200,208,216,224,232,240,250,260,275,290, 305,325,345,370,400};
	TH1F *p_LepPt = new TH1F("p_LepPt","p_LepPt",46,plotPtBins);
	TH1F *p_LepEta = new TH1F("p_LepEta","p_LepEta",60,-3,3);
	TH1F *p_HT = new TH1F("p_HT","HT; HT (GeV);",100,0,1200); 
	TH1F *p_MET = new TH1F("p_MET","MET; MET (GeV);",100,0,1000);
	TH1F *p_Mt = new TH1F("p_Mt","M_{T}; M_{T} (GeV);",200,0,1000); 
	TH1F *p_dPhiEleMET = new TH1F("p_dPhiEleMET","dPhiEleMET",32,0,3.2); 
	TH1F *p_PU = new TH1F("p_PU","",100,0,100);

	TH1F *jesup_MET = new TH1F("jesup_MET","MET; MET (GeV);",100,0,1000);
	TH1F *jesup_Mt = new TH1F("jesup_Mt","M_{T}; M_{T} (GeV);",200,0,1000);
	TH1F *jesup_HT = new TH1F("jesup_HT","HT; HT (GeV);",100,0,1200); 
	TH1F *jesup_dPhiEleMET = new TH1F("jesup_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	TH1F *jesdo_MET = new TH1F("jesdo_MET","MET; MET (GeV);",100,0,1000);
	TH1F *jesdo_Mt = new TH1F("jesdo_Mt","M_{T}; M_{T} (GeV);",200,0,1000); 
	TH1F *jesdo_HT = new TH1F("jesdo_HT","HT; HT (GeV);",100,0,1200); 
	TH1F *jesdo_dPhiEleMET = new TH1F("jesdo_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	TH1F *jerup_MET = new TH1F("jerup_MET","MET; MET (GeV);",100,0,1000);
	TH1F *jerup_Mt = new TH1F("jerup_Mt","M_{T}; M_{T} (GeV);",200,0,1000); 
	TH1F *jerup_dPhiEleMET = new TH1F("jerup_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	TH1F *jerdo_MET = new TH1F("jerdo_MET","MET; MET (GeV);",100,0,1000);
	TH1F *jerdo_Mt = new TH1F("jerdo_Mt","M_{T}; M_{T} (GeV);",200,0,1000); 
	TH1F *jerdo_dPhiEleMET = new TH1F("jerdo_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	TH1F *scaleup_PhoEt = new TH1F("scaleup_PhoEt","#gamma E_{T}; E_{T} (GeV)",44,plotEtBins);
	TH1F *scaleup_PhoEta = new TH1F("scaleup_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1F *scaleup_LepPt = new TH1F("scaleup_LepPt","scaleup_LepPt",46,plotPtBins);
	TH1F *scaleup_LepEta = new TH1F("scaleup_LepEta","scaleup_LepEta",60,-3,3);
	TH1F *scaleup_MET = new TH1F("scaleup_MET","MET; MET (GeV);",100,0,1000);
	TH1F *scaleup_Mt = new TH1F("scaleup_Mt","M_{T}; M_{T} (GeV);",200,0,1000); 
	TH1F *scaleup_HT = new TH1F("scaleup_HT","HT; HT (GeV);",100,0,1200); 
	TH1F *scaleup_dPhiEleMET = new TH1F("scaleup_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	TH1F *scaledo_PhoEt = new TH1F("scaledo_PhoEt","#gamma E_{T}; E_{T} (GeV)",44,plotEtBins);
	TH1F *scaledo_PhoEta = new TH1F("scaledo_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1F *scaledo_LepPt = new TH1F("scaledo_LepPt","scaledo_LepPt",46,plotPtBins);
	TH1F *scaledo_LepEta = new TH1F("scaledo_LepEta","scaledo_LepEta",60,-3,3);
	TH1F *scaledo_MET = new TH1F("scaledo_MET","MET; MET (GeV);",100,0,1000);
	TH1F *scaledo_Mt = new TH1F("scaledo_Mt","M_{T}; M_{T} (GeV);",200,0,1000); 
	TH1F *scaledo_HT = new TH1F("scaledo_HT","HT; HT (GeV);",100,0,1200); 
	TH1F *scaledo_dPhiEleMET = new TH1F("scaledo_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	TH1F *normup_PhoEt = new TH1F("normup_PhoEt","#gamma E_{T}; E_{T} (GeV)",44,plotEtBins);
	TH1F *normup_PhoEta = new TH1F("normup_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1F *normup_LepPt = new TH1F("normup_LepPt","normup_LepPt",46,plotPtBins);
	TH1F *normup_LepEta = new TH1F("normup_LepEta","normup_LepEta",60,-3,3);
	TH1F *normup_MET = new TH1F("normup_MET","MET; MET (GeV);",100,0,1000);
	TH1F *normup_Mt = new TH1F("normup_Mt","M_{T}; M_{T} (GeV);",200,0,1000); 
	TH1F *normup_HT = new TH1F("normup_HT","HT; HT (GeV);",100,0,1200); 
	TH1F *normup_dPhiEleMET = new TH1F("normup_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	TH1F *normdo_PhoEt = new TH1F("normdo_PhoEt","#gamma E_{T}; E_{T} (GeV)",44,plotEtBins);
	TH1F *normdo_PhoEta = new TH1F("normdo_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1F *normdo_LepPt = new TH1F("normdo_LepPt","normdo_LepPt",46,plotPtBins);
	TH1F *normdo_LepEta = new TH1F("normdo_LepEta","normdo_LepEta",60,-3,3);
	TH1F *normdo_MET = new TH1F("normdo_MET","MET; MET (GeV);",100,0,1000);
	TH1F *normdo_Mt = new TH1F("normdo_Mt","M_{T}; M_{T} (GeV);",200,0,1000); 
	TH1F *normdo_HT = new TH1F("normdo_HT","HT; HT (GeV);",100,0,1200); 
	TH1F *normdo_dPhiEleMET = new TH1F("normdo_dPhiEleMET","dPhiEleMET",32,0,3.2); 
// ********  MC *************************//
  TChain *mctree;
	if(channelType == 1)mctree = new TChain("egTree","egTree");
  else if(channelType == 2)mctree = new TChain("mgTree","mgTree");
  //mctree->Add("/uscms_data/d3/mengleis/FullStatusData/resTree_VGamma_WGToLNuG.root");
  mctree->Add("/uscms_data/d3/mengleis/resTree_VGamma_WGToLNuG.root");
  mctree->Add("/uscms_data/d3/mengleis/FullStatusData/resTree_VGamma_WGToLNuG_130.root");
	float crosssection(0);
	float ntotalevent(0);
	float PUweight(1);
  float phoEt(0);
  float phoEta(0);
  float phoPhi(0);
  float lepPt(0);
  float lepEta(0);
  float lepPhi(0);
  float sigMT(0);
  float sigMET(0);
  float sigMETPhi(0);
  float dPhiLepMET(0);
  int   nVertex(0);
  float dRPhoLep(0);
  float HT(0);
  float nJet(0);
  float invmass(0);  
	float sigMETJESup(0);
	float sigMETJESdo(0);
	float sigMETJERup(0);
	float sigMETJERdo(0);
	float sigMTJESup(0);
	float sigMTJESdo(0);
	float sigMTJERup(0);
	float sigMTJERdo(0);
	float HTJESup(0);
	float HTJESdo(0);
	float dPhiLepMETJESup(0);
	float dPhiLepMETJESdo(0);
	float dPhiLepMETJERup(0);
	float dPhiLepMETJERdo(0);

	mctree->SetBranchAddress("crosssection",&crosssection);
	mctree->SetBranchAddress("ntotalevent", &ntotalevent);
	mctree->SetBranchAddress("PUweight",  &PUweight);
  mctree->SetBranchAddress("phoEt",     &phoEt);
  mctree->SetBranchAddress("phoEta",    &phoEta);
  mctree->SetBranchAddress("phoPhi",    &phoPhi);
  mctree->SetBranchAddress("lepPt",     &lepPt);
  mctree->SetBranchAddress("lepEta",    &lepEta);
  mctree->SetBranchAddress("lepPhi",    &lepPhi);
  mctree->SetBranchAddress("sigMT",     &sigMT);
  mctree->SetBranchAddress("sigMET",    &sigMET);
  mctree->SetBranchAddress("sigMETPhi", &sigMETPhi);
  mctree->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
  mctree->SetBranchAddress("nVertex",   &nVertex);
  mctree->SetBranchAddress("dRPhoLep",  &dRPhoLep);
  mctree->SetBranchAddress("HT",        &HT);
  mctree->SetBranchAddress("nJet",      &nJet);
  //mctree->SetBranchAddress("invmass",   &invmass);
	mctree->SetBranchAddress("sigMETJESup",     &sigMETJESup);
	mctree->SetBranchAddress("sigMETJESdo",     &sigMETJESdo);
	mctree->SetBranchAddress("sigMETJERup",     &sigMETJERup);
	mctree->SetBranchAddress("sigMETJERdo",     &sigMETJERdo);
	mctree->SetBranchAddress("sigMTJESup",      &sigMTJESup);
	mctree->SetBranchAddress("sigMTJESdo",      &sigMTJESdo);
	mctree->SetBranchAddress("sigMTJERup",      &sigMTJERup);
	mctree->SetBranchAddress("sigMTJERdo",      &sigMTJERdo);
	mctree->SetBranchAddress("HTJESup",         &HTJESup);
	mctree->SetBranchAddress("HTJESdo",         &HTJESdo);
	mctree->SetBranchAddress("dPhiLepMETJESup", &dPhiLepMETJESup);
	mctree->SetBranchAddress("dPhiLepMETJESdo", &dPhiLepMETJESdo);
	mctree->SetBranchAddress("dPhiLepMETJERup", &dPhiLepMETJERup);
	mctree->SetBranchAddress("dPhiLepMETJERdo", &dPhiLepMETJERdo);

	for(unsigned ievt(0); ievt < mctree->GetEntries(); ievt++){
		mctree->GetEntry(ievt);
		p_PU->Fill(nVertex,PUweight);
		double scalefactor(0);
		double scalefactorup(0);
		double scalefactordo(0);
		if(channelType == 1){
			scalefactor = objectESF.getElectronESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getegPhotonTRGESF(phoEt,phoEta)*objectESF.getElectronTRGESF(lepPt,lepEta);
			double s_ele_error = objectESF.getElectronESFError(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getegPhotonTRGESF(phoEt,phoEta)*objectESF.getElectronTRGESF(lepPt,lepEta);
			double s_pho_error = objectESF.getPhotonESFError(phoEt,phoEta)*objectESF.getElectronESF(lepPt,lepEta)*objectESF.getegPhotonTRGESF(phoEt,phoEta)*objectESF.getElectronTRGESF(lepPt,lepEta);
			double s_eletrg_error = objectESF.getElectronTRGESFError(lepPt,lepEta)*objectESF.getElectronESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getegPhotonTRGESF(phoEt,phoEta);
			double s_photrg_error = objectESF.getegPhotonTRGESFError(phoEt,phoEta)*objectESF.getElectronESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getElectronTRGESF(lepPt,lepEta);
			double s_error = sqrt(pow(s_ele_error,2) + pow(s_pho_error,2) + pow(s_eletrg_error,2) + pow(s_photrg_error, 2)); 
			scalefactorup = scalefactor + s_error; 
			scalefactordo = scalefactor - s_error;
		}
		if(channelType == 2){
			scalefactor = objectESF.getMuonESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getMuonEGTRGESF(phoEt, lepPt);
			double s_mu_error = objectESF.getMuonESFError(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta)*objectESF.getMuonEGTRGESF(phoEt, lepPt);
      double s_pho_error = objectESF.getPhotonESFError(phoEt,phoEta)*objectESF.getMuonESF(lepPt,lepEta)*objectESF.getMuonEGTRGESF(phoEt, lepPt);
      double s_trg_error = objectESF.getMuonEGTRGESFError(phoEt, lepPt)*objectESF.getMuonESF(lepPt,lepEta)*objectESF.getPhotonESF(phoEt,phoEta);
			double s_error = sqrt(pow(s_mu_error,2) + pow(s_pho_error,2) + pow(s_trg_error,2));
			scalefactorup = scalefactor + s_error; 
			scalefactordo = scalefactor - s_error;
		}
		float XS_weight = 35.87*1000*crosssection/ntotalevent;
		if(crosssection < 1)XS_weight = XS_weight*ZGammaFactor;
		if(crosssection > 1 && crosssection < 2)XS_weight = XS_weight*WGammaFactor;

//		if(crosssection < 1)XS_weight = XS_weight*0.58;
//		else if(crosssection > 1 && crosssection < 2)XS_weight = XS_weight*0.62;
//		else if(crosssection > 2 && crosssection < 200)XS_weight = XS_weight*0.58;
//		else if(crosssection > 400)XS_weight = XS_weight*0.62;
//
		float weight = PUweight*XS_weight*scalefactor;
		float weight_scaleup = PUweight*XS_weight*scalefactorup;
		float weight_scaledo = PUweight*XS_weight*scalefactordo;


		bool doprocess(true);
		if(crosssection < 1 && phoEt >=140)p_PhoEt_high_ZG->Fill(phoEt, XS_weight);
		else if(crosssection > 2 && crosssection < 200 && phoEt < 140)p_PhoEt_low_ZG->Fill(phoEt, XS_weight);
		if(crosssection > 1 && crosssection < 2 && phoEt >=140)p_PhoEt_high_WG->Fill(phoEt, XS_weight);
		else if(crosssection > 400 && phoEt < 140)p_PhoEt_low_WG->Fill(phoEt, XS_weight);

		if(crosssection < 2 && phoEt < 140)doprocess = false;
		if(crosssection > 10 && phoEt >= 140)doprocess = false;
		if(!doprocess)continue;

		if(sigMET > met_upcut || sigMET < met_lowcut)continue;
		if(phoEt < 35 || lepPt < 25)continue;
		if(fabs(phoEta) > 1.4442 || fabs(lepEta) > 2.5)continue;

		if(lepPt < leplow || lepPt > lephigh)continue;
		//if(fabs(lepEta) > 1.4442)continue;
		//if(phoEt < 80)continue;
		//if(sigMT > 100)continue;

//		std::cout << "weight = " << weight*factorMCE << " pt = " << phoEt <<  " XS = " << XS_weight << " scale = " << scalefactor << " PU= " << PUweight << std::endl;

		p_PhoEt->Fill(phoEt, weight*factorMCE);
		p_PhoEta->Fill(phoEta,weight*factorMCE);
		p_LepPt->Fill(lepPt, weight*factorMCE);
		p_LepEta->Fill(lepEta,weight*factorMCE);
		p_MET->Fill(sigMET, weight*factorMCE);
		p_Mt->Fill(sigMT, weight*factorMCE);
		p_HT->Fill(HT, weight*factorMCE);
		p_dPhiEleMET->Fill(fabs(dPhiLepMET), weight*factorMCE);

		jesup_MET->Fill(sigMETJESup, weight*factorMCE);
		jesup_Mt->Fill(sigMTJESup, weight*factorMCE);
		jesup_HT->Fill(HTJESup, weight*factorMCE);
		jesup_dPhiEleMET->Fill(fabs(dPhiLepMETJESup), weight*factorMCE);

		jesdo_MET->Fill(sigMETJESdo, weight*factorMCE);
		jesdo_Mt->Fill(sigMTJESdo, weight*factorMCE);
		jesdo_HT->Fill(HTJESdo, weight*factorMCE);
		jesdo_dPhiEleMET->Fill(fabs(dPhiLepMETJESdo), weight*factorMCE);

		jerup_MET->Fill(sigMETJERup, weight*factorMCE);
		jerup_Mt->Fill(sigMTJERup, weight*factorMCE);
		jerup_dPhiEleMET->Fill(fabs(dPhiLepMETJERup), weight*factorMCE);

		jerdo_MET->Fill(sigMETJERdo, weight*factorMCE);
		jerdo_Mt->Fill(sigMTJERdo, weight*factorMCE);
		jerdo_dPhiEleMET->Fill(fabs(dPhiLepMETJERdo), weight*factorMCE);

		scaleup_PhoEt->Fill(phoEt, weight_scaleup*factorMCE);
		scaleup_PhoEta->Fill(phoEta,weight_scaleup*factorMCE);
		scaleup_LepPt->Fill(lepPt, weight_scaleup*factorMCE);
		scaleup_LepEta->Fill(lepEta,weight_scaleup*factorMCE);
		scaleup_MET->Fill(sigMET, weight_scaleup*factorMCE);
		scaleup_Mt->Fill(sigMT, weight_scaleup*factorMCE);
		scaleup_HT->Fill(HT, weight_scaleup*factorMCE);
		scaleup_dPhiEleMET->Fill(fabs(dPhiLepMET), weight_scaleup*factorMCE);

		scaledo_PhoEt->Fill(phoEt, weight_scaledo*factorMCE);
		scaledo_PhoEta->Fill(phoEta,weight_scaledo*factorMCE);
		scaledo_LepPt->Fill(lepPt, weight_scaledo*factorMCE);
		scaledo_LepEta->Fill(lepEta,weight_scaledo*factorMCE);
		scaledo_MET->Fill(sigMET, weight_scaledo*factorMCE);
		scaledo_Mt->Fill(sigMT, weight_scaledo*factorMCE);
		scaledo_HT->Fill(HT, weight_scaledo*factorMCE);
		scaledo_dPhiEleMET->Fill(fabs(dPhiLepMET), weight_scaledo*factorMCE);

		normup_PhoEt->Fill(phoEt, weight*factorMCEUP);
		normup_PhoEta->Fill(phoEta,weight*factorMCEUP);
		normup_LepPt->Fill(lepPt, weight*factorMCEUP);
		normup_LepEta->Fill(lepEta,weight*factorMCEUP);
		normup_MET->Fill(sigMET, weight*factorMCEUP);
		normup_Mt->Fill(sigMT, weight*factorMCEUP);
		normup_HT->Fill(HT, weight*factorMCEUP);
		normup_dPhiEleMET->Fill(fabs(dPhiLepMET), weight*factorMCEUP);

		normdo_PhoEt->Fill(phoEt, weight*factorMCEDO);
		normdo_PhoEta->Fill(phoEta,weight*factorMCEDO);
		normdo_LepPt->Fill(lepPt, weight*factorMCEDO);
		normdo_LepEta->Fill(lepEta,weight*factorMCEDO);
		normdo_MET->Fill(sigMET, weight*factorMCEDO);
		normdo_Mt->Fill(sigMT, weight*factorMCEDO);
		normdo_HT->Fill(HT, weight*factorMCEDO);
		normdo_dPhiEleMET->Fill(fabs(dPhiLepMET), weight*factorMCEDO);
	}
	for(int ibin(1); ibin < p_PhoEt->GetSize(); ibin++){
		double syserror(0);
		syserror += p_PhoEt->GetBinError(ibin)* p_PhoEt->GetBinError(ibin);
		syserror += pow((scaleup_PhoEt->GetBinContent(ibin)-p_PhoEt->GetBinContent(ibin)),2);
		p_PhoEt->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_LepPt->GetSize(); ibin++){
		double syserror(0);
		syserror += p_LepPt->GetBinError(ibin)* p_LepPt->GetBinError(ibin);
		syserror += pow((scaleup_LepPt->GetBinContent(ibin)-p_LepPt->GetBinContent(ibin)),2);
		p_LepPt->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_MET->GetSize(); ibin++){
		double syserror(0);
		syserror += p_MET->GetBinError(ibin)* p_MET->GetBinError(ibin);
		syserror += pow((scaleup_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)),2);
		p_MET->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_Mt->GetSize(); ibin++){
		double syserror(0);
		syserror += p_Mt->GetBinError(ibin)* p_Mt->GetBinError(ibin);
		syserror += pow((scaleup_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)),2);
		p_Mt->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_HT->GetSize(); ibin++){
		double syserror(0);
		syserror += p_HT->GetBinError(ibin)* p_HT->GetBinError(ibin);
		syserror += pow((scaleup_HT->GetBinContent(ibin)-p_HT->GetBinContent(ibin)),2);
		p_HT->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_dPhiEleMET->GetSize(); ibin++){
		double syserror(0);
		syserror += p_dPhiEleMET->GetBinError(ibin)* p_dPhiEleMET->GetBinError(ibin);
		syserror += pow((scaleup_dPhiEleMET->GetBinContent(ibin)-p_dPhiEleMET->GetBinContent(ibin)),2);
		p_dPhiEleMET->SetBinError(ibin,sqrt(syserror));
	}	

	outputfile->Write();
	outputfile->Close();
}


