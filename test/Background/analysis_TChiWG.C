#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TF3.h"
#include "TH1D.h"
#include "TH2D.h"
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

void analysis_TChiWG(){//main  

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

	std::string conftype;
	double confvalue;

	std::ifstream binfile("binConfig.txt");
	float METbin1(200), METbin2(300);
	if(binfile.is_open()){
		for(int i(0); i<2; i++){
			binfile >> conftype >> confvalue;
			if(conftype.find("METbin1")!=std::string::npos)METbin1= confvalue;
			if(conftype.find("METbin2")!=std::string::npos)METbin2= confvalue;
		}
	}

	TFile xSecFile("../cross/susyCrossSection.root");
	TH1D *p_crosssection_tchiwg = (TH1D*)xSecFile.Get("p_charginoSec");
	TH1D *p_crosssection_t5wg = (TH1D*)xSecFile.Get("p_gluinoxSec");

	TH2D *p_char_neu = new TH2D("p_char_neu","p_char_neu",80, 12.5, 2012.5, 80, 12.5, 2012.5);

	std::ostringstream histname;
	//**************   T5WG  ***************************//


  TFile *file_t5wg = TFile::Open("/uscms_data/d3/mengleis/Sep1/resTree_T5WG.root");
  TTree *tree_t5wg = (TTree*)file_t5wg->Get("SUSYtree");
	float Mgluino_t5wg(0);
  float Mchargino_t5wg(0);
  float Mneutralino_t5wg(0);
	tree_t5wg->SetBranchAddress("Mgluino",    &Mgluino_t5wg);  
  tree_t5wg->SetBranchAddress("Mchargino",  &Mchargino_t5wg);
  tree_t5wg->SetBranchAddress("Mneutralino",&Mneutralino_t5wg);

	TFile *outputfile_t5wg = TFile::Open("signalTree_T5WG.root","RECREATE");
	outputfile_t5wg->cd();

	TH2D *p_T5WGMASS = new TH2D("SUSYMass","",27, 775.0, 2125.0, 80, 12.5, 2012.5);
	TH2D *p_T5WGselect = new TH2D("p_T5WGselect","",27, 775.0, 2125.0, 80, 12.5, 2012.5); 

	double WGevent(0);	
	double WWevent(0);	
	double GGevent(0);
	for(unsigned ievt(0); ievt < tree_t5wg->GetEntries(); ievt++){
		tree_t5wg->GetEntry(ievt);
    if (ievt%1000000==0) std::cout << " -- Processing event " << ievt << std::endl;
		if(Mgluino_t5wg >0 && Mchargino_t5wg >0 && Mneutralino_t5wg >0){
			p_T5WGMASS->Fill(Mgluino_t5wg, Mchargino_t5wg);
		}
		p_char_neu->Fill(Mchargino_t5wg, Mneutralino_t5wg);

		if(Mgluino_t5wg >0 && Mchargino_t5wg >0 && Mneutralino_t5wg <= 0)WWevent += 1;
		if(Mgluino_t5wg >0 && Mchargino_t5wg <= 0 && Mneutralino_t5wg > 0)GGevent += 1;
		if(Mgluino_t5wg >0 && Mchargino_t5wg >0 && Mneutralino_t5wg > 0)WGevent += 1;
		
	}
	std::cout << "WW " << WWevent << "  GG " << GGevent << "   WG " << WGevent << std::endl;

  //p_T5WGMASS->Draw();
	TH2D *t5wg_h_chan_rate_nom[18]; 
	TH2D *t5wg_h_chan_rate_jesUp[18]; 
	TH2D *t5wg_h_chan_rate_jesDown[18];    
	TH2D *t5wg_h_chan_rate_jerUp[18];    
	TH2D *t5wg_h_chan_rate_jerDown[18];    
	TH2D *t5wg_h_chan_rate_xsUp[18];       
                                  
	TH2D *t5wg_h_chan_syserr_jes[18];      
	TH2D *t5wg_h_chan_syserr_jer[18];     
	TH2D *t5wg_h_chan_syserr_esf[18];     
	TH2D *t5wg_h_chan_syserr_scale[18];   
	TH2D *t5wg_h_chan_syserr_eleshape[18]; 
	TH2D *t5wg_h_chan_syserr_jetshape[18];
	TH2D *t5wg_h_chan_syserr_xs[18];
	TH2D *t5wg_h_chan_syserr_lumi[18];     

	for(unsigned i(0); i < 18; i++){
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_nom";
		t5wg_h_chan_rate_nom[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_jesUp";
		t5wg_h_chan_rate_jesUp[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_jesDown";
		t5wg_h_chan_rate_jesDown[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_jerUp";
		t5wg_h_chan_rate_jerUp[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_jerDown";
		t5wg_h_chan_rate_jerDown[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_xsUp";
		t5wg_h_chan_rate_xsUp[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
															
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_jes";
		t5wg_h_chan_syserr_jes[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_jer";
		t5wg_h_chan_syserr_jer[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_esf";
		t5wg_h_chan_syserr_esf[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_scale";
		t5wg_h_chan_syserr_scale[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_eleshape";
		t5wg_h_chan_syserr_eleshape[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_jetshape";
		t5wg_h_chan_syserr_jetshape[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_xs";
		t5wg_h_chan_syserr_xs[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_lumi";
		t5wg_h_chan_syserr_lumi[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
	}
	
	TH1D *p_t5wg_MET_valid[12][20];
	TH1D *p_t5wg_MET_signal[12][20];
	TH1D *p_t5wg_HT_signal[12][20];
	for(unsigned i(0); i < 12; i++){
		for(unsigned j(0); j < 20; j++){
			histname.str("");
			histname << "p_t5wg_MET_valid_" << int(800+i*100) << "_" << int(100+j*100);
			p_t5wg_MET_valid[i][j] = new TH1D(histname.str().c_str(),histname.str().c_str(),nBkgMETBins, bkgMETBins);
			histname.str("");
			histname << "p_t5wg_MET_signal_" << int(800+i*100) << "_" << int(100+j*100);
			p_t5wg_MET_signal[i][j] = new TH1D(histname.str().c_str(),histname.str().c_str(),nSigMETBins, sigMETBins);
			histname.str("");
			histname << "p_t5wg_HT_signal_" << int(800+i*100) << "_" << int(100+j*100);
			p_t5wg_HT_signal[i][j] = new TH1D(histname.str().c_str(),histname.str().c_str(),nSigHTBins, sigHTBins);
		}
	}
		
  TChain *mgtree_t5wg;
  mgtree_t5wg = new TChain("mgTree","mgTree");
  mgtree_t5wg->Add("/uscms_data/d3/mengleis/Sep1/resTree_T5WG.root");
  float phoEt_t5wg_mg(0);
  float phoEta_t5wg_mg(0);
  float lepPt_t5wg_mg(0);
  float lepEta_t5wg_mg(0);
  float sigMT_t5wg_mg(0);
  float sigMET_t5wg_mg(0);
  float HT_t5wg_mg(0);
	float sigMETJESup_t5wg_mg(0);
	float sigMETJESdo_t5wg_mg(0);
	float sigMETJERup_t5wg_mg(0);
	float sigMETJERdo_t5wg_mg(0);
	float sigMTJESup_t5wg_mg(0);
	float sigMTJESdo_t5wg_mg(0);
	float sigMTJERup_t5wg_mg(0);
	float sigMTJERdo_t5wg_mg(0);
	float HTJESup_t5wg_mg(0);
	float HTJESdo_t5wg_mg(0);
	float gluinoMass_t5wg_mg(0);
  float charginoMass_t5wg_mg(0);
  float neutralinoMass_t5wg_mg(0);
  mgtree_t5wg->SetBranchAddress("phoEt",      &phoEt_t5wg_mg);
  mgtree_t5wg->SetBranchAddress("phoEta",     &phoEta_t5wg_mg);
  mgtree_t5wg->SetBranchAddress("lepPt",      &lepPt_t5wg_mg);
  mgtree_t5wg->SetBranchAddress("lepEta",     &lepEta_t5wg_mg);
  mgtree_t5wg->SetBranchAddress("sigMT",      &sigMT_t5wg_mg);
  mgtree_t5wg->SetBranchAddress("sigMET",     &sigMET_t5wg_mg);
  mgtree_t5wg->SetBranchAddress("HT",         &HT_t5wg_mg);
	mgtree_t5wg->SetBranchAddress("Mgluino",    &gluinoMass_t5wg_mg);
  mgtree_t5wg->SetBranchAddress("Mchargino",  &charginoMass_t5wg_mg);
  mgtree_t5wg->SetBranchAddress("Mneutralino",&neutralinoMass_t5wg_mg);
	mgtree_t5wg->SetBranchAddress("sigMETJESup",&sigMETJESup_t5wg_mg);
	mgtree_t5wg->SetBranchAddress("sigMETJESdo",&sigMETJESdo_t5wg_mg);
	mgtree_t5wg->SetBranchAddress("sigMETJERup",&sigMETJERup_t5wg_mg);
	mgtree_t5wg->SetBranchAddress("sigMETJERdo",&sigMETJERdo_t5wg_mg);
	mgtree_t5wg->SetBranchAddress("sigMTJESup", &sigMTJESup_t5wg_mg);
	mgtree_t5wg->SetBranchAddress("sigMTJESdo", &sigMTJESdo_t5wg_mg);
	mgtree_t5wg->SetBranchAddress("sigMTJERup", &sigMTJERup_t5wg_mg);
	mgtree_t5wg->SetBranchAddress("sigMTJERdo", &sigMTJERdo_t5wg_mg);
	mgtree_t5wg->SetBranchAddress("HTJESup",    &HTJESup_t5wg_mg);
	mgtree_t5wg->SetBranchAddress("HTJESdo",    &HTJESdo_t5wg_mg);

	for(unsigned ievt(0); ievt < mgtree_t5wg->GetEntries(); ievt++){
		mgtree_t5wg->GetEntry(ievt);

		if(charginoMass_t5wg_mg <=0 || neutralinoMass_t5wg_mg <=0)continue;
    p_T5WGselect->Fill(gluinoMass_t5wg_mg, charginoMass_t5wg_mg);
		/** cut flow *****/
		if(phoEt_t5wg_mg < 35 || lepPt_t5wg_mg < 25)continue;
		if(fabs(phoEta_t5wg_mg) > 1.4442 || fabs(lepEta_t5wg_mg) > 2.5)continue;

		for(unsigned i(0); i < 12; i++){
			for(unsigned j(0); j < 20; j++){
				if( fabs(gluinoMass_t5wg_mg - (800+i*100) ) < 2 && fabs(charginoMass_t5wg_mg - (100+j*100) ) < 2){
					if(sigMT_t5wg_mg < 100)p_t5wg_MET_valid[i][j]->Fill(sigMET_t5wg_mg);
					else{
					 p_t5wg_MET_signal[i][j]->Fill(sigMET_t5wg_mg);
					 p_t5wg_HT_signal[i][j]->Fill(HT_t5wg_mg);
					}
				}
			}
		}

		if(sigMET_t5wg_mg < 120 || sigMT_t5wg_mg < 100)continue;

		int SigBinIndex(-1);
		SigBinIndex = findSignalBin(sigMET_t5wg_mg, HT_t5wg_mg, METbin1, METbin2);
		if(SigBinIndex >=0){
			t5wg_h_chan_rate_nom[SigBinIndex]->Fill( gluinoMass_t5wg_mg, charginoMass_t5wg_mg, 1);
			t5wg_h_chan_rate_xsUp[SigBinIndex]->Fill( gluinoMass_t5wg_mg, charginoMass_t5wg_mg, 1); 
		}
		int jesupBinIndex(-1);	
		jesupBinIndex = findSignalBin(sigMETJESup_t5wg_mg, HTJESup_t5wg_mg, METbin1, METbin2);
		if(jesupBinIndex >=0){
			t5wg_h_chan_rate_jesUp[jesupBinIndex]->Fill( gluinoMass_t5wg_mg, charginoMass_t5wg_mg, 1); 
		}
		int jesdoBinIndex(-1);
		jesdoBinIndex = findSignalBin(sigMETJESdo_t5wg_mg, HTJESdo_t5wg_mg, METbin1, METbin2);
		if(jesdoBinIndex >=0){
			t5wg_h_chan_rate_jesDown[jesdoBinIndex]->Fill( gluinoMass_t5wg_mg, charginoMass_t5wg_mg, 1);   
		}	
		int jerupBinIndex(-1);
		jerupBinIndex = findSignalBin(sigMETJERup_t5wg_mg, HT_t5wg_mg, METbin1, METbin2);
		if( jerupBinIndex >=0){
			t5wg_h_chan_rate_jerUp[jerupBinIndex]->Fill( gluinoMass_t5wg_mg, charginoMass_t5wg_mg, 1);
		}
		int jerdoBinIndex(-1);
		jerdoBinIndex = findSignalBin(sigMETJERdo_t5wg_mg, HT_t5wg_mg, METbin1, METbin2);	
		if(jerdoBinIndex >= 0){
			t5wg_h_chan_rate_jerDown[jerdoBinIndex]->Fill( gluinoMass_t5wg_mg, charginoMass_t5wg_mg, 1);
		}  
	}

  TChain *egtree_t5wg;
  egtree_t5wg = new TChain("egTree","egTree");
  egtree_t5wg->Add("/uscms_data/d3/mengleis/Sep1/resTree_T5WG.root");
  float phoEt_t5wg_eg(0);
  float phoEta_t5wg_eg(0);
  float lepPt_t5wg_eg(0);
  float lepEta_t5wg_eg(0);
  float sigMT_t5wg_eg(0);
  float sigMET_t5wg_eg(0);
  float HT_t5wg_eg(0);
	float sigMETJESup_t5wg_eg(0);
	float sigMETJESdo_t5wg_eg(0);
	float sigMETJERup_t5wg_eg(0);
	float sigMETJERdo_t5wg_eg(0);
	float sigMTJESup_t5wg_eg(0);
	float sigMTJESdo_t5wg_eg(0);
	float sigMTJERup_t5wg_eg(0);
	float sigMTJERdo_t5wg_eg(0);
	float HTJESup_t5wg_eg(0);
	float HTJESdo_t5wg_eg(0);
	float gluinoMass_t5wg_eg(0);
  float charginoMass_t5wg_eg(0);
  float neutralinoMass_t5wg_eg(0);
  egtree_t5wg->SetBranchAddress("phoEt",      &phoEt_t5wg_eg);
  egtree_t5wg->SetBranchAddress("phoEta",     &phoEta_t5wg_eg);
  egtree_t5wg->SetBranchAddress("lepPt",      &lepPt_t5wg_eg);
  egtree_t5wg->SetBranchAddress("lepEta",     &lepEta_t5wg_eg);
  egtree_t5wg->SetBranchAddress("sigMT",      &sigMT_t5wg_eg);
  egtree_t5wg->SetBranchAddress("sigMET",     &sigMET_t5wg_eg);
  egtree_t5wg->SetBranchAddress("HT",         &HT_t5wg_eg);
	egtree_t5wg->SetBranchAddress("Mgluino",    &gluinoMass_t5wg_eg);
  egtree_t5wg->SetBranchAddress("Mchargino",  &charginoMass_t5wg_eg);
  egtree_t5wg->SetBranchAddress("Mneutralino",&neutralinoMass_t5wg_eg);
	egtree_t5wg->SetBranchAddress("sigMETJESup",&sigMETJESup_t5wg_eg);
	egtree_t5wg->SetBranchAddress("sigMETJESdo",&sigMETJESdo_t5wg_eg);
	egtree_t5wg->SetBranchAddress("sigMETJERup",&sigMETJERup_t5wg_eg);
	egtree_t5wg->SetBranchAddress("sigMETJERdo",&sigMETJERdo_t5wg_eg);
	egtree_t5wg->SetBranchAddress("sigMTJESup", &sigMTJESup_t5wg_eg);
	egtree_t5wg->SetBranchAddress("sigMTJESdo", &sigMTJESdo_t5wg_eg);
	egtree_t5wg->SetBranchAddress("sigMTJERup", &sigMTJERup_t5wg_eg);
	egtree_t5wg->SetBranchAddress("sigMTJERdo", &sigMTJERdo_t5wg_eg);
	egtree_t5wg->SetBranchAddress("HTJESup",    &HTJESup_t5wg_eg);
	egtree_t5wg->SetBranchAddress("HTJESdo",    &HTJESdo_t5wg_eg);

	for(unsigned ievt(0); ievt < egtree_t5wg->GetEntries(); ievt++){
		egtree_t5wg->GetEntry(ievt);

		/** cut flow *****/
		if(phoEt_t5wg_eg < 35 || lepPt_t5wg_eg < 25)continue;
		if(fabs(phoEta_t5wg_eg) > 1.4442 || fabs(lepEta_t5wg_eg) > 2.5)continue;
		if(sigMET_t5wg_eg < 120 || sigMT_t5wg_eg < 100)continue;

		int SigBinIndex(-1);
		SigBinIndex = findSignalBin(sigMET_t5wg_eg, HT_t5wg_eg, METbin1, METbin2) + 9;
		if(SigBinIndex >=0){
			t5wg_h_chan_rate_nom[SigBinIndex]->Fill( gluinoMass_t5wg_eg, charginoMass_t5wg_eg, 1);
			t5wg_h_chan_rate_xsUp[SigBinIndex]->Fill( gluinoMass_t5wg_eg, charginoMass_t5wg_eg, 1); 
		}
		int jesupBinIndex(-1);	
		jesupBinIndex = findSignalBin(sigMETJESup_t5wg_eg, HTJESup_t5wg_eg, METbin1, METbin2) + 9;
		if(jesupBinIndex >=0){
			t5wg_h_chan_rate_jesUp[jesupBinIndex]->Fill( gluinoMass_t5wg_eg, charginoMass_t5wg_eg, 1); 
		}
		int jesdoBinIndex(-1);
		jesdoBinIndex = findSignalBin(sigMETJESdo_t5wg_eg, HTJESdo_t5wg_eg, METbin1, METbin2) + 9;
		if(jesdoBinIndex >=0){
			t5wg_h_chan_rate_jesDown[jesdoBinIndex]->Fill( gluinoMass_t5wg_eg, charginoMass_t5wg_eg, 1);   
		}	
		int jerupBinIndex(-1);
		jerupBinIndex = findSignalBin(sigMETJERup_t5wg_eg, HT_t5wg_eg, METbin1, METbin2) + 9;
		if( jerupBinIndex >=0){
			t5wg_h_chan_rate_jerUp[jerupBinIndex]->Fill( gluinoMass_t5wg_eg, charginoMass_t5wg_eg, 1);
		}
		int jerdoBinIndex(-1);
		jerdoBinIndex = findSignalBin(sigMETJERdo_t5wg_eg, HT_t5wg_eg, METbin1, METbin2) + 9;	
		if(jerdoBinIndex >= 0){
			t5wg_h_chan_rate_jerDown[jerdoBinIndex]->Fill( gluinoMass_t5wg_eg, charginoMass_t5wg_eg, 1);
		}  
	}


	for(unsigned ih(0); ih < 18; ih++){
		for(unsigned i(1); i < t5wg_h_chan_rate_nom[ih]->GetXaxis()->GetNbins() + 1; i++){
			for(unsigned j(1); j < t5wg_h_chan_rate_nom[ih]->GetYaxis()->GetNbins() + 1; j++){
				if(p_T5WGMASS->GetBinContent(i,j) < 1000){
					p_T5WGselect->SetBinContent(i,j,-1);
					p_T5WGMASS->SetBinContent(i,j,-1);
				}
		
				if(p_T5WGMASS->GetBinContent(i,j) <= 0){
					t5wg_h_chan_rate_nom[ih]->SetBinContent(i,j, -1); 
					t5wg_h_chan_rate_jesUp[ih]->SetBinContent(i,j,-1);
					t5wg_h_chan_rate_jesDown[ih]->SetBinContent(i,j, -1);
					t5wg_h_chan_rate_jerUp[ih]->SetBinContent(i,j, -1);
					t5wg_h_chan_rate_jerDown[ih]->SetBinContent(i,j,-1);
					t5wg_h_chan_rate_xsUp[ih]->SetBinContent(i,j, -1);
							
					t5wg_h_chan_syserr_jes[ih]->SetBinContent(i,j, -1); 
					t5wg_h_chan_syserr_jer[ih]->SetBinContent(i,j, -1); 
					t5wg_h_chan_syserr_esf[ih]->SetBinContent(i,j, -1); 
					t5wg_h_chan_syserr_scale[ih]->SetBinContent(i,j, -1); 
					t5wg_h_chan_syserr_eleshape[ih]->SetBinContent(i,j, -1); 
					t5wg_h_chan_syserr_jetshape[ih]->SetBinContent(i,j, -1); 
					t5wg_h_chan_syserr_xs[ih]->SetBinContent(i,j, -1); 
					t5wg_h_chan_syserr_lumi[ih]->SetBinContent(i,j, -1); 
				}
				else{ 
					float noe = p_T5WGMASS->GetBinContent(i,j);
					float sparticleMass = p_T5WGMASS->GetXaxis()->GetBinCenter(i);
					float crosssection = 0.5*p_crosssection_t5wg->GetBinContent( p_crosssection_t5wg->FindBin(sparticleMass) );
					float crosssectionUp = (crosssection+0.5*p_crosssection_t5wg->GetBinError( p_crosssection_t5wg->FindBin(sparticleMass) ) ); 

					t5wg_h_chan_rate_nom[ih]->SetBinContent(i,j, t5wg_h_chan_rate_nom[ih]->GetBinContent(i,j)*35.8*crosssection/noe); 
					t5wg_h_chan_rate_jesUp[ih]->SetBinContent(i,j, t5wg_h_chan_rate_jesUp[ih]->GetBinContent(i,j)*35.8*crosssection/noe); 
					t5wg_h_chan_rate_jesDown[ih]->SetBinContent(i,j, t5wg_h_chan_rate_jesDown[ih]->GetBinContent(i,j)*35.8*crosssection/noe);
					t5wg_h_chan_rate_jerUp[ih]->SetBinContent(i,j, t5wg_h_chan_rate_jerUp[ih]->GetBinContent(i,j)*35.8*crosssection/noe);
					t5wg_h_chan_rate_jerDown[ih]->SetBinContent(i,j, t5wg_h_chan_rate_jerDown[ih]->GetBinContent(i,j)*35.8*crosssection/noe);
					t5wg_h_chan_rate_xsUp[ih]->SetBinContent(i,j, t5wg_h_chan_rate_xsUp[ih]->GetBinContent(i,j)*35.8*crosssectionUp/noe);
									
					t5wg_h_chan_syserr_jes[ih]->SetBinContent(i,j, max( fabs(t5wg_h_chan_rate_jesUp[ih]->GetBinContent(i,j)-t5wg_h_chan_rate_nom[ih]->GetBinContent(i,j)), fabs(t5wg_h_chan_rate_jesDown[ih]->GetBinContent(i,j)-t5wg_h_chan_rate_nom[ih]->GetBinContent(i,j)) ) ); 
					t5wg_h_chan_syserr_jer[ih]->SetBinContent(i,j, max( fabs(t5wg_h_chan_rate_jerUp[ih]->GetBinContent(i,j)-t5wg_h_chan_rate_nom[ih]->GetBinContent(i,j)), fabs(t5wg_h_chan_rate_jerDown[ih]->GetBinContent(i,j)-t5wg_h_chan_rate_nom[ih]->GetBinContent(i,j)) ) ); 
					t5wg_h_chan_syserr_esf[ih]->SetBinContent(i,j, -1);
					t5wg_h_chan_syserr_scale[ih]->SetBinContent(i,j, -1);
					t5wg_h_chan_syserr_eleshape[ih]->SetBinContent(i,j, -1);
					t5wg_h_chan_syserr_jetshape[ih]->SetBinContent(i,j, -1);
					t5wg_h_chan_syserr_xs[ih]->SetBinContent(i,j, fabs( t5wg_h_chan_rate_xsUp[ih]->GetBinContent(i,j)-t5wg_h_chan_rate_nom[ih]->GetBinContent(i,j)) );
					t5wg_h_chan_syserr_lumi[ih]->SetBinContent(i,j, 0.026*t5wg_h_chan_rate_nom[ih]->GetBinContent(i,j));   

				}
			}
		} 
	} 

	for(unsigned ibin(1); ibin < p_T5WGMASS->GetXaxis()->GetNbins() + 1; ibin++){
		for(unsigned jbin(1); jbin < p_T5WGMASS->GetYaxis()->GetNbins() + 1; jbin++){
	
			if(p_T5WGMASS->GetBinContent(ibin,jbin) > 0){
				float noe = p_T5WGMASS->GetBinContent(ibin,jbin);
				float sparticleMass = p_T5WGMASS->GetXaxis()->GetBinCenter(ibin);
				float crosssection = p_crosssection_t5wg->GetBinContent( p_crosssection_t5wg->FindBin(sparticleMass) );
				for(unsigned i(0); i < 12; i++){
					for(unsigned j(0); j < 20; j++){
						if( fabs(sparticleMass - (800+i*100) ) < 10 && fabs( p_T5WGMASS->GetYaxis()->GetBinCenter(jbin) - (100+j*100) ) < 10){
							std::cout << "mass " << sparticleMass << " xs " << crosssection << " noe " << noe << std::endl;
							p_t5wg_MET_valid[i][j]->Scale(35.8*crosssection/noe);
							p_t5wg_MET_signal[i][j]->Scale(35.8*crosssection/noe);
							p_t5wg_HT_signal[i][j]->Scale(35.8*crosssection/noe);
						}
					}
				}

			}
		}
	} 

	outputfile_t5wg->Write();
	outputfile_t5wg->Close();


	//****************   TChiWG ***************************//
  TFile *file_tchiwg = TFile::Open("/uscms_data/d3/mengleis/Sep1/resTree_TChiWG.root");
  TTree *tree_tchiwg = (TTree*)file_tchiwg->Get("SUSYtree");
	float Mgluino_tchiwg(0);
  float Mchargino_tchiwg(0);
  float Mneutralino_tchiwg(0);
	tree_tchiwg->SetBranchAddress("Mgluino",    &Mgluino_tchiwg);  
  tree_tchiwg->SetBranchAddress("Mchargino",  &Mchargino_tchiwg);
  tree_tchiwg->SetBranchAddress("Mneutralino",&Mneutralino_tchiwg);

	TFile *outputfile_tchiwg = TFile::Open("signalTree_TChiWG.root","RECREATE");
	outputfile_tchiwg->cd();

	TH1D *p_TChiWGMASS = new TH1D("p_TChiWGMASS","",40,287.5,1287.5);

	TH1D *p_tchiwg_MET_valid[10];
	TH1D *p_tchiwg_MET_signal[10];
	TH1D *p_tchiwg_HT_signal[10];
	for(unsigned i(0); i < 10; i++){
			histname.str("");
			histname << "p_tchiwg_MET_valid_" << int(300+i*100);
			p_tchiwg_MET_valid[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),nBkgMETBins, bkgMETBins);
			histname.str("");
			histname << "p_tchiwg_MET_signal_" << int(300+i*100);
			p_tchiwg_MET_signal[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),nSigMETBins, sigMETBins);
			histname.str("");
			histname << "p_tchiwg_HT_signal_" << int(300+i*100);
			p_tchiwg_HT_signal[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),nSigHTBins, sigHTBins);
	}

	for(unsigned ievt(0); ievt < tree_tchiwg->GetEntries(); ievt++){
		tree_tchiwg->GetEntry(ievt);
    if (ievt%1000000==0) std::cout << " -- Processing event " << ievt << std::endl;
		if(Mchargino_tchiwg >=0){
			p_TChiWGMASS->Fill(Mchargino_tchiwg);
		}
	}

  //p_T5WGMASS->Draw();
	TH1D *tchiwg_h_chan_rate_nom[18]; 
	TH1D *tchiwg_h_chan_rate_jesUp[18]; 
	TH1D *tchiwg_h_chan_rate_jesDown[18];    
	TH1D *tchiwg_h_chan_rate_jerUp[18];    
	TH1D *tchiwg_h_chan_rate_jerDown[18];    
	TH1D *tchiwg_h_chan_rate_xsUp[18];       
                                  
	TH1D *tchiwg_h_chan_syserr_jes[18];      
	TH1D *tchiwg_h_chan_syserr_jer[18];     
	TH1D *tchiwg_h_chan_syserr_esf[18];     
	TH1D *tchiwg_h_chan_syserr_scale[18];   
	TH1D *tchiwg_h_chan_syserr_eleshape[18]; 
	TH1D *tchiwg_h_chan_syserr_jetshape[18];
	TH1D *tchiwg_h_chan_syserr_xs[18];
	TH1D *tchiwg_h_chan_syserr_lumi[18];     

	for(unsigned i(0); i < 18; i++){
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_nom";
		tchiwg_h_chan_rate_nom[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_jesUp";
		tchiwg_h_chan_rate_jesUp[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_jesDown";
		tchiwg_h_chan_rate_jesDown[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_jerUp";
		tchiwg_h_chan_rate_jerUp[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_jerDown";
		tchiwg_h_chan_rate_jerDown[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_xsUp";
		tchiwg_h_chan_rate_xsUp[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
															
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_jes";
		tchiwg_h_chan_syserr_jes[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_jer";
		tchiwg_h_chan_syserr_jer[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_esf";
		tchiwg_h_chan_syserr_esf[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_scale";
		tchiwg_h_chan_syserr_scale[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_eleshape";
		tchiwg_h_chan_syserr_eleshape[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_jetshape";
		tchiwg_h_chan_syserr_jetshape[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_xs";
		tchiwg_h_chan_syserr_xs[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_syserr_lumi";
		tchiwg_h_chan_syserr_lumi[i] = new TH1D(histname.str().c_str(),histname.str().c_str(),40,287.5,1287.5);
	}

  TChain *mgtree_tchiwg;
  mgtree_tchiwg = new TChain("mgTree","mgTree");
  mgtree_tchiwg->Add("/uscms_data/d3/mengleis/Sep1/resTree_TChiWG.root");
  float phoEt_tchiwg_mg(0);
  float phoEta_tchiwg_mg(0);
  float lepPt_tchiwg_mg(0);
  float lepEta_tchiwg_mg(0);
  float sigMT_tchiwg_mg(0);
  float sigMET_tchiwg_mg(0);
  float HT_tchiwg_mg(0);
	float sigMETJESup_tchiwg_mg(0);
	float sigMETJESdo_tchiwg_mg(0);
	float sigMETJERup_tchiwg_mg(0);
	float sigMETJERdo_tchiwg_mg(0);
	float sigMTJESup_tchiwg_mg(0);
	float sigMTJESdo_tchiwg_mg(0);
	float sigMTJERup_tchiwg_mg(0);
	float sigMTJERdo_tchiwg_mg(0);
	float HTJESup_tchiwg_mg(0);
	float HTJESdo_tchiwg_mg(0);
	float gluinoMass_tchiwg_mg(0);
  float charginoMass_tchiwg_mg(0);
  float neutralinoMass_tchiwg_mg(0);
  mgtree_tchiwg->SetBranchAddress("phoEt",      &phoEt_tchiwg_mg);
  mgtree_tchiwg->SetBranchAddress("phoEta",     &phoEta_tchiwg_mg);
  mgtree_tchiwg->SetBranchAddress("lepPt",      &lepPt_tchiwg_mg);
  mgtree_tchiwg->SetBranchAddress("lepEta",     &lepEta_tchiwg_mg);
  mgtree_tchiwg->SetBranchAddress("sigMT",      &sigMT_tchiwg_mg);
  mgtree_tchiwg->SetBranchAddress("sigMET",     &sigMET_tchiwg_mg);
  mgtree_tchiwg->SetBranchAddress("HT",         &HT_tchiwg_mg);
	mgtree_tchiwg->SetBranchAddress("Mgluino",    &gluinoMass_tchiwg_mg);
  mgtree_tchiwg->SetBranchAddress("Mchargino",  &charginoMass_tchiwg_mg);
  mgtree_tchiwg->SetBranchAddress("Mneutralino",&neutralinoMass_tchiwg_mg);
	mgtree_tchiwg->SetBranchAddress("sigMETJESup",&sigMETJESup_tchiwg_mg);
	mgtree_tchiwg->SetBranchAddress("sigMETJESdo",&sigMETJESdo_tchiwg_mg);
	mgtree_tchiwg->SetBranchAddress("sigMETJERup",&sigMETJERup_tchiwg_mg);
	mgtree_tchiwg->SetBranchAddress("sigMETJERdo",&sigMETJERdo_tchiwg_mg);
	mgtree_tchiwg->SetBranchAddress("sigMTJESup", &sigMTJESup_tchiwg_mg);
	mgtree_tchiwg->SetBranchAddress("sigMTJESdo", &sigMTJESdo_tchiwg_mg);
	mgtree_tchiwg->SetBranchAddress("sigMTJERup", &sigMTJERup_tchiwg_mg);
	mgtree_tchiwg->SetBranchAddress("sigMTJERdo", &sigMTJERdo_tchiwg_mg);
	mgtree_tchiwg->SetBranchAddress("HTJESup",    &HTJESup_tchiwg_mg);
	mgtree_tchiwg->SetBranchAddress("HTJESdo",    &HTJESdo_tchiwg_mg);
		
	for(unsigned ievt(0); ievt < mgtree_tchiwg->GetEntries(); ievt++){
		mgtree_tchiwg->GetEntry(ievt);

		/** cut flow *****/
		if(phoEt_tchiwg_mg < 35 || lepPt_tchiwg_mg < 25)continue;
		if(fabs(phoEta_tchiwg_mg) > 1.4442 || fabs(lepEta_tchiwg_mg) > 2.5)continue;

		for(unsigned i(0); i < 10; i++){
			if( fabs(charginoMass_tchiwg_mg - (300+i*100) ) < 2){
				if(sigMT_tchiwg_mg < 100)p_tchiwg_MET_valid[i]->Fill( sigMET_tchiwg_mg );
				else{
				  p_tchiwg_MET_signal[i]->Fill( sigMET_tchiwg_mg );
				  p_tchiwg_HT_signal[i]->Fill( HT_tchiwg_mg );
				}
			}
		}

		if(sigMET_tchiwg_mg < 120 || sigMT_tchiwg_mg < 100)continue;

		int SigBinIndex(-1);
		SigBinIndex = findSignalBin(sigMET_tchiwg_mg, HT_tchiwg_mg, METbin1, METbin2);
		if(SigBinIndex >=0){
			tchiwg_h_chan_rate_nom[SigBinIndex]->Fill( charginoMass_tchiwg_mg, 1);
			tchiwg_h_chan_rate_xsUp[SigBinIndex]->Fill( charginoMass_tchiwg_mg, 1); 
		}
		int jesupBinIndex(-1);	
		jesupBinIndex = findSignalBin(sigMETJESup_tchiwg_mg, HTJESup_tchiwg_mg, METbin1, METbin2);
		if(jesupBinIndex >=0){
			tchiwg_h_chan_rate_jesUp[jesupBinIndex]->Fill( charginoMass_tchiwg_mg, 1); 
		}
		int jesdoBinIndex(-1);
		jesdoBinIndex = findSignalBin(sigMETJESdo_tchiwg_mg, HTJESdo_tchiwg_mg, METbin1, METbin2);
		if(jesdoBinIndex >=0){
			tchiwg_h_chan_rate_jesDown[jesdoBinIndex]->Fill( charginoMass_tchiwg_mg, 1);   
		}	
		int jerupBinIndex(-1);
		jerupBinIndex = findSignalBin(sigMETJERup_tchiwg_mg, HT_tchiwg_mg, METbin1, METbin2);
		if( jerupBinIndex >=0){
			tchiwg_h_chan_rate_jerUp[jerupBinIndex]->Fill(  charginoMass_tchiwg_mg, 1);
		}
		int jerdoBinIndex(-1);
		jerdoBinIndex = findSignalBin(sigMETJERdo_tchiwg_mg, HT_tchiwg_mg, METbin1, METbin2);	
		if(jerdoBinIndex >= 0){
			tchiwg_h_chan_rate_jerDown[jerdoBinIndex]->Fill( charginoMass_tchiwg_mg, 1);
		}  
	}


  TChain *egtree_tchiwg;
  egtree_tchiwg = new TChain("egTree","egTree");
  egtree_tchiwg->Add("/uscms_data/d3/mengleis/Sep1/resTree_TChiWG.root");
  float phoEt_tchiwg_eg(0);
  float phoEta_tchiwg_eg(0);
  float lepPt_tchiwg_eg(0);
  float lepEta_tchiwg_eg(0);
  float sigMT_tchiwg_eg(0);
  float sigMET_tchiwg_eg(0);
  float HT_tchiwg_eg(0);
	float sigMETJESup_tchiwg_eg(0);
	float sigMETJESdo_tchiwg_eg(0);
	float sigMETJERup_tchiwg_eg(0);
	float sigMETJERdo_tchiwg_eg(0);
	float sigMTJESup_tchiwg_eg(0);
	float sigMTJESdo_tchiwg_eg(0);
	float sigMTJERup_tchiwg_eg(0);
	float sigMTJERdo_tchiwg_eg(0);
	float HTJESup_tchiwg_eg(0);
	float HTJESdo_tchiwg_eg(0);
	float gluinoMass_tchiwg_eg(0);
  float charginoMass_tchiwg_eg(0);
  float neutralinoMass_tchiwg_eg(0);
  egtree_tchiwg->SetBranchAddress("phoEt",      &phoEt_tchiwg_eg);
  egtree_tchiwg->SetBranchAddress("phoEta",     &phoEta_tchiwg_eg);
  egtree_tchiwg->SetBranchAddress("lepPt",      &lepPt_tchiwg_eg);
  egtree_tchiwg->SetBranchAddress("lepEta",     &lepEta_tchiwg_eg);
  egtree_tchiwg->SetBranchAddress("sigMT",      &sigMT_tchiwg_eg);
  egtree_tchiwg->SetBranchAddress("sigMET",     &sigMET_tchiwg_eg);
  egtree_tchiwg->SetBranchAddress("HT",         &HT_tchiwg_eg);
	egtree_tchiwg->SetBranchAddress("Mgluino",    &gluinoMass_tchiwg_eg);
  egtree_tchiwg->SetBranchAddress("Mchargino",  &charginoMass_tchiwg_eg);
  egtree_tchiwg->SetBranchAddress("Mneutralino",&neutralinoMass_tchiwg_eg);
	egtree_tchiwg->SetBranchAddress("sigMETJESup",&sigMETJESup_tchiwg_eg);
	egtree_tchiwg->SetBranchAddress("sigMETJESdo",&sigMETJESdo_tchiwg_eg);
	egtree_tchiwg->SetBranchAddress("sigMETJERup",&sigMETJERup_tchiwg_eg);
	egtree_tchiwg->SetBranchAddress("sigMETJERdo",&sigMETJERdo_tchiwg_eg);
	egtree_tchiwg->SetBranchAddress("sigMTJESup", &sigMTJESup_tchiwg_eg);
	egtree_tchiwg->SetBranchAddress("sigMTJESdo", &sigMTJESdo_tchiwg_eg);
	egtree_tchiwg->SetBranchAddress("sigMTJERup", &sigMTJERup_tchiwg_eg);
	egtree_tchiwg->SetBranchAddress("sigMTJERdo", &sigMTJERdo_tchiwg_eg);
	egtree_tchiwg->SetBranchAddress("HTJESup",    &HTJESup_tchiwg_eg);
	egtree_tchiwg->SetBranchAddress("HTJESdo",    &HTJESdo_tchiwg_eg);
		
	for(unsigned ievt(0); ievt < egtree_tchiwg->GetEntries(); ievt++){
		egtree_tchiwg->GetEntry(ievt);

		/** cut flow *****/
		if(phoEt_tchiwg_eg < 35 || lepPt_tchiwg_eg < 25)continue;
		if(fabs(phoEta_tchiwg_eg) > 1.4442 || fabs(lepEta_tchiwg_eg) > 2.5)continue;
		if(sigMET_tchiwg_eg < 120 || sigMT_tchiwg_eg < 100)continue;

		int SigBinIndex(-1);
		SigBinIndex = findSignalBin(sigMET_tchiwg_eg, HT_tchiwg_eg, METbin1, METbin2)+9;
		if(SigBinIndex >=0){
			tchiwg_h_chan_rate_nom[SigBinIndex]->Fill( charginoMass_tchiwg_eg, 1);
			tchiwg_h_chan_rate_xsUp[SigBinIndex]->Fill( charginoMass_tchiwg_eg, 1); 
		}
		int jesupBinIndex(-1);	
		jesupBinIndex = findSignalBin(sigMETJESup_tchiwg_eg, HTJESup_tchiwg_eg, METbin1, METbin2)+9;
		if(jesupBinIndex >=0){
			tchiwg_h_chan_rate_jesUp[jesupBinIndex]->Fill( charginoMass_tchiwg_eg, 1); 
		}
		int jesdoBinIndex(-1);
		jesdoBinIndex = findSignalBin(sigMETJESdo_tchiwg_eg, HTJESdo_tchiwg_eg, METbin1, METbin2)+9;
		if(jesdoBinIndex >=0){
			tchiwg_h_chan_rate_jesDown[jesdoBinIndex]->Fill( charginoMass_tchiwg_eg, 1);   
		}	
		int jerupBinIndex(-1);
		jerupBinIndex = findSignalBin(sigMETJERup_tchiwg_eg, HT_tchiwg_eg, METbin1, METbin2)+9;
		if( jerupBinIndex >=0){
			tchiwg_h_chan_rate_jerUp[jerupBinIndex]->Fill(  charginoMass_tchiwg_eg, 1);
		}
		int jerdoBinIndex(-1);
		jerdoBinIndex = findSignalBin(sigMETJERdo_tchiwg_eg, HT_tchiwg_eg, METbin1, METbin2)+9;	
		if(jerdoBinIndex >= 0){
			tchiwg_h_chan_rate_jerDown[jerdoBinIndex]->Fill( charginoMass_tchiwg_eg, 1);
		}  
	}

	for(unsigned ih(0); ih < 18; ih++){
		for(unsigned i(1); i < tchiwg_h_chan_rate_nom[ih]->GetXaxis()->GetNbins() + 1; i++){
				
			if(p_TChiWGMASS->GetBinContent(i) <= 0){
				tchiwg_h_chan_rate_nom[ih]->SetBinContent(i, -1); 
				tchiwg_h_chan_rate_jesUp[ih]->SetBinContent(i,-1);
				tchiwg_h_chan_rate_jesDown[ih]->SetBinContent(i, -1);
				tchiwg_h_chan_rate_jerUp[ih]->SetBinContent(i, -1);
				tchiwg_h_chan_rate_jerDown[ih]->SetBinContent(i,-1);
				tchiwg_h_chan_rate_xsUp[ih]->SetBinContent(i, -1);
						
				tchiwg_h_chan_syserr_jes[ih]->SetBinContent(i, -1); 
				tchiwg_h_chan_syserr_jer[ih]->SetBinContent(i, -1); 
				tchiwg_h_chan_syserr_esf[ih]->SetBinContent(i, -1); 
				tchiwg_h_chan_syserr_scale[ih]->SetBinContent(i, -1); 
				tchiwg_h_chan_syserr_eleshape[ih]->SetBinContent(i, -1); 
				tchiwg_h_chan_syserr_jetshape[ih]->SetBinContent(i, -1); 
				tchiwg_h_chan_syserr_xs[ih]->SetBinContent(i, -1); 
				tchiwg_h_chan_syserr_lumi[ih]->SetBinContent(i, -1); 
			}
			else{ 
				float noe = p_TChiWGMASS->GetBinContent(i);
				float sparticleMass = p_TChiWGMASS->GetXaxis()->GetBinCenter(i);
				float crosssection = p_crosssection_tchiwg->GetBinContent( p_crosssection_tchiwg->FindBin(sparticleMass) );
				float crosssectionUp = (crosssection+p_crosssection_tchiwg->GetBinError( p_crosssection_tchiwg->FindBin(sparticleMass) ) ); 

				tchiwg_h_chan_rate_nom[ih]->SetBinContent(i, tchiwg_h_chan_rate_nom[ih]->GetBinContent(i)*35.8*crosssection/noe); 
				tchiwg_h_chan_rate_jesUp[ih]->SetBinContent(i, tchiwg_h_chan_rate_jesUp[ih]->GetBinContent(i)*35.8*crosssection/noe); 
				tchiwg_h_chan_rate_jesDown[ih]->SetBinContent(i, tchiwg_h_chan_rate_jesDown[ih]->GetBinContent(i)*35.8*crosssection/noe);
				tchiwg_h_chan_rate_jerUp[ih]->SetBinContent(i, tchiwg_h_chan_rate_jerUp[ih]->GetBinContent(i)*35.8*crosssection/noe);
				tchiwg_h_chan_rate_jerDown[ih]->SetBinContent(i, tchiwg_h_chan_rate_jerDown[ih]->GetBinContent(i)*35.8*crosssection/noe);
				tchiwg_h_chan_rate_xsUp[ih]->SetBinContent(i, tchiwg_h_chan_rate_xsUp[ih]->GetBinContent(i)*35.8*crosssectionUp/noe);
								
				tchiwg_h_chan_syserr_jes[ih]->SetBinContent(i, max( fabs(tchiwg_h_chan_rate_jesUp[ih]->GetBinContent(i)-tchiwg_h_chan_rate_nom[ih]->GetBinContent(i)), fabs(tchiwg_h_chan_rate_jesDown[ih]->GetBinContent(i)-tchiwg_h_chan_rate_nom[ih]->GetBinContent(i)) ) ); 
				tchiwg_h_chan_syserr_jer[ih]->SetBinContent(i, max( fabs(tchiwg_h_chan_rate_jerUp[ih]->GetBinContent(i)-tchiwg_h_chan_rate_nom[ih]->GetBinContent(i)), fabs(tchiwg_h_chan_rate_jerDown[ih]->GetBinContent(i)-tchiwg_h_chan_rate_nom[ih]->GetBinContent(i)) ) ); 
				tchiwg_h_chan_syserr_esf[ih]->SetBinContent(i, -1);
				tchiwg_h_chan_syserr_scale[ih]->SetBinContent(i, -1);
				tchiwg_h_chan_syserr_eleshape[ih]->SetBinContent(i, -1);
				tchiwg_h_chan_syserr_jetshape[ih]->SetBinContent(i, -1);
				tchiwg_h_chan_syserr_xs[ih]->SetBinContent(i, fabs( tchiwg_h_chan_rate_xsUp[ih]->GetBinContent(i)-tchiwg_h_chan_rate_nom[ih]->GetBinContent(i)) );
				tchiwg_h_chan_syserr_lumi[ih]->SetBinContent(i, 0.026*tchiwg_h_chan_rate_nom[ih]->GetBinContent(i));   
			}
		} 
	} 

		for(unsigned i(1); i < p_TChiWGMASS->GetXaxis()->GetNbins() + 1; i++){
				
			if(p_TChiWGMASS->GetBinContent(i) > 0){
				float noe = p_TChiWGMASS->GetBinContent(i);
				float sparticleMass = p_TChiWGMASS->GetXaxis()->GetBinCenter(i);
				float crosssection = p_crosssection_tchiwg->GetBinContent( p_crosssection_tchiwg->FindBin(sparticleMass) );

				for(unsigned j(0); j < 10; j++){
					if( fabs(sparticleMass - (300+j*100) ) < 2){
						p_tchiwg_MET_valid[j]->Scale(35.8*crosssection/noe);
						p_tchiwg_MET_signal[j]->Scale(35.8*crosssection/noe);
						p_tchiwg_HT_signal[j]->Scale(35.8*crosssection/noe);
					}
				}
			}
		} 

	outputfile_tchiwg->Write();
	outputfile_tchiwg->Close();



  xSecFile.Close();
}


