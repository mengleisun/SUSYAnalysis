#include "analysis_PU.C"
#include "analysis_SUSY.h"

void analysis_T5WG(){//main  

	SetSignalConfig();
	binning Bin(NBIN, METbin1, METbin2, HTbin1, HTbin2, PHOETbin);
	esfScaleFactor  objectESF;

  gSystem->Load("/uscms/homes/t/tmishra/work/CMSSW_10_2_22/src/SUSYAnalysis/lib/libAnaClasses.so");

	TFile xSecFile("../cross/susyCrossSection.root");
	TH1D *p_crosssection_susy   = (TH1D*)xSecFile.Get("p_gluinoxSec");

	TChain *datachain = new TChain("signalTree");
	datachain->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_mgsignal_MuonEG_FullEcal.root");
	TH1D *p_PU_data = new TH1D("p_PU_data",";N_{vtx};",100,0,100); 
  datachain->Draw("nVertex >> p_PU_data");
	p_PU_data->Scale(1.0/p_PU_data->Integral(1,101));

	std::ostringstream histname;
	//**************   T5WG  ***************************//
  //TFile *file_susy = TFile::Open("/uscms_data/d3/mengleis/Sep1/resTree_T5WG.root");
  TFile *file_susy = TFile::Open("/uscms_data/d3/mengleis/FullStatusOct/resTree_T5WG.root");
  TTree *tree_susy = (TTree*)file_susy->Get("SUSYtree");
	float Mgluino_susy(0);
  float Mchargino_susy(0);
  float Mneutralino_susy(0);
	float Mass1_susy(0);
	float Mass2_susy(0);
	int   nVertex(0);
	tree_susy->SetBranchAddress("MsGsQ",      &Mgluino_susy);  
  tree_susy->SetBranchAddress("Mchargino",  &Mchargino_susy);
  tree_susy->SetBranchAddress("Mneutralino",&Mneutralino_susy);
  tree_susy->SetBranchAddress("Mass1",      &Mass1_susy);
  tree_susy->SetBranchAddress("Mass2",      &Mass2_susy);
  tree_susy->SetBranchAddress("nVertex",    &nVertex);

	init_histo("T5WG","test_T5WG.root",2,27, 775.0, 2125.0, 420, 2.5, 2102.5);
	for(unsigned ievt(0); ievt < tree_susy->GetEntries(); ievt++){
		tree_susy->GetEntry(ievt);
    if (ievt%1000000==0) std::cout << " -- Processing event " << ievt << std::endl;
	
	//	double NLSPMass(0);
  //  if(Mchargino_susy >0)NLSPMass = Mchargino_susy;
	//	else if(Mneutralino_susy >0)NLSPMass = Mneutralino_susy;
	//	p_SUSYMass->Fill(Mgluino_susy, NLSPMass);

		p_SUSYMass->Fill(Mass1_susy, Mass2_susy);
		if(nVertex < 20)p_lowPU_susy_all->Fill(Mass1_susy, Mass2_susy, 1);
		else p_highPU_susy_all->Fill(Mass1_susy, Mass2_susy, 1);
	}

  TChain *mgtree_susy;
  mgtree_susy = new TChain("mgTree","mgTree");
  //mgtree_susy->Add("/uscms_data/d3/mengleis/Sep1/resTree_T5WG.root");
  mgtree_susy->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_T5WG.root");
  float phoEt_susy_mg(0);
  float phoEta_susy_mg(0);
  float lepPt_susy_mg(0);
  float lepEta_susy_mg(0);
  float sigMT_susy_mg(0);
  float sigMET_susy_mg(0);
  float HT_susy_mg(0);
	int   nVertex_susy_mg(0);
	float sigMETJESup_susy_mg(0);
	float sigMETJESdo_susy_mg(0);
	float sigMETJERup_susy_mg(0);
	float sigMETJERdo_susy_mg(0);
	float sigMTJESup_susy_mg(0);
	float sigMTJESdo_susy_mg(0);
	float sigMTJERup_susy_mg(0);
	float sigMTJERdo_susy_mg(0);
	float HTJESup_susy_mg(0);
	float HTJESdo_susy_mg(0);
	float gluinoMass_susy_mg(0);
  float charginoMass_susy_mg(0);
  float neutralinoMass_susy_mg(0);
  float Mass1_susy_mg(0);
  float Mass2_susy_mg(0);
  mgtree_susy->SetBranchAddress("phoEt",      &phoEt_susy_mg);
  mgtree_susy->SetBranchAddress("phoEta",     &phoEta_susy_mg);
  mgtree_susy->SetBranchAddress("lepPt",      &lepPt_susy_mg);
  mgtree_susy->SetBranchAddress("lepEta",     &lepEta_susy_mg);
  mgtree_susy->SetBranchAddress("sigMT",      &sigMT_susy_mg);
  mgtree_susy->SetBranchAddress("sigMET",     &sigMET_susy_mg);
  mgtree_susy->SetBranchAddress("HT",         &HT_susy_mg);
	mgtree_susy->SetBranchAddress("nVertex",    &nVertex_susy_mg);
	mgtree_susy->SetBranchAddress("MsGsQ",      &gluinoMass_susy_mg);
  mgtree_susy->SetBranchAddress("Mchargino",  &charginoMass_susy_mg);
  mgtree_susy->SetBranchAddress("Mneutralino",&neutralinoMass_susy_mg);
	mgtree_susy->SetBranchAddress("Mass1",      &Mass1_susy_mg);
	mgtree_susy->SetBranchAddress("Mass2",      &Mass2_susy_mg);
	mgtree_susy->SetBranchAddress("sigMETJESup",&sigMETJESup_susy_mg);
	mgtree_susy->SetBranchAddress("sigMETJESdo",&sigMETJESdo_susy_mg);
	mgtree_susy->SetBranchAddress("sigMETJERup",&sigMETJERup_susy_mg);
	mgtree_susy->SetBranchAddress("sigMETJERdo",&sigMETJERdo_susy_mg);
	mgtree_susy->SetBranchAddress("sigMTJESup", &sigMTJESup_susy_mg);
	mgtree_susy->SetBranchAddress("sigMTJESdo", &sigMTJESdo_susy_mg);
	mgtree_susy->SetBranchAddress("sigMTJERup", &sigMTJERup_susy_mg);
	mgtree_susy->SetBranchAddress("sigMTJERdo", &sigMTJERdo_susy_mg);
	mgtree_susy->SetBranchAddress("HTJESup",    &HTJESup_susy_mg);
	mgtree_susy->SetBranchAddress("HTJESdo",    &HTJESdo_susy_mg);

	for(unsigned ievt(0); ievt < mgtree_susy->GetEntries(); ievt++){
		mgtree_susy->GetEntry(ievt);
	
		/** cut flow *****/
		if(phoEt_susy_mg < 35 || lepPt_susy_mg < 25)continue;
		if(fabs(phoEta_susy_mg) > 1.4442 || fabs(lepEta_susy_mg) > 2.5)continue;
    p_SUSYselect->Fill(Mass1_susy_mg, Mass2_susy_mg);

		double scalefactor(0);
		double scalefactorup(0);
		scalefactor = objectESF.getFastMuonESF(lepPt_susy_mg,lepEta_susy_mg)*objectESF.getPhotonESF(phoEt_susy_mg, phoEta_susy_mg)*objectESF.getFastMuonEGTRGESF(phoEt_susy_mg, lepPt_susy_mg);
		double s_mu_error = objectESF.getFastMuonESFError(lepPt_susy_mg,lepEta_susy_mg)*objectESF.getPhotonESF(phoEt_susy_mg, phoEta_susy_mg)*objectESF.getFastMuonEGTRGESF(phoEt_susy_mg, lepPt_susy_mg);
    double s_pho_error = objectESF.getPhotonESFError(phoEt_susy_mg, phoEta_susy_mg)*objectESF.getFastMuonESF(lepPt_susy_mg,lepEta_susy_mg)*objectESF.getFastMuonEGTRGESF(phoEt_susy_mg, lepPt_susy_mg);
    double s_trg_error = objectESF.getFastMuonEGTRGESFError(phoEt_susy_mg, lepPt_susy_mg)*objectESF.getFastMuonESF(lepPt_susy_mg,lepEta_susy_mg)*objectESF.getPhotonESF(phoEt_susy_mg, phoEta_susy_mg);
		double s_error = sqrt(pow(s_mu_error,2) + pow(s_pho_error,2) + pow(s_trg_error,2));
		scalefactorup = scalefactor + s_error; 

		if( fabs(Mass1_susy_mg - 1700) < 10 && fabs(Mass2_susy_mg - 1000) < 5){
			if(sigMT_susy_mg > 100){
				p_susy_MET_signal_mg->Fill(sigMET_susy_mg);
				if(sigMET_susy_mg > 120)p_susy_HT_signal_mg->Fill(HT_susy_mg);
				if(sigMET_susy_mg > 120)p_susy_PhoEt_signal_mg->Fill(phoEt_susy_mg);
			}
		}

		if(sigMET_susy_mg < 120 || sigMT_susy_mg < 100)continue;

		if(charginoMass_susy_mg > 0 && neutralinoMass_susy_mg > 0){
 
			int SigBinIndex(-1);
			SigBinIndex = Bin.findSignalBin(sigMET_susy_mg, HT_susy_mg, phoEt_susy_mg);
			if(SigBinIndex >=0){
				susy_wg_chan_rate_nom[SigBinIndex]->Fill( Mass1_susy_mg, Mass2_susy_mg, scalefactor);
				susy_wg_chan_rate_xsUp[SigBinIndex]->Fill( Mass1_susy_mg, Mass2_susy_mg, scalefactor); 
				susy_wg_chan_rate_esfUp[SigBinIndex]->Fill( Mass1_susy_mg, Mass2_susy_mg,scalefactorup); 
	
				if(nVertex_susy_mg < 20)p_lowPU_susy_pass->Fill( Mass1_susy_mg, Mass2_susy_mg, 1);
				else p_highPU_susy_pass->Fill( Mass1_susy_mg, Mass2_susy_mg, 1);
			}
			int jesupBinIndex(-1);	
			jesupBinIndex = Bin.findSignalBin(sigMETJESup_susy_mg, HTJESup_susy_mg, phoEt_susy_mg); 
			if(jesupBinIndex >=0){
				susy_wg_chan_rate_jesUp[jesupBinIndex]->Fill( Mass1_susy_mg, Mass2_susy_mg, scalefactor); 
			}
			int jesdoBinIndex(-1);
			jesdoBinIndex = Bin.findSignalBin(sigMETJESdo_susy_mg, HTJESdo_susy_mg, phoEt_susy_mg);
			if(jesdoBinIndex >=0){
				susy_wg_chan_rate_jesDown[jesdoBinIndex]->Fill( Mass1_susy_mg, Mass2_susy_mg, scalefactor);   
			}	
			int jerupBinIndex(-1);
			jerupBinIndex = Bin.findSignalBin(sigMETJERup_susy_mg, HT_susy_mg, phoEt_susy_mg); 
			if( jerupBinIndex >=0){
				susy_wg_chan_rate_jerUp[jerupBinIndex]->Fill( Mass1_susy_mg, Mass2_susy_mg, scalefactor);
			}
			int jerdoBinIndex(-1);
			jerdoBinIndex = Bin.findSignalBin(sigMETJERdo_susy_mg, HT_susy_mg, phoEt_susy_mg);
			if(jerdoBinIndex >= 0){
				susy_wg_chan_rate_jerDown[jerdoBinIndex]->Fill( Mass1_susy_mg, Mass2_susy_mg, scalefactor);
			} 
		} 
		else if(charginoMass_susy_mg <= 0 && neutralinoMass_susy_mg > 0){
 
			int SigBinIndex(-1);
			SigBinIndex = Bin.findSignalBin(sigMET_susy_mg, HT_susy_mg, phoEt_susy_mg);
			if(SigBinIndex >=0){
				susy_gg_chan_rate_nom[SigBinIndex]->Fill( Mass1_susy_mg, Mass2_susy_mg, scalefactor);
				susy_gg_chan_rate_xsUp[SigBinIndex]->Fill( Mass1_susy_mg, Mass2_susy_mg, scalefactor); 
				susy_gg_chan_rate_esfUp[SigBinIndex]->Fill( Mass1_susy_mg, Mass2_susy_mg,scalefactorup); 
	
				if(nVertex_susy_mg < 20)p_lowPU_susy_pass->Fill( Mass1_susy_mg, Mass2_susy_mg, 1);
				else p_highPU_susy_pass->Fill( Mass1_susy_mg, Mass2_susy_mg, 1);
			}
			int jesupBinIndex(-1);	
			jesupBinIndex = Bin.findSignalBin(sigMETJESup_susy_mg, HTJESup_susy_mg, phoEt_susy_mg); 
			if(jesupBinIndex >=0){
				susy_gg_chan_rate_jesUp[jesupBinIndex]->Fill( Mass1_susy_mg, Mass2_susy_mg, scalefactor); 
			}
			int jesdoBinIndex(-1);
			jesdoBinIndex = Bin.findSignalBin(sigMETJESdo_susy_mg, HTJESdo_susy_mg, phoEt_susy_mg);
			if(jesdoBinIndex >=0){
				susy_gg_chan_rate_jesDown[jesdoBinIndex]->Fill( Mass1_susy_mg, Mass2_susy_mg, scalefactor);   
			}	
			int jerupBinIndex(-1);
			jerupBinIndex = Bin.findSignalBin(sigMETJERup_susy_mg, HT_susy_mg, phoEt_susy_mg); 
			if( jerupBinIndex >=0){
				susy_gg_chan_rate_jerUp[jerupBinIndex]->Fill( Mass1_susy_mg, Mass2_susy_mg, scalefactor);
			}
			int jerdoBinIndex(-1);
			jerdoBinIndex = Bin.findSignalBin(sigMETJERdo_susy_mg, HT_susy_mg, phoEt_susy_mg);
			if(jerdoBinIndex >= 0){
				susy_gg_chan_rate_jerDown[jerdoBinIndex]->Fill( Mass1_susy_mg, Mass2_susy_mg, scalefactor);
			} 
		} 
	}

  TChain *egtree_susy;
  egtree_susy = new TChain("egTree","egTree");
  egtree_susy->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_T5WG.root");
  float phoEt_susy_eg(0);
  float phoEta_susy_eg(0);
  float lepPt_susy_eg(0);
  float lepEta_susy_eg(0);
  float sigMT_susy_eg(0);
  float sigMET_susy_eg(0);
  float HT_susy_eg(0);
  int   nVertex_susy_eg(0); 
	float sigMETJESup_susy_eg(0);
	float sigMETJESdo_susy_eg(0);
	float sigMETJERup_susy_eg(0);
	float sigMETJERdo_susy_eg(0);
	float sigMTJESup_susy_eg(0);
	float sigMTJESdo_susy_eg(0);
	float sigMTJERup_susy_eg(0);
	float sigMTJERdo_susy_eg(0);
	float HTJESup_susy_eg(0);
	float HTJESdo_susy_eg(0);
	float gluinoMass_susy_eg(0);
  float charginoMass_susy_eg(0);
  float neutralinoMass_susy_eg(0);
	float Mass1_susy_eg(0);
  float Mass2_susy_eg(0);
  egtree_susy->SetBranchAddress("phoEt",      &phoEt_susy_eg);
  egtree_susy->SetBranchAddress("phoEta",     &phoEta_susy_eg);
  egtree_susy->SetBranchAddress("lepPt",      &lepPt_susy_eg);
  egtree_susy->SetBranchAddress("lepEta",     &lepEta_susy_eg);
  egtree_susy->SetBranchAddress("sigMT",      &sigMT_susy_eg);
  egtree_susy->SetBranchAddress("sigMET",     &sigMET_susy_eg);
	egtree_susy->SetBranchAddress("nVertex",    &nVertex_susy_eg);
  egtree_susy->SetBranchAddress("HT",         &HT_susy_eg);
	egtree_susy->SetBranchAddress("MsGsQ",      &gluinoMass_susy_eg);
  egtree_susy->SetBranchAddress("Mchargino",  &charginoMass_susy_eg);
  egtree_susy->SetBranchAddress("Mneutralino",&neutralinoMass_susy_eg);
	egtree_susy->SetBranchAddress("Mass1",      &Mass1_susy_eg);
	egtree_susy->SetBranchAddress("Mass2",      &Mass2_susy_eg);
	egtree_susy->SetBranchAddress("sigMETJESup",&sigMETJESup_susy_eg);
	egtree_susy->SetBranchAddress("sigMETJESdo",&sigMETJESdo_susy_eg);
	egtree_susy->SetBranchAddress("sigMETJERup",&sigMETJERup_susy_eg);
	egtree_susy->SetBranchAddress("sigMETJERdo",&sigMETJERdo_susy_eg);
	egtree_susy->SetBranchAddress("sigMTJESup", &sigMTJESup_susy_eg);
	egtree_susy->SetBranchAddress("sigMTJESdo", &sigMTJESdo_susy_eg);
	egtree_susy->SetBranchAddress("sigMTJERup", &sigMTJERup_susy_eg);
	egtree_susy->SetBranchAddress("sigMTJERdo", &sigMTJERdo_susy_eg);
	egtree_susy->SetBranchAddress("HTJESup",    &HTJESup_susy_eg);
	egtree_susy->SetBranchAddress("HTJESdo",    &HTJESdo_susy_eg);

	for(unsigned ievt(0); ievt < egtree_susy->GetEntries(); ievt++){
		egtree_susy->GetEntry(ievt);

		/** cut flow *****/
		if(phoEt_susy_eg < 35 || lepPt_susy_eg < 25)continue;
		if(fabs(phoEta_susy_eg) > 1.4442 || fabs(lepEta_susy_eg) > 2.5)continue;

		double scalefactor(1);
		double scalefactorup(0);
		scalefactor = objectESF.getFastElectronESF(lepPt_susy_eg,lepEta_susy_eg)*objectESF.getPhotonESF(phoEt_susy_eg,phoEta_susy_eg)*objectESF.getFastegPhotonTRGESF(phoEt_susy_eg,phoEta_susy_eg)*objectESF.getFastElectronTRGESF(lepPt_susy_eg,lepEta_susy_eg);
		double s_ele_error = objectESF.getFastElectronESFError(lepPt_susy_eg,lepEta_susy_eg)*objectESF.getPhotonESF(phoEt_susy_eg,phoEta_susy_eg)*objectESF.getFastegPhotonTRGESF(phoEt_susy_eg,phoEta_susy_eg)*objectESF.getFastElectronTRGESF(lepPt_susy_eg,lepEta_susy_eg);
		double s_pho_error = objectESF.getPhotonESFError(phoEt_susy_eg,phoEta_susy_eg)*objectESF.getFastElectronESF(lepPt_susy_eg,lepEta_susy_eg)*objectESF.getFastegPhotonTRGESF(phoEt_susy_eg,phoEta_susy_eg)*objectESF.getFastElectronTRGESF(lepPt_susy_eg,lepEta_susy_eg);
		double s_eletrg_error = objectESF.getFastElectronTRGESFError(lepPt_susy_eg,lepEta_susy_eg)*objectESF.getFastElectronESF(lepPt_susy_eg,lepEta_susy_eg)*objectESF.getPhotonESF(phoEt_susy_eg,phoEta_susy_eg)*objectESF.getFastegPhotonTRGESF(phoEt_susy_eg,phoEta_susy_eg);
		double s_photrg_error = objectESF.getFastegPhotonTRGESFError(phoEt_susy_eg,phoEta_susy_eg)*objectESF.getFastElectronESF(lepPt_susy_eg,lepEta_susy_eg)*objectESF.getPhotonESF(phoEt_susy_eg,phoEta_susy_eg)*objectESF.getFastElectronTRGESF(lepPt_susy_eg,lepEta_susy_eg);
		double s_error = sqrt(pow(s_ele_error,2) + pow(s_pho_error,2) + pow(s_eletrg_error,2) + pow(s_photrg_error, 2)); 
		scalefactorup = scalefactor + s_error; 

		if( fabs(Mass1_susy_eg - 1700) < 10 && fabs(Mass2_susy_eg - 1000) < 5){
			if(sigMT_susy_eg > 100){
				p_susy_MET_signal_eg->Fill(sigMET_susy_eg);
				if(sigMET_susy_eg > 120)p_susy_HT_signal_eg->Fill(HT_susy_eg);
				if(sigMET_susy_eg > 120)p_susy_PhoEt_signal_eg->Fill(phoEt_susy_eg);
			}
		}

		if(sigMET_susy_eg < 120 || sigMT_susy_eg < 100)continue;

		if(charginoMass_susy_eg > 0 && neutralinoMass_susy_eg > 0){
			int SigBinIndex(-1);
			SigBinIndex = Bin.findSignalBin(sigMET_susy_eg, HT_susy_eg, phoEt_susy_eg) + NBIN;
			if(SigBinIndex >=0){
				susy_wg_chan_rate_nom[SigBinIndex]->Fill( Mass1_susy_eg, Mass2_susy_eg, scalefactor);
				susy_wg_chan_rate_xsUp[SigBinIndex]->Fill( Mass1_susy_eg, Mass2_susy_eg, scalefactor); 
				susy_wg_chan_rate_esfUp[SigBinIndex]->Fill( Mass1_susy_eg, Mass2_susy_eg,scalefactorup); 
	
				if(nVertex_susy_eg < 20)p_lowPU_susy_pass->Fill( Mass1_susy_eg, Mass2_susy_eg, 1);
				else p_highPU_susy_pass->Fill( Mass1_susy_eg, Mass2_susy_eg, 1);
			}
			int jesupBinIndex(-1);	
			jesupBinIndex = Bin.findSignalBin(sigMETJESup_susy_eg, HTJESup_susy_eg, phoEt_susy_eg) + NBIN;
			if(jesupBinIndex >=0){
				susy_wg_chan_rate_jesUp[jesupBinIndex]->Fill( Mass1_susy_eg, Mass2_susy_eg, scalefactor); 
			}
			int jesdoBinIndex(-1);
			jesdoBinIndex = Bin.findSignalBin(sigMETJESdo_susy_eg, HTJESdo_susy_eg, phoEt_susy_eg) + NBIN;
			if(jesdoBinIndex >=0){
				susy_wg_chan_rate_jesDown[jesdoBinIndex]->Fill( Mass1_susy_eg, Mass2_susy_eg, scalefactor);   
			}	
			int jerupBinIndex(-1);
			jerupBinIndex = Bin.findSignalBin(sigMETJERup_susy_eg, HT_susy_eg, phoEt_susy_eg) + NBIN;
			if( jerupBinIndex >=0){
				susy_wg_chan_rate_jerUp[jerupBinIndex]->Fill( Mass1_susy_eg, Mass2_susy_eg, scalefactor);
			}
			int jerdoBinIndex(-1);
			jerdoBinIndex = Bin.findSignalBin(sigMETJERdo_susy_eg, HT_susy_eg, phoEt_susy_eg) + NBIN;	
			if(jerdoBinIndex >= 0){
				susy_wg_chan_rate_jerDown[jerdoBinIndex]->Fill( Mass1_susy_eg, Mass2_susy_eg, scalefactor);
			} 
		} 
		else if(charginoMass_susy_eg <= 0 && neutralinoMass_susy_eg > 0){
			int SigBinIndex(-1);
			SigBinIndex = Bin.findSignalBin(sigMET_susy_eg, HT_susy_eg, phoEt_susy_eg) + NBIN;
			if(SigBinIndex >=0){
				susy_gg_chan_rate_nom[SigBinIndex]->Fill( Mass1_susy_eg, Mass2_susy_eg, scalefactor);
				susy_gg_chan_rate_xsUp[SigBinIndex]->Fill( Mass1_susy_eg, Mass2_susy_eg, scalefactor); 
				susy_gg_chan_rate_esfUp[SigBinIndex]->Fill( Mass1_susy_eg, Mass2_susy_eg,scalefactorup); 
	
				if(nVertex_susy_eg < 20)p_lowPU_susy_pass->Fill( Mass1_susy_eg, Mass2_susy_eg, 1);
				else p_highPU_susy_pass->Fill( Mass1_susy_eg, Mass2_susy_eg, 1);
			}
			int jesupBinIndex(-1);	
			jesupBinIndex = Bin.findSignalBin(sigMETJESup_susy_eg, HTJESup_susy_eg, phoEt_susy_eg) + NBIN;
			if(jesupBinIndex >=0){
				susy_gg_chan_rate_jesUp[jesupBinIndex]->Fill( Mass1_susy_eg, Mass2_susy_eg, scalefactor); 
			}
			int jesdoBinIndex(-1);
			jesdoBinIndex = Bin.findSignalBin(sigMETJESdo_susy_eg, HTJESdo_susy_eg, phoEt_susy_eg) + NBIN;
			if(jesdoBinIndex >=0){
				susy_gg_chan_rate_jesDown[jesdoBinIndex]->Fill( Mass1_susy_eg, Mass2_susy_eg, scalefactor);   
			}	
			int jerupBinIndex(-1);
			jerupBinIndex = Bin.findSignalBin(sigMETJERup_susy_eg, HT_susy_eg, phoEt_susy_eg) + NBIN;
			if( jerupBinIndex >=0){
				susy_gg_chan_rate_jerUp[jerupBinIndex]->Fill( Mass1_susy_eg, Mass2_susy_eg, scalefactor);
			}
			int jerdoBinIndex(-1);
			jerdoBinIndex = Bin.findSignalBin(sigMETJERdo_susy_eg, HT_susy_eg, phoEt_susy_eg) + NBIN;	
			if(jerdoBinIndex >= 0){
				susy_gg_chan_rate_jerDown[jerdoBinIndex]->Fill( Mass1_susy_eg, Mass2_susy_eg, scalefactor);
			} 
		} 
	}

	for(int ih(0); ih < NBIN*2; ih++){
		for(int i(1); i < susy_wg_chan_rate_nom[ih]->GetXaxis()->GetNbins() + 1; i++){
			for(int j(1); j < susy_wg_chan_rate_nom[ih]->GetYaxis()->GetNbins() + 1; j++){
				if(p_SUSYMass->GetBinContent(i,j) < 1000){
					p_SUSYselect->SetBinContent(i,j,-1);
					p_SUSYMass->SetBinContent(i,j,-1);
				}
	
				if(p_SUSYMass->GetBinContent(i,j) <= 0){
					susy_wg_chan_rate_nom[ih]->SetBinContent(i,j, -1); 
					susy_wg_chan_rate_jesUp[ih]->SetBinContent(i,j,-1);
					susy_wg_chan_rate_jesDown[ih]->SetBinContent(i,j, -1);
					susy_wg_chan_rate_jerUp[ih]->SetBinContent(i,j, -1);
					susy_wg_chan_rate_jerDown[ih]->SetBinContent(i,j,-1);
					susy_wg_chan_rate_xsUp[ih]->SetBinContent(i,j, -1);
							
					susy_wg_syserr_PU->SetBinContent(i,j, -1);
					susy_wg_chan_syserr_jes[ih]->SetBinContent(i,j, -1); 
					susy_wg_chan_syserr_jer[ih]->SetBinContent(i,j, -1); 
					susy_wg_chan_syserr_esf[ih]->SetBinContent(i,j, -1); 
					susy_wg_chan_syserr_scale[ih]->SetBinContent(i,j, -1); 
					susy_wg_chan_syserr_eleshape[ih]->SetBinContent(i,j, -1); 
					susy_wg_chan_syserr_jetshape[ih]->SetBinContent(i,j, -1); 
					susy_wg_chan_syserr_qcdshape[ih]->SetBinContent(i,j, -1); 
					susy_wg_chan_syserr_xs[ih]->SetBinContent(i,j, -1); 
					susy_wg_chan_syserr_lumi[ih]->SetBinContent(i,j, -1); 
					susy_wg_chan_syserr_isr[ih]->SetBinContent(i,j, -1); 

					susy_gg_chan_rate_nom[ih]->SetBinContent(i,j, -1); 
					susy_gg_chan_rate_jesUp[ih]->SetBinContent(i,j,-1);
					susy_gg_chan_rate_jesDown[ih]->SetBinContent(i,j, -1);
					susy_gg_chan_rate_jerUp[ih]->SetBinContent(i,j, -1);
					susy_gg_chan_rate_jerDown[ih]->SetBinContent(i,j,-1);
					susy_gg_chan_rate_xsUp[ih]->SetBinContent(i,j, -1);
							
					susy_gg_syserr_PU->SetBinContent(i,j, -1);
					susy_gg_chan_syserr_jes[ih]->SetBinContent(i,j, -1); 
					susy_gg_chan_syserr_jer[ih]->SetBinContent(i,j, -1); 
					susy_gg_chan_syserr_esf[ih]->SetBinContent(i,j, -1); 
					susy_gg_chan_syserr_scale[ih]->SetBinContent(i,j, -1); 
					susy_gg_chan_syserr_eleshape[ih]->SetBinContent(i,j, -1); 
					susy_gg_chan_syserr_jetshape[ih]->SetBinContent(i,j, -1); 
					susy_gg_chan_syserr_qcdshape[ih]->SetBinContent(i,j, -1); 
					susy_gg_chan_syserr_xs[ih]->SetBinContent(i,j, -1); 
					susy_gg_chan_syserr_lumi[ih]->SetBinContent(i,j, -1); 
					susy_gg_chan_syserr_isr[ih]->SetBinContent(i,j, -1); 
				}

				else{
					float noe = p_SUSYMass->GetBinContent(i,j);
					float sparticleMass = p_SUSYMass->GetXaxis()->GetBinCenter(i);
					float crosssection = p_crosssection_susy->GetBinContent( p_crosssection_susy->FindBin(sparticleMass) );
					float crosssectionUp = (crosssection+ p_crosssection_susy->GetBinError( p_crosssection_susy->FindBin(sparticleMass) ) );

					susy_wg_chan_rate_nom[ih]->SetBinError(i,j, sqrt(susy_wg_chan_rate_nom[ih]->GetBinContent(i,j))*35.9*crosssection/noe);
					susy_wg_chan_rate_nom[ih]->SetBinContent(i,j, susy_wg_chan_rate_nom[ih]->GetBinContent(i,j)*35.9*crosssection/noe); 
					susy_wg_chan_rate_jesUp[ih]->SetBinContent(i,j, susy_wg_chan_rate_jesUp[ih]->GetBinContent(i,j)*35.9*crosssection/noe); 
					susy_wg_chan_rate_jesDown[ih]->SetBinContent(i,j, susy_wg_chan_rate_jesDown[ih]->GetBinContent(i,j)*35.9*crosssection/noe);
					susy_wg_chan_rate_jerUp[ih]->SetBinContent(i,j, susy_wg_chan_rate_jerUp[ih]->GetBinContent(i,j)*35.9*crosssection/noe);
					susy_wg_chan_rate_jerDown[ih]->SetBinContent(i,j, susy_wg_chan_rate_jerDown[ih]->GetBinContent(i,j)*35.9*crosssection/noe);
					susy_wg_chan_rate_xsUp[ih]->SetBinContent(i,j, susy_wg_chan_rate_xsUp[ih]->GetBinContent(i,j)*35.9*crosssectionUp/noe);
					susy_wg_chan_rate_esfUp[ih]->SetBinContent(i,j, susy_wg_chan_rate_esfUp[ih]->GetBinContent(i,j)*35.9*crosssection/noe);
									
					susy_wg_chan_syserr_jes[ih]->SetBinContent(i,j, max( fabs(susy_wg_chan_rate_jesUp[ih]->GetBinContent(i,j)-susy_wg_chan_rate_nom[ih]->GetBinContent(i,j)), fabs(susy_wg_chan_rate_jesDown[ih]->GetBinContent(i,j)-susy_wg_chan_rate_nom[ih]->GetBinContent(i,j)) ) ); 
					susy_wg_chan_syserr_jer[ih]->SetBinContent(i,j, max( fabs(susy_wg_chan_rate_jerUp[ih]->GetBinContent(i,j)-susy_wg_chan_rate_nom[ih]->GetBinContent(i,j)), fabs(susy_wg_chan_rate_jerDown[ih]->GetBinContent(i,j)-susy_wg_chan_rate_nom[ih]->GetBinContent(i,j)) ) ); 
					susy_wg_chan_syserr_esf[ih]->SetBinContent(i,j, fabs( susy_wg_chan_rate_esfUp[ih]->GetBinContent(i,j)-susy_wg_chan_rate_nom[ih]->GetBinContent(i,j)) );
					susy_wg_chan_syserr_scale[ih]->SetBinContent(i,j, -1);
					susy_wg_chan_syserr_eleshape[ih]->SetBinContent(i,j, -1);
					susy_wg_chan_syserr_jetshape[ih]->SetBinContent(i,j, -1);
					susy_wg_chan_syserr_qcdshape[ih]->SetBinContent(i,j, -1);
					susy_wg_chan_syserr_xs[ih]->SetBinContent(i,j, -1); 
					susy_wg_chan_syserr_lumi[ih]->SetBinContent(i,j, 0.026*susy_wg_chan_rate_nom[ih]->GetBinContent(i,j));   
					susy_wg_chan_syserr_isr[ih]->SetBinContent(i,j, -1);
					if(ih==0)susy_wg_syserr_PU->SetBinContent(i,j, analysis_PU(i,j, p_lowPU_susy_pass, p_lowPU_susy_all, p_highPU_susy_pass, p_highPU_susy_all, p_PU_data));  


					susy_gg_chan_rate_nom[ih]->SetBinError(i,j, sqrt(susy_gg_chan_rate_nom[ih]->GetBinContent(i,j))*35.9*crosssection/noe);
					susy_gg_chan_rate_nom[ih]->SetBinContent(i,j, susy_gg_chan_rate_nom[ih]->GetBinContent(i,j)*35.9*crosssection/noe); 
					susy_gg_chan_rate_jesUp[ih]->SetBinContent(i,j, susy_gg_chan_rate_jesUp[ih]->GetBinContent(i,j)*35.9*crosssection/noe); 
					susy_gg_chan_rate_jesDown[ih]->SetBinContent(i,j, susy_gg_chan_rate_jesDown[ih]->GetBinContent(i,j)*35.9*crosssection/noe);
					susy_gg_chan_rate_jerUp[ih]->SetBinContent(i,j, susy_gg_chan_rate_jerUp[ih]->GetBinContent(i,j)*35.9*crosssection/noe);
					susy_gg_chan_rate_jerDown[ih]->SetBinContent(i,j, susy_gg_chan_rate_jerDown[ih]->GetBinContent(i,j)*35.9*crosssection/noe);
					susy_gg_chan_rate_xsUp[ih]->SetBinContent(i,j, susy_gg_chan_rate_xsUp[ih]->GetBinContent(i,j)*35.9*crosssectionUp/noe);
					susy_gg_chan_rate_esfUp[ih]->SetBinContent(i,j, susy_gg_chan_rate_esfUp[ih]->GetBinContent(i,j)*35.9*crosssection/noe);
									
					susy_gg_chan_syserr_jes[ih]->SetBinContent(i,j, max( fabs(susy_gg_chan_rate_jesUp[ih]->GetBinContent(i,j)-susy_gg_chan_rate_nom[ih]->GetBinContent(i,j)), fabs(susy_gg_chan_rate_jesDown[ih]->GetBinContent(i,j)-susy_gg_chan_rate_nom[ih]->GetBinContent(i,j)) ) ); 
					susy_gg_chan_syserr_jer[ih]->SetBinContent(i,j, max( fabs(susy_gg_chan_rate_jerUp[ih]->GetBinContent(i,j)-susy_gg_chan_rate_nom[ih]->GetBinContent(i,j)), fabs(susy_gg_chan_rate_jerDown[ih]->GetBinContent(i,j)-susy_gg_chan_rate_nom[ih]->GetBinContent(i,j)) ) ); 
					susy_gg_chan_syserr_esf[ih]->SetBinContent(i,j, fabs( susy_gg_chan_rate_esfUp[ih]->GetBinContent(i,j)-susy_gg_chan_rate_nom[ih]->GetBinContent(i,j)) );
					susy_gg_chan_syserr_scale[ih]->SetBinContent(i,j, -1);
					susy_gg_chan_syserr_eleshape[ih]->SetBinContent(i,j, -1);
					susy_gg_chan_syserr_jetshape[ih]->SetBinContent(i,j, -1);
					susy_gg_chan_syserr_qcdshape[ih]->SetBinContent(i,j, -1);
					susy_gg_chan_syserr_xs[ih]->SetBinContent(i,j, -1); 
					susy_gg_chan_syserr_lumi[ih]->SetBinContent(i,j, 0.026*susy_gg_chan_rate_nom[ih]->GetBinContent(i,j));   
					susy_gg_chan_syserr_isr[ih]->SetBinContent(i,j, -1);
					if(ih==0)susy_gg_syserr_PU->SetBinContent(i,j, analysis_PU(i,j, p_lowPU_susy_pass, p_lowPU_susy_all, p_highPU_susy_pass, p_highPU_susy_all, p_PU_data));  
				}
			}
		} 
	} 

		for(int i(1); i < p_SUSYMass->GetXaxis()->GetNbins() + 1; i++){
			for(int j(1); j < p_SUSYMass->GetYaxis()->GetNbins() + 1; j++){
				float sparticleMass = p_SUSYMass->GetXaxis()->GetBinCenter(i);
				if( fabs(sparticleMass - 1700 ) < 10 && fabs( p_SUSYMass->GetYaxis()->GetBinCenter(j) - 1000) < 5){
					float noe = p_SUSYMass->GetBinContent(i,j);
					float crosssection = 0.5*p_crosssection_susy->GetBinContent( p_crosssection_susy->FindBin(sparticleMass) );
					p_susy_MET_signal_mg->Scale(35.8*crosssection/noe);
					p_susy_HT_signal_mg->Scale(35.8*crosssection/noe);
					p_susy_PhoEt_signal_mg->Scale(35.8*crosssection/noe);
					p_susy_MET_signal_eg->Scale(35.8*crosssection/noe);
					p_susy_HT_signal_eg->Scale(35.8*crosssection/noe);
					p_susy_PhoEt_signal_eg->Scale(35.8*crosssection/noe);
			}
		}
	}

	write_histo();
  xSecFile.Close();
}


