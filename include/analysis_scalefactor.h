#ifndef ANALYSIS_SCALEFACTOR 
#define ANALYSIS_SCALEFACTOR

#ifndef ROOT_TFile
#include "TFile.h"
#endif
#ifndef ROOT_TH2F
#include "TH2F.h"
#endif
#ifndef ROOT_TAxis
#include "TAxis.h"
#endif
#include<iostream>

  class esfScaleFactor{
  public:
    esfScaleFactor()
  {
		TFile *electronIDFile = TFile::Open("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/sf/ele_scaleFactors.root");
		TFile *photonIDFile   = TFile::Open("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/sf/egammaEffi.txt_EGM2D.root");
		TFile *muonIDFile     = TFile::Open("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/sf/TnP_NUM_MediumID_DENOM_generalTracks_VAR_map_pt_eta.root");
		TFile *muonIsoFile    = TFile::Open("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/sf/TnP_NUM_MiniIsoTight_DENOM_MediumID_VAR_map_pt_eta.root");
		//TFile *muonIsoFile    = TFile::Open("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/sf/ratio_NUM_RelIsoVTight_DENOM_MediumID_VAR_map_pt_eta.root");
		TFile *muonIPFile     = TFile::Open("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/sf/TnP_NUM_TightIP2D_DENOM_MediumID_VAR_map_pt_eta.root");
		TFile *photonTRGFile  = TFile::Open("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/sf/diphoton_pholeg.root");
		TFile *electronTRGFile  = TFile::Open("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/sf/diphoton_eleg.root");
		TFile *muonTRGFile    = TFile::Open("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/sf/muonphoton_trigger.root");

		if(electronIDFile->IsZombie()){
			std::cout << "no electron ESF file" << std::endl;
		}
		else{
			electronIDESF = (TH2F*)electronIDFile->Get("GsfElectronToCutBasedSpring15M");
			electronISOESF= (TH2F*)electronIDFile->Get("MVAVLooseElectronToMini"); 
		}

		if(photonIDFile->IsZombie()){
			std::cout << "no photon ESF file" << std::endl;
		}
		else{
			photonIDESF   = (TH2F*)photonIDFile->Get("EGamma_SF2D");
		}

		if(muonIDFile->IsZombie()){
      std::cout << "no muon ESF file" << std::endl;
    }
		else{
			muonIDESF    = (TH2F*)muonIDFile->Get("SF");
		}

		if(muonIsoFile->IsZombie()){
      std::cout << "no muon Iso ESF file" << std::endl;
    }
		else{
			muonISOESF    = (TH2F*)muonIsoFile->Get("SF");
			//muonISOESF    = (TH2F*)muonIsoFile->Get("cr");
		}

		if(muonIPFile->IsZombie()){
      std::cout << "no muon ESF file" << std::endl;
    }
		else{
			muonIPESF    = (TH2F*)muonIPFile->Get("SF");
		}


		if(photonTRGFile->IsZombie()){
			std::cout << "no photon TRG ESF file" << std::endl;
		}
		else{
			photonTRGESF   = (TH2F*)photonTRGFile->Get("ESF_Lead");
		}

		if(electronTRGFile->IsZombie()){
			std::cout << "no electron TRG ESF file" << std::endl;
		}
		else{
			electronTRGESF   = (TH2F*)electronTRGFile->Get("p_TriggerScale");
		}

		if(muonTRGFile->IsZombie()){
			std::cout << "no muon TRG ESF file" << std::endl;
		}
		else{
			muonTRGESF   = (TH2F*)muonTRGFile->Get("p_mgESF");
		}
	}

	~esfScaleFactor(){
	};

	float getElectronESF(float pt, float eta);
	float getPhotonESF(float pt, float eta);
	float getMuonESF(float pt, float eta);
	float getElectronESFError(float pt, float eta);
	float getPhotonESFError(float pt, float eta);
	float getMuonESFError(float pt, float eta);
	float getElectronTRGESF(float pt, float eta);
	float getegPhotonTRGESF(float pt, float eta);
	float getMuonEGTRGESF(float pt1, float pt2);
	float getegPhotonTRGESFError(float pt, float eta);
	float getElectronTRGESFError(float pt, float eta);
	float getMuonEGTRGESFError(float pt1, float pt2);
	
	private:
		TH2F *electronIDESF;
		TH2F *electronISOESF;
		TH2F *photonIDESF;
		TH2F *muonIDESF;
		TH2F *muonISOESF;
		TH2F *muonIPESF;
		TH2F *electronTRGESF;
		TH2F *photonTRGESF;
		TH2F *muonTRGESF;
	};

#endif

