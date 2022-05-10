#include "../include/analysis_rawData.h"

bool rawData::passHLT(){
  bool passDataHLT(false);
  
  switch(runtype_){
    case MC: passDataHLT = true; break;
    case DoubleEG2015: if(((HLTPho >> 14) &1) !=0)passDataHLT = true; break;
    case MuonEG2015: if(((HLTEleMuX >> 45) &1) !=0)passDataHLT = true; break; 
    case SingleElectron2015: if(((HLTEleMuX >> 6) &1) !=0)passDataHLT = true; break;
    case SingleMuon2015: if(((HLTEleMuX >> 25) &1)!=0)passDataHLT = true; break;
    case DoubleMuon2015: if(((HLTEleMuX >>20 ) &1) !=0 || ((HLTEleMuX >>21 ) &1)!=0)passDataHLT = true; break;

    // HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v
    case DoubleEG2016: if(((HLTPho >> 14) &1) !=0)passDataHLT = true; break;
    // HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v
    case DoubleEG2017: if(((HLTPho >> 14) &1) !=0)passDataHLT = true; break;
    case DoubleEG2018: if(((HLTPho >> 14) &1) !=0)passDataHLT = true; break;

    case MCDoubleEG2016: if(((HLTPho >> 14) &1) !=0)passDataHLT = true; break;
    case MCDoubleEG2017: if(((HLTPho >> 14) &1) !=0)passDataHLT = true; break;
    case MCDoubleEG2018: if(((HLTPho >> 14) &1) !=0)passDataHLT = true; break;
    
    // HLT_Mu17_Photon30_CaloIdL_L1ISO_v || HLT_Mu38NoFiltersNoVtx_Photon38_CaloIdL_v
    case MuonEG2016: if(((HLTEleMuX >> 8) &1) !=0 || ((HLTEleMuX >> 51) &1) !=0)passDataHLT = true; break;
    // HLT_Mu17_Photon30_IsoCaloId_v || HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_v
    case MuonEG2017: if(((HLTEleMuX >> 8) &1) !=0 || ((HLTEleMuX >> 57) &1) !=0)passDataHLT = true; break;
    case MuonEG2018: if(((HLTEleMuX >> 8) &1) !=0 || ((HLTEleMuX >> 57) &1) !=0)passDataHLT = true; break;

    case MCMuonEG2016: if(((HLTEleMuX >> 8) &1) !=0 || ((HLTEleMuX >> 51) &1) !=0)passDataHLT = true; break;
    case MCMuonEG2017: if(((HLTEleMuX >> 8) &1) !=0 || ((HLTEleMuX >> 57) &1) !=0)passDataHLT = true; break;
    case MCMuonEG2018: if(((HLTEleMuX >> 8) &1) !=0 || ((HLTEleMuX >> 57) &1) !=0)passDataHLT = true; break;

    case SingleElectron2016:if(((HLTEleMuX >> 4) &1) !=0)passDataHLT = true; break;  // HLT_Ele27_WPTight_Gsf_v
    case SingleElectron2017:if(((HLTEleMuX >> 3) &1) !=0)passDataHLT = true; break;  // HLT_Ele35_WPTight_Gsf_v
    case SingleElectron2018:if(((HLTEleMuX >> 55) &1) !=0)passDataHLT = true; break; // HLT_Ele32_WPTight_Gsf_v

    case SingleMuon2016: if(((HLTEleMuX >> 19) &1)!=0 || ((HLTEleMuX >> 20) &1)!=0)passDataHLT = true; break;// HLT_IsoMu24_v || HLT_IsoTkMu24_v
    case SingleMuon2017: if(((HLTEleMuX >> 56) &1)!=0)passDataHLT = true; break;	// HLT_IsoMu27_v
    case SingleMuon2018: if(((HLTEleMuX >> 19) &1)!=0)passDataHLT = true; break;	// HLT_IsoMu24_v

    case DoubleMuon2016: if(((HLTEleMuX >>14 ) &1) !=0 || ((HLTEleMuX >>15 ) &1)!=0 || ((HLTEleMuX >>16 ) &1)!=0 )passDataHLT = true; break;
    case DoubleMuon2017: if(((HLTEleMuX >>14 ) &1) !=0 || ((HLTEleMuX >>15 ) &1)!=0 )passDataHLT = true; break;
    case DoubleMuon2018: if(((HLTEleMuX >>15 ) &1) !=0 )passDataHLT = true; break;

    case MET2016: if(((HLTJet >>8 ) &1) !=0 || ((HLTJet >>24 ) &1)!=0 || ((HLTJet >>25 ) &1)!=0)passDataHLT = true; break;
    case MCSingleElectron2016:if(((HLTEleMuX >> 2) &1) !=0)passDataHLT = true; break;
    case MCSingleElectron2017:if(((HLTEleMuX >> 59) &1) !=0)passDataHLT = true; break;
    case MCSingleElectron2018:if(((HLTEleMuX >> 59) &1) !=0)passDataHLT = true; break;
    case MCSingleMuon2016: if(((HLTEleMuX >> 19) &1)!=0 || ((HLTEleMuX >> 20) &1)!=0)passDataHLT = true; break;
    case MCSingleMuon2017: if(((HLTEleMuX >> 56) &1)!=0)passDataHLT = true; break;	// HLT_IsoMu27_v
    case MCSingleMuon2018: if(((HLTEleMuX >> 19) &1)!=0)passDataHLT = true; break;	// HLT_IsoMu24_v
    case MCDoubleMuon2016: if(((HLTEleMuX >>14 ) &1) !=0 || ((HLTEleMuX >>15 ) &1)!=0 || ((HLTEleMuX >>16 ) &1)!=0 )passDataHLT = true; break;
    case MCMET2016: if(((HLTJet >>8 ) &1) !=0 || ((HLTJet >>24 ) &1)!=0 || ((HLTJet >>25 ) &1)!=0)passDataHLT = true; break;
    default: break;
  }
  
  return passDataHLT;
    
}
//string filterNamesToCheck[9] = {
//      "Flag_HBHENoiseFilter",
//      "Flag_HBHENoiseIsoFilter",
//      "Flag_globalSuperTightHalo2016Filter",
//      "Flag_goodVertices",
//      "Flag_eeBadScFilter",
//      "Flag_EcalDeadCellTriggerPrimitiveFilter",
//      "Flag_BadPFMuonFilter",
//      "Flag_BadChargedCandidateFilter", // not recommended
//      "Flag_BadPFMuonDzFilter"

bool rawData::passMETFilter(int filter){
  bool passfilter(true);
  for(int im(1); im <= 9; im++){
		if(im != 8){ // no ee Bad Sc filter
    	if(((filter >> im)&1)!=0){passfilter = false; return passfilter;}
		}
	}

//  switch(runtype_){
//    case DoubleEG2016: case MuonEG2016: case SingleElectron2016:case SingleMuon2016: case DoubleMuon2016: case MET2016: 
//									if(((filter >> 9)&1)!=1){passfilter = false; return passfilter;}
//									if(((filter >> 10)&1)!=1){passfilter = false; return passfilter;}
//									break; 
//    default: return passfilter; break;
//  }

  return passfilter;
}

int rawData::failFilterStep(int filter){
		int failStep(-1);
		
  	for(int im(1); im <= 8; im++){
			if(im != 5){ // no ee Bad Sc filter
    		if(((filter >> im)&1)!=0){failStep = im; return failStep;}
			}
		}
		
		return failStep;
}
