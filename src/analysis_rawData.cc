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
    case DoubleEG2016: if(((HLTPho >> 14) &1) !=0)passDataHLT = true; break; 
    case MuonEG2016: if(((HLTEleMuX >> 8) &1) !=0 || ((HLTEleMuX >> 41) &1) !=0)passDataHLT = true; break;
    case SingleElectron2016:if(((HLTEleMuX >> 2) &1) !=0)passDataHLT = true; break;
    case SingleMuon2016: if(((HLTEleMuX >> 19) &1)!=0 || ((HLTEleMuX >> 20) &1)!=0)passDataHLT = true; break; 
    case DoubleMuon2016: if(((HLTEleMuX >>14 ) &1) !=0 || ((HLTEleMuX >>15 ) &1)!=0 || ((HLTEleMuX >>16 ) &1)!=0 )passDataHLT = true; break;
    case MET2016: if(((HLTJet >>8 ) &1) !=0 || ((HLTJet >>24 ) &1)!=0 || ((HLTJet >>25 ) &1)!=0)passDataHLT = true; break;
    case MCDoubleEG: if(((HLTPho >> 14) &1) !=0)passDataHLT = true; break; 
    case MCMuonEG: if(((HLTEleMuX >> 8) &1) !=0 || ((HLTEleMuX >> 41) &1) !=0)passDataHLT = true; break;
    case MCSingleElectron:if(((HLTEleMuX >> 2) &1) !=0)passDataHLT = true; break;
    case MCSingleMuon: if(((HLTEleMuX >> 19) &1)!=0 || ((HLTEleMuX >> 20) &1)!=0)passDataHLT = true; break; 
    case MCDoubleMuon: if(((HLTEleMuX >>14 ) &1) !=0 || ((HLTEleMuX >>15 ) &1)!=0 || ((HLTEleMuX >>16 ) &1)!=0 )passDataHLT = true; break;
    case MCMET: if(((HLTJet >>8 ) &1) !=0 || ((HLTJet >>24 ) &1)!=0 || ((HLTJet >>25 ) &1)!=0)passDataHLT = true; break;
    default: break;
  }
  
  return passDataHLT;
    
}

bool rawData::passMETFilter(int filter){
  bool passfilter(true);
  for(int im(1); im <= 8; im++)
    if(((filter >> im)&1)!=0){passfilter = false; return passfilter;}

  switch(runtype_){
    case DoubleEG2016: case MuonEG2016: case SingleElectron2016:case SingleMuon2016: case DoubleMuon2016: case MET2016: 
									if(((filter >> 9)&1)!=1){passfilter = false; return passfilter;}
									if(((filter >> 10)&1)!=1){passfilter = false; return passfilter;}
									break; 
    default: return passfilter; break;
  }

  return passfilter;
}

