#include "../include/analysis_muon.h"

bool recoMuon::passSignalSelection(){
  
  bool passCut(true);
  
  if(p4_.Pt() < 25.0){passCut = false; return passCut;} 
  if(fabs(p4_.Eta())>2.40){passCut = false; return passCut;}
  if(!isMedium()){passCut = false; return passCut;}
  if(getMiniIso()>0.2){passCut = false; return passCut;}
  //if(getRelIso()>0.15){passCut = false; return passCut;}
	if(fabs(getD0()) > 0.05 || fabs(getDz()) > 0.1){passCut = false; return passCut;}
  else passCut = true;

  return passCut; 
}

bool recoMuon::passHLTSelection(){
  bool passHLT(false);

  switch(runtype_){
    case MC: passHLT = true; break;
    case DoubleEG2015: passHLT = false; break;
    case MuonEG2015:   passHLT = false; break;
    case SingleElectron2015: passHLT = false; break;
    case SingleMuon2015: passHLT = false; break;
    case DoubleMuon2015: passHLT = false; break;
    case DoubleEG2016: passHLT = true; break;
    // HLT_Mu17_Photon || HLT_Mu17_Photon || HLTMu12_DIPho20
    case MuonEG2016: if(fireSingleTrg(2) || fireSingleTrg(21) || fireSingleTrg(22))passHLT = true; break;
    case MuonEG2017: if(fireSingleTrg(2) || fireSingleTrg(21) || fireSingleTrg(22))passHLT = true; break;
    case MuonEG2018: if(fireSingleTrg(2) || fireSingleTrg(21) || fireSingleTrg(22))passHLT = true; break;
    case SingleElectron2016: passHLT = true; break;
    case SingleElectron2017: passHLT = true; break;
    case SingleElectron2018: passHLT = true; break;
    // HLT_IsoMu24_v1 || HLT_IsoTkMu24
    case SingleMuon2016: if(fireSingleTrg(1) || fireSingleTrg(19))passHLT = true; break;
    case SingleMuon2017: if(fireSingleTrg(31))passHLT = true; break;
    case SingleMuon2018: if(fireSingleTrg(30))passHLT = true; break;
    case DoubleMuon2016: passHLT = true; break;
    case MCDoubleEG2016: passHLT = true; break;
    case MCMuonEG2016: if(fireSingleTrg(2) || fireSingleTrg(21) || fireSingleTrg(22))passHLT = true; break;
    case MCSingleElectron2016:passHLT = true; break;
    case MCSingleMuon2016: if(fireSingleTrg(1) || fireSingleTrg(19))passHLT = true; break;
    case MCDoubleMuon2016: passHLT = true; break;
    case MCMET2016: passHLT = true; break;
    default: break;
  }
  
  return passHLT;
}

bool recoMuon::isFakeProxy(){
	bool passCut(true);

  if(p4_.Pt() < 25.0){passCut = false; return passCut;}
  if(fabs(p4_.Eta())>2.40){passCut = false; return passCut;}
  //if(!isMedium()){passCut = false; return passCut;}
  if(!isLoose()){passCut = false; return passCut;}
  if(getD0() > 0.05 || getDz() > 0.1){passCut = false; return passCut;}
  if(isMedium() && getMiniIso() <= 0.2){passCut = false; return passCut;}
  else passCut = true;

  return passCut;
}
 
