#include "../include/analysis_ele.h"

bool recoEle::passSignalSelection(){
  
  bool passCut(true);
  
  if(Calibp4_.Pt() < 25.0){passCut = false; return passCut;} 
  if(!isEB() && !isEE()){passCut = false; return passCut;}
  if(!isMiniMedium()){passCut = false; return passCut;}
  else passCut = true;

  return passCut; 
}

bool recoEle::passHLTSelection(){
  bool passHLT(false);

  switch(runtype_){
    case MC: passHLT = true; break;
    case DoubleEG2015: case MuonEG2015: case SingleMuon2015: case DoubleMuon2015: case SingleElectron2015:
      passHLT = false; break;
    case DoubleEG2016:  if(fireTrgs(21) || fireTrgs(22))passHLT = true; break;
    case MuonEG2016: passHLT = true; break;
    case SingleElectron2016: if(fireTrgs(10))passHLT = true; break;
    case SingleMuon2016: passHLT = true; break;
    case DoubleMuon2016: passHLT = true; break;
    case MCDoubleEG: if(fireTrgs(21) || fireTrgs(22))passHLT = true; break; 
    case MCMuonEG:  passHLT = true; break; 
    case MCSingleElectron: if(fireTrgs(10))passHLT = true; break;
    case MCSingleMuon: passHLT = true; break;
    case MCDoubleMuon: passHLT = true; break;
    case MCMET: passHLT = true; break;
    default: break;
  }
  
  return passHLT;
}

bool recoEle::passBasicID(){
  bool passBasic(true);
 
  if(!isEB() && !isEE()) {passBasic=false; return passBasic; }
  if(isEB()){
		if(getHoverE() > 0.253){passBasic=false; return passBasic;}
		if(fabs(getEoverPInv()) > 0.134){passBasic=false; return passBasic;}
		if(getMissHits() > 1){passBasic=false; return passBasic;}
		if(getConvVeto() == 0){passBasic=false; return passBasic;}
  }
  else if(isEE()){
		if(getHoverE() > 0.0878){passBasic=false; return passBasic;}
		if(fabs(getEoverPInv()) > 0.13){passBasic=false; return passBasic;}
		if(getMissHits() > 1){passBasic=false; return passBasic;}
		if(getConvVeto() == 0){passBasic=false; return passBasic;}
  }
  return passBasic;
}

bool recoEle::isMiniMedium(){
  bool passMiniMedium(true);
 
  if(!isEB() && !isEE()) {passMiniMedium=false; return passMiniMedium; }
  if(isEB()){
	if(getSigma() > 0.00998) {passMiniMedium=false; return passMiniMedium; }
	if(fabs(getdEtaIn()) > 0.00311){passMiniMedium=false; return passMiniMedium;}
	if(fabs(getdPhiIn()) > 0.103){passMiniMedium=false; return passMiniMedium;}
	if(getHoverE() > 0.253){passMiniMedium=false; return passMiniMedium;}
	if(fabs(getEoverPInv()) > 0.134){passMiniMedium=false; return passMiniMedium;}
//	if(fabs(getD0()) > 0.05){passMiniMedium=false; return passMiniMedium;}
//	if(fabs(getDz()) > 0.10){passMiniMedium=false; return passMiniMedium; }
	if(getMissHits() > 1){passMiniMedium=false; return passMiniMedium;}
	if(getConvVeto() == 0){passMiniMedium=false; return passMiniMedium;}
	if(getMiniIso() > 0.1){passMiniMedium=false; return passMiniMedium;}
  }
  else if(isEE()){
	if(getSigma() > 0.0298) {passMiniMedium=false; return passMiniMedium; }
	if(fabs(getdEtaIn()) > 0.00609){passMiniMedium=false; return passMiniMedium;}
	if(fabs(getdPhiIn()) > 0.045){passMiniMedium=false; return passMiniMedium;}
	if(getHoverE() > 0.0878){passMiniMedium=false; return passMiniMedium;}
	if(fabs(getEoverPInv()) > 0.13){passMiniMedium=false; return passMiniMedium;}
//	if(fabs(getD0()) > 0.10){passMiniMedium=false; return passMiniMedium;}
//	if(fabs(getDz()) > 0.20){passMiniMedium=false; return passMiniMedium; }
	if(getMissHits() > 1){passMiniMedium=false; return passMiniMedium;}
	if(getConvVeto() == 0){passMiniMedium=false; return passMiniMedium;}
	if(getMiniIso() > 0.1){passMiniMedium=false; return passMiniMedium;}
  }
  return passMiniMedium;
}


bool recoEle::isMiniLoose(){
  bool passMiniLoose(true);
 
  if(!isEB() && !isEE()) {passMiniLoose=false; return passMiniLoose; }
  if(isEB()){
	if(getSigma() > 0.011) {passMiniLoose=false; return passMiniLoose; }
	if(fabs(getdEtaIn()) > 0.00477){passMiniLoose=false; return passMiniLoose;}
	if(fabs(getdPhiIn()) > 0.222){passMiniLoose=false; return passMiniLoose;}
	if(getHoverE() > 0.298){passMiniLoose=false; return passMiniLoose;}
	if(fabs(getEoverPInv()) > 0.241){passMiniLoose=false; return passMiniLoose;}
//	if(fabs(getD0()) > 0.05){passMiniLoose=false; return passMiniLoose;}
//	if(fabs(getDz()) > 0.10){passMiniLoose=false; return passMiniLoose; }
	if(getMissHits() > 1){passMiniLoose=false; return passMiniLoose;}
	if(!getConvVeto()){passMiniLoose=false; return passMiniLoose;}
	if(getMiniIso() > 0.2){passMiniLoose=false; return passMiniLoose;}
  }
  else if(isEE()){
	if(getSigma() > 0.0314) {passMiniLoose=false; return passMiniLoose; }
	if(fabs(getdEtaIn()) > 0.00868){passMiniLoose=false; return passMiniLoose;}
	if(fabs(getdPhiIn()) > 0.213){passMiniLoose=false; return passMiniLoose;}
	if(getHoverE() > 0.101){passMiniLoose=false; return passMiniLoose;}
	if(fabs(getEoverPInv()) > 0.14){passMiniLoose=false; return passMiniLoose;}
//	if(fabs(getD0()) > 0.10){passMiniLoose=false; return passMiniLoose;}
//	if(fabs(getDz()) > 0.20){passMiniLoose=false; return passMiniLoose; }
	if(getMissHits() > 1){passMiniLoose=false; return passMiniLoose;}
	if(!getConvVeto()){passMiniLoose=false; return passMiniLoose;}
	if(getMiniIso() > 0.2){passMiniLoose=false; return passMiniLoose;}
  }
  return passMiniLoose;
}

bool recoEle::isFakeProxy(){
	bool passFakeProxy(true);

  if(!isEB() && !isEE()) {passFakeProxy=false; return passFakeProxy; }
  if(isEB()){
	if(getHoverE() > 0.253){passFakeProxy=false; return passFakeProxy;}
	if(fabs(getEoverPInv()) > 0.134){passFakeProxy=false; return passFakeProxy;}
	if(getMissHits() > 1){passFakeProxy=false; return passFakeProxy;}
	if(!getConvVeto()){passFakeProxy=false; return passFakeProxy;}
  }
  else if(isEE()){
	if(getHoverE() > 0.0878){passFakeProxy=false; return passFakeProxy;}
	if(fabs(getEoverPInv()) > 0.13){passFakeProxy=false; return passFakeProxy;}
	if(getMissHits() > 1){passFakeProxy=false; return passFakeProxy;}
	if(!getConvVeto()){passFakeProxy=false; return passFakeProxy;}
  }
	
  if(isEB()){
  	if(getSigma() > 0.00998 || fabs(getdEtaIn()) > 0.00311 || fabs(getdPhiIn()) > 0.103 || getMiniIso() > 0.1)passFakeProxy=true; 
  	else passFakeProxy=false;
  }
  else if(isEE()){
  	if(getSigma() > 0.0298 || fabs(getdEtaIn()) > 0.00609 || fabs(getdPhiIn()) > 0.045 || getMiniIso() > 0.1)passFakeProxy=true;
  	else passFakeProxy=false;
  }

//Upper bound:  Veto
  if(isEB()){
  	if(getSigma() > 0.0115 || fabs(getdEtaIn()) > 0.00749 || fabs(getdPhiIn()) > 0.228)passFakeProxy=false; 
  }
  else if(isEE()){
  	if(getSigma() > 0.037  || fabs(getdEtaIn()) > 0.00895 || fabs(getdPhiIn()) > 0.213)passFakeProxy=false;
  }
	return passFakeProxy;
}


bool recoEle::isLooseFakeProxy(){
	bool passFakeProxy(true);

  if(!isEB() && !isEE()) {passFakeProxy=false; return passFakeProxy; }
  if(isEB()){
	if(getHoverE() > 0.253){passFakeProxy=false; return passFakeProxy;}
	if(fabs(getEoverPInv()) > 0.134){passFakeProxy=false; return passFakeProxy;}
	if(getMissHits() > 1){passFakeProxy=false; return passFakeProxy;}
	if(!getConvVeto()){passFakeProxy=false; return passFakeProxy;}
  }
  else if(isEE()){
	if(getHoverE() > 0.0878){passFakeProxy=false; return passFakeProxy;}
	if(fabs(getEoverPInv()) > 0.13){passFakeProxy=false; return passFakeProxy;}
	if(getMissHits() > 1){passFakeProxy=false; return passFakeProxy;}
	if(!getConvVeto()){passFakeProxy=false; return passFakeProxy;}
  }
	
  if(isEB()){
  	if(getSigma() > 0.00998 || fabs(getdEtaIn()) > 0.00311 || fabs(getdPhiIn()) > 0.103 || getMiniIso() > 0.1)passFakeProxy=true; 
  	else passFakeProxy=false;
  }
  else if(isEE()){
  	if(getSigma() > 0.0298 || fabs(getdEtaIn()) > 0.00609 || fabs(getdPhiIn()) > 0.045 || getMiniIso() > 0.1)passFakeProxy=true;
  	else passFakeProxy=false;
  }

	return passFakeProxy;

}
