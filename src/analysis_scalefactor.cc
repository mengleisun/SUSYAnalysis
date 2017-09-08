#include "../include/analysis_scalefactor.h"

float esfScaleFactor::getElectronESF(float inputpt, float inputeta){
	
	bool  ptOnX(true);
	TAxis *ptaxis;
	TAxis *etaaxis;
	if(electronIDESF->GetXaxis()->GetXmax() > 3.0){
		ptaxis = electronIDESF->GetXaxis();
		etaaxis= electronIDESF->GetYaxis();
		ptOnX = true;
	}
	else{
		ptaxis = electronIDESF->GetYaxis();
		etaaxis= electronIDESF->GetXaxis();
		ptOnX = false;
	}	

	float pt = inputpt;
	if(pt > ptaxis->GetXmax())pt = ptaxis->GetXmax() - 1.0;
	float eta= inputeta;
	if(etaaxis->GetXmin() > -1.0)eta =fabs(eta);
	if(fabs(eta) > 2.5) eta = 2.49;

	float eleESF(1);
	if(ptOnX){
		eleESF = electronIDESF->GetBinContent(ptaxis->FindBin(pt), etaaxis->FindBin(eta))*electronISOESF->GetBinContent(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
	}
	else eleESF = electronIDESF->GetBinContent(etaaxis->FindBin(eta), ptaxis->FindBin(pt))*electronISOESF->GetBinContent(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
		
  return eleESF;
}


float esfScaleFactor::getPhotonESF(float inputpt, float inputeta){
	
	bool  ptOnX(true);
	TAxis *ptaxis;
	TAxis *etaaxis;
	if(photonIDESF->GetXaxis()->GetXmax() > 3.0){
		ptaxis = photonIDESF->GetXaxis();
		etaaxis= photonIDESF->GetYaxis();
		ptOnX = true;
	}
	else{
		ptaxis = photonIDESF->GetYaxis();
		etaaxis= photonIDESF->GetXaxis();
		ptOnX = false;
	}	

	float pt = inputpt;
	if(pt > ptaxis->GetXmax())pt = ptaxis->GetXmax() - 1.0;
	float eta= inputeta;
	if(etaaxis->GetXmin() > -1.0)eta =fabs(eta);
	if(fabs(eta) > 2.5) eta = 2.49;

	float phoESF(1);
	if(ptOnX){
		phoESF = photonIDESF->GetBinContent(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
	}
	else phoESF = photonIDESF->GetBinContent(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
		
  return phoESF;
}


float esfScaleFactor::getMuonESF(float inputpt, float inputeta){
	
	bool  ptOnX(true);
	TAxis *ptaxis;
	TAxis *etaaxis;
	if(muonIDESF->GetXaxis()->GetXmax() > 5.0){
		ptaxis = muonIDESF->GetXaxis();
		etaaxis= muonIDESF->GetYaxis();
		ptOnX = true;
	}
	else{
		ptaxis = muonIDESF->GetYaxis();
		etaaxis= muonIDESF->GetXaxis();
		ptOnX = false;
	}	

	float pt = inputpt;
	if(pt > ptaxis->GetXmax())pt = ptaxis->GetXmax() - 1.0;
	float eta= inputeta;
	if(etaaxis->GetXmin() > -1.0)eta =fabs(eta);
	if(fabs(eta) > 2.5) eta = 2.49;

	float muonESF(1);
	float muonIDESFvalue(1), muonISOESFvalue(1), muonIPESFvalue(1);
	if(ptOnX){
		muonIDESFvalue = muonIDESF->GetBinContent(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
		muonISOESFvalue= muonISOESF->GetBinContent(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
		muonIPESFvalue = muonIPESF->GetBinContent(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
		muonESF = muonIDESFvalue*muonISOESFvalue*muonIPESFvalue;
	}
	else{
	 muonIDESFvalue = muonIDESF->GetBinContent(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
	 muonISOESFvalue = muonISOESF->GetBinContent(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
	 muonIPESFvalue = muonIPESF->GetBinContent(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
	 muonESF = muonIDESFvalue*muonISOESFvalue*muonIPESFvalue;
	}
		
  return muonESF;
}

float esfScaleFactor::getegPhotonTRGESF(float inputpt, float inputeta){
	
	bool  ptOnX(true);
	TAxis *ptaxis;
	TAxis *etaaxis;
	if(photonTRGESF->GetXaxis()->GetXmax() > 3.0){
		ptaxis = photonTRGESF->GetXaxis();
		etaaxis= photonTRGESF->GetYaxis();
		ptOnX = true;
	}
	else{
		ptaxis = photonTRGESF->GetYaxis();
		etaaxis= photonTRGESF->GetXaxis();
		ptOnX = false;
	}	

	float pt = inputpt;
	if(pt > ptaxis->GetXmax())pt = ptaxis->GetXmax() - 1.0;
	float eta= inputeta;
	if(etaaxis->GetXmin() > -1.0)eta =fabs(eta);
	if(fabs(eta) > 2.5) eta = 2.49;

	float phoTRGESF(1);
	if(ptOnX){
		phoTRGESF = photonTRGESF->GetBinContent(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
	}
	else phoTRGESF = photonTRGESF->GetBinContent(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
		
  return phoTRGESF;
}

float esfScaleFactor::getElectronTRGESF(float inputpt, float inputeta){
	
	bool  ptOnX(true);
	TAxis *ptaxis;
	TAxis *etaaxis;
	if(electronTRGESF->GetXaxis()->GetXmax() > 3.0){
		ptaxis = electronTRGESF->GetXaxis();
		etaaxis= electronTRGESF->GetYaxis();
		ptOnX = true;
	}
	else{
		ptaxis = electronTRGESF->GetYaxis();
		etaaxis= electronTRGESF->GetXaxis();
		ptOnX = false;
	}	

	float pt = inputpt;
	if(pt > ptaxis->GetXmax())pt = ptaxis->GetXmax() - 1.0;
	float eta= inputeta;
	if(etaaxis->GetXmin() > -1.0)eta =fabs(eta);
	if(fabs(eta) > 2.5) eta = 2.49;

	float eleTRGESF(1);
	if(ptOnX){
		eleTRGESF = electronTRGESF->GetBinContent(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
	}
	else eleTRGESF = electronTRGESF->GetBinContent(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
		
  return eleTRGESF;
}

float esfScaleFactor::getMuonEGTRGESF(float pt1, float pt2){
	
	TAxis *phoptaxis = muonTRGESF->GetXaxis();
	TAxis *muptaxis  = muonTRGESF->GetYaxis();

	float phopt = pt1;
	if(phopt > phoptaxis->GetXmax())phopt = phoptaxis->GetXmax() - 1.0;
	float mupt = pt2;
	if(mupt > muptaxis->GetXmax())mupt = muptaxis->GetXmax() - 1.0;

	float muonegTRGESF = muonTRGESF->GetBinContent(phoptaxis->FindBin(phopt), muptaxis->FindBin(mupt));
		
  return muonegTRGESF; 
}


float esfScaleFactor::getElectronESFError(float inputpt, float inputeta){
	
	bool  ptOnX(true);
	TAxis *ptaxis;
	TAxis *etaaxis;
	if(electronIDESF->GetXaxis()->GetXmax() > 3.0){
		ptaxis = electronIDESF->GetXaxis();
		etaaxis= electronIDESF->GetYaxis();
		ptOnX = true;
	}
	else{
		ptaxis = electronIDESF->GetYaxis();
		etaaxis= electronIDESF->GetXaxis();
		ptOnX = false;
	}	

	float pt = inputpt;
	if(pt > ptaxis->GetXmax())pt = ptaxis->GetXmax() - 1.0;
	float eta= inputeta;
	if(etaaxis->GetXmin() > -1.0)eta =fabs(eta);
	if(fabs(eta) > 2.5) eta = 2.49;

	float eleESFError(1);
	if(ptOnX){
		eleESFError = electronIDESF->GetBinError(ptaxis->FindBin(pt), etaaxis->FindBin(eta))+electronISOESF->GetBinError(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
	}
	else eleESFError = electronIDESF->GetBinError(etaaxis->FindBin(eta), ptaxis->FindBin(pt))+electronISOESF->GetBinError(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
		
  return eleESFError;
}


float esfScaleFactor::getPhotonESFError(float inputpt, float inputeta){
	
	bool  ptOnX(true);
	TAxis *ptaxis;
	TAxis *etaaxis;
	if(photonIDESF->GetXaxis()->GetXmax() > 3.0){
		ptaxis = photonIDESF->GetXaxis();
		etaaxis= photonIDESF->GetYaxis();
		ptOnX = true;
	}
	else{
		ptaxis = photonIDESF->GetYaxis();
		etaaxis= photonIDESF->GetXaxis();
		ptOnX = false;
	}	

	float pt = inputpt;
	if(pt > ptaxis->GetXmax())pt = ptaxis->GetXmax() - 1.0;
	float eta= inputeta;
	if(etaaxis->GetXmin() > -1.0)eta =fabs(eta);
	if(fabs(eta) > 2.5) eta = 2.49;

	float phoESFError(1);
	if(ptOnX){
		phoESFError = photonIDESF->GetBinError(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
	}
	else phoESFError = photonIDESF->GetBinError(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
		
  return phoESFError;
}


float esfScaleFactor::getMuonESFError(float inputpt, float inputeta){
	
	bool  ptOnX(true);
	TAxis *ptaxis;
	TAxis *etaaxis;
	if(muonIDESF->GetXaxis()->GetXmax() > 5.0){
		ptaxis = muonIDESF->GetXaxis();
		etaaxis= muonIDESF->GetYaxis();
		ptOnX = true;
	}
	else{
		ptaxis = muonIDESF->GetYaxis();
		etaaxis= muonIDESF->GetXaxis();
		ptOnX = false;
	}	

	float pt = inputpt;
	if(pt > ptaxis->GetXmax())pt = ptaxis->GetXmax() - 1.0;
	float eta= inputeta;
	if(etaaxis->GetXmin() > -1.0)eta =fabs(eta);
	if(fabs(eta) > 2.5) eta = 2.49;

	float muonESFError(1);
	float muonIDESFError(1), muonISOESFError(1), muonIPESFError(1);
	if(ptOnX){
		muonIDESFError = muonIDESF->GetBinError(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
		muonISOESFError = muonISOESF->GetBinError(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
		muonIPESFError = muonIPESF->GetBinError(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
		muonESFError = muonIDESFError+muonISOESFError+muonIPESFError;
	}
	else{
	  muonIDESFError = muonIDESF->GetBinError(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
	  muonISOESFError = muonISOESF->GetBinError(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
	  muonIPESFError = muonIPESF->GetBinError(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
		muonESFError = muonIDESFError+muonISOESFError+muonIPESFError;
	}
		
  return muonESFError;
}

float esfScaleFactor::getegPhotonTRGESFError(float inputpt, float inputeta){
	
	bool  ptOnX(true);
	TAxis *ptaxis;
	TAxis *etaaxis;
	if(photonTRGESF->GetXaxis()->GetXmax() > 3.0){
		ptaxis = photonTRGESF->GetXaxis();
		etaaxis= photonTRGESF->GetYaxis();
		ptOnX = true;
	}
	else{
		ptaxis = photonTRGESF->GetYaxis();
		etaaxis= photonTRGESF->GetXaxis();
		ptOnX = false;
	}	

	float pt = inputpt;
	if(pt > ptaxis->GetXmax())pt = ptaxis->GetXmax() - 1.0;
	float eta= inputeta;
	if(etaaxis->GetXmin() > -1.0)eta =fabs(eta);
	if(fabs(eta) > 2.5) eta = 2.49;

	float phoTRGESFError(1);
	if(ptOnX){
		phoTRGESFError = photonTRGESF->GetBinError(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
	}
	else phoTRGESFError = photonTRGESF->GetBinError(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
		
  return phoTRGESFError;
}

float esfScaleFactor::getElectronTRGESFError(float inputpt, float inputeta){
	
	bool  ptOnX(true);
	TAxis *ptaxis;
	TAxis *etaaxis;
	if(electronTRGESF->GetXaxis()->GetXmax() > 3.0){
		ptaxis = electronTRGESF->GetXaxis();
		etaaxis= electronTRGESF->GetYaxis();
		ptOnX = true;
	}
	else{
		ptaxis = electronTRGESF->GetYaxis();
		etaaxis= electronTRGESF->GetXaxis();
		ptOnX = false;
	}	

	float pt = inputpt;
	if(pt > ptaxis->GetXmax())pt = ptaxis->GetXmax() - 1.0;
	float eta= inputeta;
	if(etaaxis->GetXmin() > -1.0)eta =fabs(eta);
	if(fabs(eta) > 2.5) eta = 2.49;

	float electronTRGESFError(1);
	if(ptOnX){
		electronTRGESFError = electronTRGESF->GetBinError(ptaxis->FindBin(pt), etaaxis->FindBin(eta));
	}
	else electronTRGESFError = electronTRGESF->GetBinError(etaaxis->FindBin(eta), ptaxis->FindBin(pt));
		
  return electronTRGESFError;
}

float esfScaleFactor::getMuonEGTRGESFError(float pt1, float pt2){
	
	TAxis *phoptaxis = muonTRGESF->GetXaxis();
	TAxis *muptaxis  = muonTRGESF->GetYaxis();

	float phopt = pt1;
	if(phopt > phoptaxis->GetXmax())phopt = phoptaxis->GetXmax() - 1.0;
	float mupt = pt2;
	if(mupt > muptaxis->GetXmax())mupt = muptaxis->GetXmax() - 1.0;

	float muonegTRGESFError = muonTRGESF->GetBinError(phoptaxis->FindBin(phopt), muptaxis->FindBin(mupt));
		
  return muonegTRGESFError; 
}
