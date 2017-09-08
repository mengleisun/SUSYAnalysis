#include "../include/analysis_photon.h"

bool recoPhoton::fireSingleTrg(int itrg) {
            if(((phoFiredSingleTrgs_ >> itrg)&1) == 1)return true;
            else{ return false;}
}

bool recoPhoton::fireDoubleTrg(int itrg) {
     if(((phoFiredDoubleTrgs_ >> itrg)&1) == 1)return true;
     else{ return false;}
}

bool recoPhoton::fireL1Trg(int itrg) {
     if(((phoFiredL1Trgs_ >> itrg)&1) == 1)return true;
     else{ return false;}
}

bool recoPhoton::isLoose() {
     if(((phoIDbit_ >> 0)&1)==1)return true;
     else return false;
}

bool recoPhoton::isMedium() {
     if(((phoIDbit_ >> 1)&1)==1)return true;
     else return false;
}

bool recoPhoton::isTight() {
     if(((phoIDbit_ >> 2)&1)==1)return true;
     else return false;
}

bool recoPhoton::isEB(){
     if(fabs(p4_.Eta())<1.444)return true;
     else return false;
}

bool recoPhoton::isEE(){
     if(fabs(p4_.Eta())>1.56 && fabs(p4_.Eta()) < 2.5)return true;
     else return false;
}

bool recoPhoton::passHoverE(int WP){
     bool pass(false);
     if(WP<=3)pass = phoHoverE_ < 0.05? true : false;
     return pass;
}

bool recoPhoton::passSigma(int WP){
     bool pass(false);
     if(fabs(p4_.Eta())<1.444){
       if(WP==1)pass = phoSigmaIEtaIEta_ < 0.01031? true:false;
       else if(WP ==2)pass = phoSigmaIEtaIEta_ < 0.01022?true:false;
       else if(WP ==3)pass = phoSigmaIEtaIEta_ < 0.00994? true:false;
     }
     else if(fabs(p4_.Eta())>1.56 && fabs(p4_.Eta()) < 2.5){
       if(WP==1)pass = phoSigmaIEtaIEta_ < 0.03013? true:false;
       else if(WP==2)pass = phoSigmaIEtaIEta_ < 0.03001? true:false;
       else if(WP ==3)pass = phoSigmaIEtaIEta_ <0.03000? true:false;
     } 
     return pass;
}


bool recoPhoton::passChIso(int WP){
     bool pass(false);
     if(fabs(p4_.Eta())<1.444){
       if(WP==1)pass = phoPFChIso_ < 1.295? true:false;
       else if(WP ==2)pass = phoPFChIso_ < 0.441? true:false;
       else if(WP ==3)pass = phoPFChIso_ < 0.202? true:false;
     }
     else if(fabs(p4_.Eta())>1.56 && fabs(p4_.Eta()) < 2.5){
       if(WP==1)pass = phoPFChIso_ < 1.011? true:false;
       else if(WP ==2)pass = phoPFChIso_ < 0.442? true:false;
       else if(WP ==3)pass = phoPFChIso_ < 0.034? true:false;
     } 
     return pass;
}

bool recoPhoton::passNeuIso(int WP){
     bool pass(false);
     if(fabs(p4_.Eta())<1.444){
       if(WP==1)pass = phoPFNeuIso_ < 10.910+0.0148*p4_.Et()+0.000017*p4_.Et()*p4_.Et()? true:false;
       else if(WP ==2)pass = phoPFNeuIso_ < 2.725+0.0148*p4_.Et()+0.000017*p4_.Et()*p4_.Et()? true:false;
       else if(WP ==3)pass = phoPFNeuIso_ < 0.264+0.0148*p4_.Et()+0.000017*p4_.Et()*p4_.Et()? true:false;
     }
     else if(fabs(p4_.Eta())>1.56 && fabs(p4_.Eta()) < 2.5){
       if(WP==1)pass = phoPFNeuIso_ < 5.931+0.0163*p4_.Et()+0.000014*p4_.Et()*p4_.Et()? true:false;
       else if(WP ==2)pass = phoPFNeuIso_ < 1.715+0.0163*p4_.Et()+0.000014*p4_.Et()*p4_.Et()? true:false;
       else if(WP ==3)pass = phoPFNeuIso_ < 0.586+0.0163*p4_.Et()+0.000014*p4_.Et()*p4_.Et()? true:false;
     } 
     return pass;
}

bool recoPhoton::passPhoIso(int WP){
     bool pass(false);
     if(fabs(p4_.Eta())<1.444){
       if(WP==1)pass = phoPFPhoIso_ < 3.630+0.0047*p4_.Et()? true:false;
       else if(WP ==2)pass = phoPFPhoIso_ < 2.571+0.0047*p4_.Et()? true:false;
       else if(WP ==3)pass = phoPFPhoIso_ < 2.362+0.0047*p4_.Et()? true:false;
     }
     else if(fabs(p4_.Eta())>1.56 && fabs(p4_.Eta()) < 2.5){
       if(WP==1)pass = phoPFPhoIso_ < 6.641+0.0034*p4_.Et()? true:false;
       else if(WP ==2)pass = phoPFPhoIso_ < 3.863+0.0034*p4_.Et()? true:false;
       else if(WP ==3)pass = phoPFPhoIso_ < 2.617+0.0034*p4_.Et()? true:false;
     } 
     return pass;
}



bool recoPhoton::passSignalSelection(){
  
  bool passCut(true);
  
  if(Calibp4_.Et() < 35.0){passCut = false; return passCut;} 
  if(!isLoose()){passCut = false; return passCut;}
  else passCut = true;

  return passCut; 
}

bool recoPhoton::passBasicSelection(){
  
  bool passCut = (passHoverE(1) && passNeuIso(1) && passPhoIso(1));
	if(getChIso() > 20)passCut = false;
  return passCut;
}

bool recoPhoton::passHLTSelection(){
  bool passHLT(false);

  switch(runtype_){
    case MC: passHLT = true; break;
    case DoubleEG2015: case MuonEG2015: case SingleElectron2015: case SingleMuon2015: case DoubleMuon2015:
      passHLT = false; break;
    case DoubleEG2016: if(fireDoubleTrg(5) || fireDoubleTrg(6))passHLT = true; break;
    case MuonEG2016: if(fireDoubleTrg(28)|| fireDoubleTrg(29) || fireDoubleTrg(30))passHLT = true; break;
    case SingleElectron2016: passHLT = true; break;
    case SingleMuon2016: passHLT = true; break;
    case DoubleMuon2016: passHLT = true; break;
    default: break;
  }
  
  return passHLT;
}

