#include "TLorentzVector.h"


  class recoJet{
  public:
  recoJet(rawData& raw, int ij){
    p4_.SetPtEtaPhiE((*raw.jetPt)[ij],(*raw.jetEta)[ij],(*raw.jetPhi)[ij],(*raw.jetEn)[ij]);
    jetArea_ = (*raw.jetArea)[ij];
		jetJECUnc_ = (*raw.jetJECUnc)[ij];
  }

  ~recoJet(){
  };
   
inline float getE() { return p4_.E();}
inline float getEt() { return p4_.Et();}
inline float getPt() { return p4_.Pt();}
inline float getEta() { return p4_.Eta();}
inline float getPhi() {return p4_.Phi();}
inline TLorentzVector getP4() {
         TLorentzVector p = p4_;
          return p;
}
inline float getArea() { return jetArea_;}
inline float getPtUnc(){ return jetJECUnc_;}

inline bool passSignalSelection(){
	bool pass(true);
	if(getPt() < 30.0)pass = false;
	if(fabs(getEta()) > 2.5)pass = false;
	return pass;
}

private:
    TLorentzVector p4_;
    float jetArea_;
		float jetJECUnc_;
};
   
