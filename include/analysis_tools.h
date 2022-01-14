enum MCType{
  NOMC = 0,
  WGJetInclusive = 1,
  WGJet40 = 2,
  WGJet130 = 3,
	ZGInclusive = 4,
	DYLL50 = 5,
	TTG = 6,
	WWG = 7,
	WZG = 8,
	DYLL10 = 9, // not used
	ZG130 = 10, // not used
	TT = 11,
	WW = 12,
	WZ = 13,
  W  = 14,
	QCDEM30 = 15, // not used
	QCDEM40 = 16, // not used 
	QCDMU = 17, // not used 
	GJet = 18, // not used 

	T5WgmG1800mLSP800 = 19,
	TChiWgmChi1000mLSP1 = 20,
	T5WgmG1800mLSP1600 = 21,
	T5WgmG1500mLSP1 = 22,
  generalMC = 23,
  nMCType = 24
};

// XSec for T5Wg_mG-1800_mLSP-800  = .001027
// XSec for TChiWg_mChi-1000_mLSP-1 = .001019
// update for WGJet40
double MC_XS[25] = {1, 412.7/*WGToLNuG*/, 19.75 /*WGJet40*/, 0.8099/*WGJet130*/, 55.48	/*ZGInclusive*/, 6404.0/*DY*/, 3.757/*TTG*/, 0.2147/*WWG*/, 0.04345/*WZG*/, 18610, 0.143, 750.5/*TTBar*/, 75.95/*WW*/, 27.59/*WZ*/, 61526/*W*/, 108000000/*QCDEM30*/, 54120000/*QCDEM40*/, 1/*MU*/, 16792/*GJet*/, 0.001027 /*T5Wg_mG-1800_mLSP-800*/, 0.001019 /*TChiWg_mChi-1000_mLSP-1*/, 0.001057/*T5Wg_mG-1800_mLSP-1600*/,   0.005763 /*T5Wg_mG-1500_mLSP-1*/, 1};

//0.143 // LO ZG130
//0.1404 // NLO ZG130

double getEvtWeight(int mcType, int year, double nEvents_MC){
    
    float luminosity;
    if(year==2016) luminosity = 35921.875595;
    if(year==2017) luminosity = 41529.548819;
    if(year==2018) luminosity = 59740.565202;
    
    double evtWeight = -1.;
    if(mcType==0) {evtWeight = 1.;}
    else evtWeight = MC_XS[mcType] * luminosity / nEvents_MC;
    
    //cout << "Using event weight " << evtWeight << endl;
    //cout << "XS = " << evtWeight/luminosity*nEvents_MC << endl;
    //cout << "lumi = " << luminosity << endl;
    //cout << "nEvents_MC = " << nEvents_MC << endl;
    
    return evtWeight;
}


float DeltaPhi(float phi1, float phi2){
  float deltaPhi = phi1 - phi2;
  if(fabs(deltaPhi) > TMath::Pi()){
    if(deltaPhi > 0)deltaPhi = -1.0*(TMath::TwoPi() - fabs(deltaPhi));
    else deltaPhi = TMath::TwoPi() - fabs(deltaPhi);
  }
  return deltaPhi;
}

float DeltaR(float eta1,float phi1,float eta2,float phi2)
{
	float deltaPhi = TMath::Abs(phi1-phi2);
	float deltaEta = eta1-eta2;
	if(deltaPhi > TMath::Pi())
	deltaPhi = TMath::TwoPi() - deltaPhi;
		return TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
}
// prompt photon if momID
// 1-6 for quarks, 11, 13, 15 leptons, 21, 22, 24 for gluons, photons, W+
// check for photon (PID 22) and electron(PID 11)
bool isHad(int PID, int momID){
   bool isFakePho;
   if(fabs(PID) == 22 || fabs(PID) == 11){
       switch(momID){
         case 1: isFakePho = false; break;
         case 2: isFakePho = false; break;
         case 3: isFakePho = false; break;
         case 4: isFakePho = false; break;
         case 5: isFakePho = false; break;
         case 6: isFakePho = false; break;
         case 999: isFakePho = false; break;
         case 11: isFakePho = false; break;
         case 13: isFakePho = false; break;
         case 15: isFakePho = false; break;
         case 22: isFakePho = false; break;
         case 24: isFakePho = false; break;
         case 21: isFakePho = false; break;
         default: isFakePho = true; break;
       }
  }
  else isFakePho = true;

  return isFakePho;
}

template<class recoParticle>
bool Veto(recoParticle reco, std::vector<recoPhoton>::iterator itpho, float dRcut){
  bool passVeto(true);
  typename recoParticle::iterator ire;
  for(ire = reco.begin(); ire != reco.end(); ire++)
    if(ire->getEt() > 2.0 && DeltaR(itpho->getEta(), itpho->getPhi(), ire->getEta(), ire->getPhi()) < dRcut)passVeto = false;
  return passVeto;
}


Double_t mybw(Double_t* x, Double_t* par)
{
  Double_t arg1 = 14.0/22.0; // 2 over pi
  Double_t arg2 = par[1]*par[1]*par[2]*par[2]; //Gamma=par[1]  M=par[2]
  Double_t arg3 = ((x[0]*x[0]) - (par[2]*par[2]))*((x[0]*x[0]) - (par[2]*par[2]));
  Double_t arg4 = x[0]*x[0]*x[0]*x[0]*((par[1]*par[1])/(par[2]*par[2]));
  Double_t arg5 = exp(par[3]+par[4]*x[0]);
  return par[0]*arg1*arg2/(arg3 + arg4)+ arg5;
}


Double_t mybwpol(Double_t* x, Double_t* par)
{
  Double_t arg1 = 14.0/22.0; // 2 over pi
  Double_t arg2 = par[1]*par[1]*par[2]*par[2]; //Gamma=par[1]  M=par[2]
  Double_t arg3 = ((x[0]*x[0]) - (par[2]*par[2]))*((x[0]*x[0]) - (par[2]*par[2]));
  Double_t arg4 = x[0]*x[0]*x[0]*x[0]*((par[1]*par[1])/(par[2]*par[2]));
  Double_t arg5 = par[3]*x[0]*x[0] + par[4]*x[0]+par[5]+par[6];
  return par[0]*arg1*arg2/(arg3 + arg4)+ arg5;
}

void DrawHisto(TObjArray *obj_list,TCanvas *canvas, int iPad, bool setLog){
  std::ostringstream name;
  int ihist(0);
  TObjArrayIter itr(obj_list);
  TH1F* hist(0);

  canvas->cd(iPad);
  if(setLog)gPad->SetLogy();
  TLegend *leg =  new TLegend(0.5,0.6,0.96,0.75);
  while( (hist = static_cast<TH1F*>(itr.Next())) ){
     ihist += 1;
     hist->SetLineColor(ihist+1);
     hist->SetLineWidth(2);
     hist->SetMarkerStyle(8);
     hist->SetMarkerColor(ihist+1);
     name.str("");
     name << hist->GetName();

     if(name.str().find("loose")!= std::string::npos){hist->SetLineColor(1);hist->SetMarkerColor(1);}
     else if(name.str().find("medium")!= std::string::npos){hist->SetLineColor(8);hist->SetMarkerColor(8);}
     else if(name.str().find("tight")!= std::string::npos){hist->SetLineColor(4);hist->SetMarkerColor(4);}
     leg->AddEntry(hist,name.str().c_str());
     if(ihist==1)hist->Draw();
     else hist->Draw("same");
  }
  leg->Draw("same");
}

int StartHEM = 319077;
double passHEMVeto(int mcType, rawData raw){
    if(mcType==0 && raw.run < StartHEM) return true;
    for(int iJet(0); iJet < raw.nJet; iJet++){
        if(recoJet(raw, iJet).getPt() > 30 && recoJet(raw, iJet).getEta() > -3.2 && recoJet(raw, iJet).getEta() < -1.2 && recoJet(raw, iJet).getPhi() > -1.77 && recoJet(raw, iJet).getPhi() < -0.67 && DeltaPhi(recoJet(raw, iJet).getPhi(),raw.pfMETPhi)<0.5) 
	return false;
    }

    for(int iEle(0); iEle < raw.nEle; iEle++){
	if(recoEle(raw, iEle).getCalibEt() > 30 && recoEle(raw, iEle).getEta() > -3.0 && recoEle(raw, iEle).getEta() < -1.4 && recoEle(raw, iEle).getPhi() > -1.57 && recoEle(raw, iEle).getPhi() < -0.87) 
	return false;
    }
        return true;   
}

unsigned findIndex(float* array, float kinevar,unsigned len){
  unsigned Index(0);
  for(unsigned i(0); i< len; i++){
    if(i != len-1){
      if(kinevar >= *(array+i) && kinevar < *(array+i+1))Index = i;}
    else{
      if(kinevar >= *(array+i))Index = i;
    } 
  }
  return Index;
}

float getPUESF(int nvertex){
// 2018
        float pileupweights[100];
        pileupweights[0] = 0;
        pileupweights[1] = 6.02567;
        pileupweights[2] = 1.39774;
        pileupweights[3] = 1.27393;
        pileupweights[4] = 0.875686;
        pileupweights[5] = 0.767948;
        pileupweights[6] = 1.00778;
        pileupweights[7] = 1.31716;
        pileupweights[8] = 1.33386;
        pileupweights[9] = 1.11231;
        pileupweights[10] = 0.914535;
        pileupweights[11] = 0.826661;
        pileupweights[12] = 0.803719;
        pileupweights[13] = 0.792745;
        pileupweights[14] = 0.808028;
        pileupweights[15] = 0.833057;
        pileupweights[16] = 0.857847;
        pileupweights[17] = 0.87476;
        pileupweights[18] = 0.884213;
        pileupweights[19] = 0.896318;
        pileupweights[20] = 0.921237;
        pileupweights[21] = 0.948391;
        pileupweights[22] = 0.967559;
        pileupweights[23] = 0.98067;
        pileupweights[24] = 0.99151;
        pileupweights[25] = 0.990645;
        pileupweights[26] = 0.985888;
        pileupweights[27] = 0.977224;
        pileupweights[28] = 0.977792;
        pileupweights[29] = 0.982772;
        pileupweights[30] = 0.98911;
        pileupweights[31] = 0.995026;
        pileupweights[32] = 0.998514;
        pileupweights[33] = 1.00086;
        pileupweights[34] = 1.00468;
        pileupweights[35] = 1.00648;
        pileupweights[36] = 1.00938;
        pileupweights[37] = 1.01238;
        pileupweights[38] = 1.01795;
        pileupweights[39] = 1.02256;
        pileupweights[40] = 1.03045;
        pileupweights[41] = 1.0391;
        pileupweights[42] = 1.04941;
        pileupweights[43] = 1.06078;
        pileupweights[44] = 1.07383;
        pileupweights[45] = 1.08321;
        pileupweights[46] = 1.09764;
        pileupweights[47] = 1.11119;
        pileupweights[48] = 1.12726;
        pileupweights[49] = 1.14584;
        pileupweights[50] = 1.15753;
        pileupweights[51] = 1.17012;
        pileupweights[52] = 1.1769;
        pileupweights[53] = 1.18909;
        pileupweights[54] = 1.19875;
        pileupweights[55] = 1.19625;
        pileupweights[56] = 1.20115;
        pileupweights[57] = 1.19815;
        pileupweights[58] = 1.20998;
        pileupweights[59] = 1.22022;
        pileupweights[60] = 1.23249;
        pileupweights[61] = 1.26372;
        pileupweights[62] = 1.2323;
        pileupweights[63] = 1.17456;
        pileupweights[64] = 1.09051;
        pileupweights[65] = 1.00584;
        pileupweights[66] = 0.898653;
        pileupweights[67] = 0.865695;
        pileupweights[68] = 0.787801;
        pileupweights[69] = 0.660631;
        pileupweights[70] = 0.636328;
        pileupweights[71] = 0.599542;
        pileupweights[72] = 0.641807;
        pileupweights[73] = 0.795867;
        pileupweights[74] = 0.646388;
        pileupweights[75] = 0.482347;
        pileupweights[76] = 0.430045;
        pileupweights[77] = 0.403128;
        pileupweights[78] = 0.447764;
        pileupweights[79] = 2.38134;
        pileupweights[80] = 0;
        pileupweights[81] = 0;
        pileupweights[82] = 0;
        pileupweights[83] = 0;
        pileupweights[84] = 0;
        pileupweights[85] = 0;
        pileupweights[86] = 0;
        pileupweights[87] = 0;
        pileupweights[88] = 0;
        pileupweights[89] = 0;
        pileupweights[90] = 0;
        pileupweights[91] = 0;
        pileupweights[92] = 0;
        pileupweights[93] = 0;
        pileupweights[94] = 0;
        pileupweights[95] = 0;
        pileupweights[96] = 0;
        pileupweights[97] = 0;
        pileupweights[98] = 0;
        pileupweights[99] = 0;
	if(nvertex > 99)return 0;
	else return pileupweights[nvertex]; 
}


//Double_t bkgEtBins[]={35,40,45,50,55,60,65,70,75,80, 85,90,95,100,105,110,115,120,125,130, 135,140,146,152,158,164,170,177,184,192, 200,208,216,224,232,240,250,260,275,290, 305,325,345,370,400,500,800};
Double_t bkgEtBins[]={35,40,50,60,70,80,90,100,110,120,130,140,150,160,170,185,200,215,230,250,275,290, 305,325,345,370,400,500,800};
int nBkgEtBins= sizeof(bkgEtBins)/sizeof(bkgEtBins[0]) -1;
//Double_t bkgPtBins[]={25,30,35,40,45,50,55,60,65,70,75,80, 85,90,95,100,105,110,115,120,125,130, 135,140,146,152,158,164,170,177,184,192, 200,208,216,224,232,240,250,260,275,290, 305,325,345,370,400,500,800};
Double_t bkgPtBins[]={25,50,75,100,125,150,200,400,800};
int nBkgPtBins= sizeof(bkgPtBins)/sizeof(bkgPtBins[0])-1;
Double_t bkgMETBins[]={0,40,60,80,100,120,140,160,180,210,240,280,320,400,600,1000};
int nBkgMETBins= sizeof(bkgMETBins)/sizeof(bkgMETBins[0]) -1;
Double_t bkgMtBins[]={0,20,40,60,80,100,120,140,160,180,200,300,400,500,1000};
int nBkgMtBins= sizeof(bkgMtBins)/sizeof(bkgMtBins[0]) -1;
Double_t bkgHTBins[]={0,40,60,80,100,120,140,160,180,200,225,250,275,300,340,380,420,500,600,1000};
int nBkgHTBins= sizeof(bkgHTBins)/sizeof(bkgHTBins[0]) -1;


//Double_t sigEtBins[]={35,50,100,150,200,300,500,800};
//int nSigEtBins= sizeof(sigEtBins)/sizeof(sigEtBins[0]) -1;
//Double_t sigPtBins[]={25,50,100,150,200,300,500,800};
//int nSigPtBins= sizeof(sigPtBins)/sizeof(sigPtBins[0])-1;
//Double_t sigMETBins[]={120,200,300,400,550,1000};
//int nSigMETBins= sizeof(sigMETBins)/sizeof(sigMETBins[0]) -1;
//Double_t sigMtBins[]={100,200,300,400,600,1000};
//int nSigMtBins= sizeof(sigMtBins)/sizeof(sigMtBins[0]) -1;
//Double_t sigHTBins[]={0,400,800,1500,2000};
//int nSigHTBins= sizeof(sigHTBins)/sizeof(sigHTBins[0]) -1;
//
Double_t sigEtBins[]={35,50,100,150,200,300,500,800};
int nSigEtBins= sizeof(sigEtBins)/sizeof(sigEtBins[0]) -1;
Double_t sigPtBins[]={25,50,100,150,200,300,500,800};
int nSigPtBins= sizeof(sigPtBins)/sizeof(sigPtBins[0])-1;
Double_t sigMETBins[]={0,20,40,60,80,100,120,150,200,250,300,350,400,600};
int nSigMETBins= sizeof(sigMETBins)/sizeof(sigMETBins[0]) -1;
Double_t sigMtBins[]={100,200,300,400,600,1000};
int nSigMtBins= sizeof(sigMtBins)/sizeof(sigMtBins[0]) -1;
Double_t sigHTBins[]={0,100,200,300,400,800,1500,2000};
int nSigHTBins= sizeof(sigHTBins)/sizeof(sigHTBins[0]) -1;
