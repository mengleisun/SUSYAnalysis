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


double getEvtWeight(int year, double XSec, double nEvents_MC){
    
    float luminosity;
    if(year==2016) luminosity = 36470;
    //if(year==2017) luminosity = 41540;
    if(year==2017) luminosity = 27130;
    if(year==2018) luminosity = 59960;
    double evtWeight = 1.;
	
    evtWeight = XSec * luminosity / nEvents_MC;
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

float getPUESF16(int nvertex){
	
	float pileupweights[100];
	pileupweights[0] = 0;
        pileupweights[1] = 0;
        pileupweights[2] = 1.04404;
        pileupweights[3] = 0.671719;
        pileupweights[4] = 0.551154;
        pileupweights[5] = 0.571996;
        pileupweights[6] = 0.489252;
        pileupweights[7] = 0.600793;
        pileupweights[8] = 0.648089;
        pileupweights[9] = 0.70624;
        pileupweights[10] = 0.78362;
        pileupweights[11] = 0.775689;
        pileupweights[12] = 0.822737;
        pileupweights[13] = 0.900937;
        pileupweights[14] = 0.910491;
        pileupweights[15] = 0.918496;
        pileupweights[16] = 0.946654;
        pileupweights[17] = 0.964116;
        pileupweights[18] = 0.966068;
        pileupweights[19] = 1.00622;
        pileupweights[20] = 1.03616;
        pileupweights[21] = 1.03452;
        pileupweights[22] = 1.07318;
        pileupweights[23] = 1.05786;
        pileupweights[24] = 1.07856;
        pileupweights[25] = 1.10927;
        pileupweights[26] = 1.13568;
        pileupweights[27] = 1.15534;
        pileupweights[28] = 1.26312;
        pileupweights[29] = 1.22121;
        pileupweights[30] = 1.3091;
        pileupweights[31] = 1.38049;
        pileupweights[32] = 1.40534;
        pileupweights[33] = 1.62435;
        pileupweights[34] = 1.71987;
        pileupweights[35] = 1.69665;
        pileupweights[36] = 1.79234;
        pileupweights[37] = 2.21024;
        pileupweights[38] = 1.8312;
        pileupweights[39] = 2.66367;
        pileupweights[40] = 2.94597;
        pileupweights[41] = 2.70223;
        pileupweights[42] = 3.94469;
        pileupweights[43] = 3.60809;
        pileupweights[44] = 3.68486;
        pileupweights[45] = 3.55326;
        pileupweights[46] = 6.55086;
        pileupweights[47] = 3.37779;
        pileupweights[48] = 11.6687;
        pileupweights[49] = 0;
        pileupweights[50] = 3.5825;
        pileupweights[51] = 9.51922;
        pileupweights[52] = 1.99597;
        pileupweights[53] = 0;
        pileupweights[54] = 0;
        pileupweights[55] = 0;
        pileupweights[56] = 0;
        pileupweights[57] = 0;
        pileupweights[58] = 0;
        pileupweights[59] = 0;
        pileupweights[60] = 0;
        pileupweights[61] = 0;
        pileupweights[62] = 0.307072;
        pileupweights[63] = 0;
        pileupweights[64] = 0;
        pileupweights[65] = 0;
        pileupweights[66] = 0;
        pileupweights[67] = 0;
        pileupweights[68] = 0;
        pileupweights[69] = 0;
        pileupweights[70] = 0;
        pileupweights[71] = 0;
        pileupweights[72] = 0;
        pileupweights[73] = 0;
        pileupweights[74] = 0;
        pileupweights[75] = 0;
        pileupweights[76] = 0;
        pileupweights[77] = 0;
        pileupweights[78] = 0;
        pileupweights[79] = 0;
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
float getPUESF17(int nvertex){
	
	float pileupweights[100];	
	pileupweights[0] = 0;
        pileupweights[1] = 0;
        pileupweights[2] = 0.368589;
        pileupweights[3] = 5.52883;
        pileupweights[4] = 0.652119;
        pileupweights[5] = 0.700319;
        pileupweights[6] = 0.425295;
        pileupweights[7] = 0.455316;
        pileupweights[8] = 0.466681;
        pileupweights[9] = 0.436581;
        pileupweights[10] = 0.463855;
        pileupweights[11] = 0.483533;
        pileupweights[12] = 0.458275;
        pileupweights[13] = 0.476209;
        pileupweights[14] = 0.479196;
        pileupweights[15] = 0.528833;
        pileupweights[16] = 0.498161;
        pileupweights[17] = 0.592429;
        pileupweights[18] = 0.610549;
        pileupweights[19] = 0.616739;
        pileupweights[20] = 0.628302;
        pileupweights[21] = 0.659463;
        pileupweights[22] = 0.711438;
        pileupweights[23] = 0.738332;
        pileupweights[24] = 0.712789;
        pileupweights[25] = 0.770604;
        pileupweights[26] = 0.808979;
        pileupweights[27] = 0.822348;
        pileupweights[28] = 0.902373;
        pileupweights[29] = 0.903288;
        pileupweights[30] = 0.987908;
        pileupweights[31] = 1.03126;
        pileupweights[32] = 1.06761;
        pileupweights[33] = 1.10689;
        pileupweights[34] = 1.20031;
        pileupweights[35] = 1.2554;
        pileupweights[36] = 1.28501;
        pileupweights[37] = 1.39552;
        pileupweights[38] = 1.43056;
        pileupweights[39] = 1.47498;
        pileupweights[40] = 1.58264;
        pileupweights[41] = 1.67398;
        pileupweights[42] = 1.76129;
        pileupweights[43] = 1.65561;
        pileupweights[44] = 1.95346;
        pileupweights[45] = 2.0641;
        pileupweights[46] = 2.35761;
        pileupweights[47] = 2.31189;
        pileupweights[48] = 2.51311;
        pileupweights[49] = 2.64964;
        pileupweights[50] = 2.74274;
        pileupweights[51] = 3.17832;
        pileupweights[52] = 3.48279;
        pileupweights[53] = 4.42307;
        pileupweights[54] = 4.15277;
        pileupweights[55] = 4.2729;
        pileupweights[56] = 4.33371;
        pileupweights[57] = 5.15272;
        pileupweights[58] = 7.05585;
        pileupweights[59] = 5.57639;
        pileupweights[60] = 5.38363;
        pileupweights[61] = 7.39021;
        pileupweights[62] = 4.79166;
        pileupweights[63] = 9.70617;
        pileupweights[64] = 12.5781;
        pileupweights[65] = 15.9968;
        pileupweights[66] = 10.8734;
        pileupweights[67] = 9.32003;
        pileupweights[68] = 28.9342;
        pileupweights[69] = 45.705;
        pileupweights[70] = 9.58331;
        pileupweights[71] = 14.0064;
        pileupweights[72] = 16.2179;
        pileupweights[73] = 7.49464;
        pileupweights[74] = 14.1907;
        pileupweights[75] = 12.1634;
        pileupweights[76] = 8.29325;
        pileupweights[77] = 0;
        pileupweights[78] = 0;
        pileupweights[79] = 11.4263;
        pileupweights[80] = 0;
        pileupweights[81] = 0;
        pileupweights[82] = 0;
        pileupweights[83] = 5.16024;
        pileupweights[84] = 0;
        pileupweights[85] = 0;
        pileupweights[86] = 0;
        pileupweights[87] = 0;
        pileupweights[88] = 2.21153;
        pileupweights[89] = 0;
        pileupweights[90] = 0;
        pileupweights[91] = 0;
        pileupweights[92] = 0;
        pileupweights[93] = 0;
        pileupweights[94] = 2.58012;
        pileupweights[95] = 0;
        pileupweights[96] = 0;
        pileupweights[97] = 0;
        pileupweights[98] = 0;
        pileupweights[99] = 0;
	
	if(nvertex > 99)return 0;
	else return pileupweights[nvertex]; 
}

float getPUESF18(int nvertex){
	
	float pileupweights[100];
	pileupweights[0] = 0;
	pileupweights[1] = 0;
        pileupweights[2] = 1.28953;
        pileupweights[3] = 0.421578;
        pileupweights[4] = 0.82745;
        pileupweights[5] = 0.523873;
        pileupweights[6] = 0.800097;
        pileupweights[7] = 0.575908;
        pileupweights[8] = 0.681356;
        pileupweights[9] = 0.631048;
        pileupweights[10] = 0.735895;
        pileupweights[11] = 0.767841;
        pileupweights[12] = 0.733907;
        pileupweights[13] = 0.854597;
        pileupweights[14] = 0.71963;
        pileupweights[15] = 0.768989;
        pileupweights[16] = 0.7891;
        pileupweights[17] = 0.811136;
        pileupweights[18] = 0.795095;
        pileupweights[19] = 0.80085;
        pileupweights[20] = 0.831406;
        pileupweights[21] = 0.835328;
        pileupweights[22] = 0.875346;
        pileupweights[23] = 0.847136;
        pileupweights[24] = 0.88979;
        pileupweights[25] = 0.901946;
        pileupweights[26] = 0.902386;
        pileupweights[27] = 0.950084;
        pileupweights[28] = 0.9262;
        pileupweights[29] = 0.936081;
        pileupweights[30] = 0.968535;
        pileupweights[31] = 0.989897;
        pileupweights[32] = 1.01124;
        pileupweights[33] = 1.03438;
        pileupweights[34] = 0.998615;
        pileupweights[35] = 1.04456;
        pileupweights[36] = 1.17434;
        pileupweights[37] = 1.1385;
        pileupweights[38] = 1.20681;
        pileupweights[39] = 1.29955;
        pileupweights[40] = 1.36529;
        pileupweights[41] = 1.5202;
        pileupweights[42] = 1.49041;
        pileupweights[43] = 1.67056;
        pileupweights[44] = 1.82047;
        pileupweights[45] = 2.1055;
        pileupweights[46] = 2.1692;
        pileupweights[47] = 2.59291;
        pileupweights[48] = 2.65146;
        pileupweights[49] = 2.43148;
        pileupweights[50] = 3.33404;
        pileupweights[51] = 3.04206;
        pileupweights[52] = 3.85581;
        pileupweights[53] = 4.45841;
        pileupweights[54] = 4.90291;
        pileupweights[55] = 5.29153;
        pileupweights[56] = 5.35156;
        pileupweights[57] = 5.23522;
        pileupweights[58] = 6.20588;
        pileupweights[59] = 9.15568;
        pileupweights[60] = 6.47697;
        pileupweights[61] = 5.22013;
        pileupweights[62] = 8.3551;
        pileupweights[63] = 9.86493;
        pileupweights[64] = 19.343;
        pileupweights[65] = 11.6058;
        pileupweights[66] = 6.57662;
        pileupweights[67] = 9.83269;
        pileupweights[68] = 4.5456;
        pileupweights[69] = 0;
        pileupweights[70] = 0;
        pileupweights[71] = 6.77005;
        pileupweights[72] = 14.6684;
        pileupweights[73] = 12.4118;
        pileupweights[74] = 11.2834;
        pileupweights[75] = 0;
        pileupweights[76] = 0;
        pileupweights[77] = 0;
        pileupweights[78] = 6.12528;
        pileupweights[79] = 2.65966;
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
