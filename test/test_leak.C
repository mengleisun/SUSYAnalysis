#include "analysis_commoncode.h"

void test_leak(){//main  

	SetSignalConfig();
	binning Bin(NBIN, METbin1, METbin2, HTbin1, HTbin2, PHOETbin);

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

  TChain* es = new TChain("SUSYtree");
  //es->Add("/uscms_data/d3/mengleis/Sep1/resTree_TChiWG.root");
  //es->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_T5WG.root");
  es->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_GMSB.root");

	float Mgluino(0);
  float Mchargino(0);
  float Mneutralino(0);
  float mcPhotonEt(0);
  float mcPhotonEta(0);
  float mcPhotonPhi(0);
  float recoPhotonEt(0);
  float recoPhotonEta(0);
  float recoPhotonPhi(0);
  float PhoR9(0);
  float PhoHoverE(0);
  float PhoSigma(0);
  float PhoChIso(0);
  float PhoNeuIso(0);
  float PhoPhoIso(0);
  bool  PhoPassID;
  float mcElePt(0);
  float mcEleEta(0);
  float mcElePhi(0);
  float recoEleEt(0);
  float recoEleEta(0);
  float recoElePhi(0);
  float EleR9(0);
  float EleHoverE(0);
  float EleSigma(0);
  float EleChIso(0);
  float EleNeuIso(0);
  float ElePhoIso(0);
  float EleMiniIso(0);
  float EledEtaIn(0);
  float EledPhiIn(0);
  float EleD0(0);
  float EleDz(0);
  float EleooEmooP(0);
  bool  ElePassID;
  float dRPhoEle(0);
  float Invmass(0);
  int   nJet(0);
  float HT(0);
	float sigMET(0);
  double MT_(0), ThreeBodyMass_(0);
	int   nVertex(0);

  //es->SetBranchAddress("MsGsQ",        &Mgluino);
  es->SetBranchAddress("Mchargino",      &Mchargino);
  es->SetBranchAddress("Mneutralino",    &Mneutralino);
	es->SetBranchAddress("nVertex",        &nVertex);
  es->SetBranchAddress("MT",&MT_);
  es->SetBranchAddress("ThreeBodyMass",&ThreeBodyMass_);
  es->SetBranchAddress("mcPhotonEt"     ,&mcPhotonEt); 
  es->SetBranchAddress("mcPhotonEta"    ,&mcPhotonEta);
  es->SetBranchAddress("mcPhotonPhi"    ,&mcPhotonPhi);
  es->SetBranchAddress("recoPhotonEt",   &recoPhotonEt);
  es->SetBranchAddress("recoPhotonEta",  &recoPhotonEta);
  es->SetBranchAddress("recoPhotonPhi",  &recoPhotonPhi);
  es->SetBranchAddress("PhoR9"          ,&PhoR9);
  es->SetBranchAddress("PhoHoverE"      ,&PhoHoverE);
  es->SetBranchAddress("PhoSigma"       ,&PhoSigma);
  es->SetBranchAddress("PhoChIso"       ,&PhoChIso);
  es->SetBranchAddress("PhoNeuIso"      ,&PhoNeuIso);
  es->SetBranchAddress("PhoPhoIso"      ,&PhoPhoIso);
  es->SetBranchAddress("PhoPassID"      ,&PhoPassID);
  es->SetBranchAddress("mcElePt"        ,&mcElePt);
  es->SetBranchAddress("mcEleEta"       ,&mcEleEta);
  es->SetBranchAddress("mcElePhi"       ,&mcElePhi);
  es->SetBranchAddress("recoEleEt",     &recoEleEt);
  es->SetBranchAddress("recoEleEta",    &recoEleEta);
  es->SetBranchAddress("recoElePhi",    &recoElePhi);
  es->SetBranchAddress("EleR9"          ,&EleR9);
  es->SetBranchAddress("EleHoverE"      ,&EleHoverE);
  es->SetBranchAddress("EleSigma"       ,&EleSigma);
  es->SetBranchAddress("EleChIso"       ,&EleChIso);
  es->SetBranchAddress("EleNeuIso"      ,&EleNeuIso);
  es->SetBranchAddress("ElePhoIso"      ,&ElePhoIso);
  es->SetBranchAddress("EleMiniIso"     ,&EleMiniIso);
  es->SetBranchAddress("EledEtaIn"      ,&EledEtaIn);
  es->SetBranchAddress("EledPhiIn"      ,&EledPhiIn);
  es->SetBranchAddress("EleD0"          ,&EleD0);
  es->SetBranchAddress("EleDz"          ,&EleDz);
  es->SetBranchAddress("EleooEmooP"     ,&EleooEmooP);
  es->SetBranchAddress("ElePassID"      ,&ElePassID);
  es->SetBranchAddress("dRPhoEle"       ,&dRPhoEle);
  es->SetBranchAddress("Invmass"        ,&Invmass);    
  es->SetBranchAddress("nJet"           ,&nJet);
  es->SetBranchAddress("HT"             ,&HT);
	es->SetBranchAddress("sigMET",         &sigMET);

	TFile xSecFile("cross/susyCrossSection.root");
	TH1D *p_crosssection_tchiwg = (TH1D*)xSecFile.Get("p_charginoSec");
	TH1D *p_crosssection_t5wg   = (TH1D*)xSecFile.Get("p_gluinoxSec");
	TH1D *p_crosssection_t6wg   = (TH1D*)xSecFile.Get("p_squarkxSec");
	TH2D *p_T5WGMASS         = new TH2D("SUSYMass","",27, 775.0, 2125.0, 80, 12.5, 2012.5);
	TH2D *h_chan_rate_nom[NBIN]; 
	TH2D *h_chan_rate_fake[NBIN];
//	TCanvas *can[NBIN];
	std::ostringstream histname; 
	for(unsigned i(0); i < NBIN; i++){
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_nom";
		h_chan_rate_nom[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
		histname.str("");
		histname << "h_chan" << i+1 << "_rate_fake";
		h_chan_rate_fake[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),27, 775.0, 2125.0, 80, 12.5, 2012.5);
//		histname.str("");
//		histname << "can" << i+1;
//    can[i] = new TCanvas(histname.str().c_str(),histname.str().c_str(), 600,600);
	}


	TH1D *h1=new TH1D("h1","",100,0,2000);
	TH1D *h2=new TH1D("h2","",100,0,2000);
	for(unsigned ievt(0); ievt < es->GetEntries(); ievt++){
		es->GetEntry(ievt);

		double NLSPMass(0);
		if(Mchargino > 0)NLSPMass = Mchargino;
		else if(Mneutralino > 0)NLSPMass = Mneutralino;
		if(NLSPMass <= 0)continue;	
		p_T5WGMASS->Fill(Mgluino, NLSPMass);

		int SigBinIndex(-1);
		SigBinIndex = Bin.findSignalBin(sigMET, HT, recoPhotonEt);
    if(recoPhotonEt < 35 || fabs(recoPhotonEta) > 1.4442 || recoEleEt < 25 || ElePassID < 1)continue;
    if(dRPhoEle < 0.8 || Invmass < 100)continue;
    if(MT_ < 100)continue;
		h1->Fill(recoPhotonEt);
 		if(PhoR9 < 0.5)continue;
 		if(PhoHoverE > 0.0597)continue;
    if(PhoNeuIso > 10.910 + 0.0148*recoPhotonEt + 0.000017*recoPhotonEt*recoPhotonEt)continue;
    if(PhoPhoIso > 3.630 + 0.0047*recoPhotonEt)continue;
		if(PhoSigma < 0.01031 && PhoChIso < 1.29)h2->Fill(recoPhotonEt);
		if(SigBinIndex >=0){
 			if((PhoSigma > 0.01031 || PhoChIso > 1.29) && PhoChIso < 15)h_chan_rate_fake[SigBinIndex]->Fill(Mgluino, NLSPMass);
      else if(PhoSigma < 0.01031 && PhoChIso < 1.29)h_chan_rate_nom[SigBinIndex]->Fill(Mgluino, NLSPMass);
		} 
	}
	h2->Divide(h1);
	h2->Draw();

	for(unsigned i(1); i<100; i++)std::cout << i << " " << h2->GetBinContent(i) << std::endl;
//	double lowfake[18];
//	double highfake[18];
//
//	for(unsigned ih(0); ih < NBIN; ih++){
//		std::cout << "chan " << ih << std::endl;
//    lowfake[ih] = 100;
//    highfake[ih] = 0;
//		for(unsigned i(1); i < h_chan_rate_nom[ih]->GetXaxis()->GetNbins() + 1; i++){
//			for(unsigned j(1); j < h_chan_rate_nom[ih]->GetYaxis()->GetNbins() + 1; j++){
//				if(p_T5WGMASS->GetXaxis()->GetBinCenter(i) < 1400)continue;
//				if(p_T5WGMASS->GetBinContent(i,j) > 1000){
//					float noe = p_T5WGMASS->GetBinContent(i,j);
//					float sparticleMass = p_T5WGMASS->GetXaxis()->GetBinCenter(i);
//					float crosssection = p_crosssection_t5wg->GetBinContent( p_crosssection_t5wg->FindBin(sparticleMass) );
//
//					h_chan_rate_nom[ih]->SetBinContent(i,j, h_chan_rate_nom[ih]->GetBinContent(i,j)*35.9*crosssection/noe);
//					h_chan_rate_fake[ih]->SetBinContent(i,j, h_chan_rate_fake[ih]->GetBinContent(i,j)*35.9*crosssection/noe);
//					std::cout << int(p_T5WGMASS->GetXaxis()->GetBinCenter(i)) << " " << int(p_T5WGMASS->GetYaxis()->GetBinCenter(j)) << " " << h_chan_rate_fake[ih]->GetBinContent(i,j) << std::endl;
//					if(h_chan_rate_fake[ih]->GetBinContent(i,j) < lowfake[ih])lowfake[ih] = h_chan_rate_fake[ih]->GetBinContent(i,j);
//          if(h_chan_rate_fake[ih]->GetBinContent(i,j) > highfake[ih])highfake[ih] = h_chan_rate_fake[ih]->GetBinContent(i,j);
//        }
//				else{
//					h_chan_rate_nom[ih]->SetBinContent(i,j, -1); 
//					h_chan_rate_fake[ih]->SetBinContent(i,j,-1); 
//				}
//			}
//		}
//		std::cout << "low " << lowfake[ih] << " high " << highfake[ih] << std::endl;
//		std::cout << std::endl;
//	}

}
