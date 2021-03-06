#include "../analysis_commoncode.h"

void analysis_qcdBkg(){

	SetRunConfig();
	setTDRStyle();
	bool toDeriveScale(false);

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

  int channelType = ichannel; // eg = 1; mg =2;
	double factorQCD(1);
	double factorQCDUP = 1; 

	if(channelType == 1){
		factorQCD = factor_egQCD;
		factorQCDUP = factor_egQCD + factorerror_egQCD;
	}
	else if(channelType == 2){
		factorQCD = factor_mgQCD;
		factorQCDUP = factor_mgQCD + factorerror_mgQCD;
	}

	if(anatype == 0){
		toDeriveScale = true;
		factorQCD = 1;
		factorQCDUP = 1;
	}

	TFile *scaleFile;
	if(channelType == 1)scaleFile = TFile::Open("../script/qcd_eg_scale.root");
	else if(channelType == 2)scaleFile = TFile::Open("../script/qcd_mg_scale.root");
	TH1D *p_scale;
	if(channelType == 1)p_scale = (TH1D*)scaleFile->Get("transfer_factor");
	else if(channelType == 2)p_scale = (TH1D*)scaleFile->Get("transfer_factor");

	//*********** histo list **********************//
	std::ostringstream histname;
	TH1D *p_PhoEt = new TH1D("p_PhoEt","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *p_LepPt = new TH1D("p_LepPt","p_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *p_MET = new TH1D("p_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *p_Mt = new TH1D("p_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins);
	TH1D *p_HT = new TH1D("p_HT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *p_PhoEta = new TH1D("p_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *p_LepEta = new TH1D("p_LepEta","p_LepEta",60,-3,3);
	TH1D *p_dPhiEleMET = new TH1D("p_dPhiEleMET","dPhiEleMET",32,0,3.2); 
	TH1D *p_PU = new TH1D("p_PU","",100,0,100);
	TH1D *p_nJet = new TH1D("p_nJet","p_nJet",10,0,10);
	TH1D *p_nBJet = new TH1D("p_nBJet","p_nBJet",5,0,5);
	TH1D *p_PhoEt_TT = new TH1D("p_PhoEt_TT","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *p_MET_TT = new TH1D("p_MET_TT","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *p_Mt_TT = new TH1D("p_Mt_TT","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins);
	TH1D *p_HT_TT = new TH1D("p_HT_TT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	
	TH1D *normup_PhoEt = new TH1D("normup_PhoEt","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *normup_PhoEta = new TH1D("normup_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *normup_LepPt = new TH1D("normup_LepPt","normup_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *normup_LepEta = new TH1D("normup_LepEta","normup_LepEta",60,-3,3);
	TH1D *normup_MET = new TH1D("normup_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *normup_Mt = new TH1D("normup_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins); 
	TH1D *normup_HT = new TH1D("normup_HT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *normup_dPhiEleMET = new TH1D("normup_dPhiEleMET","dPhiEleMET",32,0,3.2); 

	TH1D *unweight_PhoEt = new TH1D("unweight_PhoEt","#gamma E_{T}; E_{T} (GeV)",nBkgEtBins,bkgEtBins);
	TH1D *unweight_PhoEta = new TH1D("unweight_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1D *unweight_LepPt = new TH1D("unweight_LepPt","unweight_LepPt",nBkgPtBins,bkgPtBins);
	TH1D *unweight_LepEta = new TH1D("unweight_LepEta","unweight_LepEta",60,-3,3);
	TH1D *unweight_MET = new TH1D("unweight_MET","MET; MET (GeV);",nBkgMETBins, bkgMETBins);
	TH1D *unweight_Mt = new TH1D("unweight_Mt","M_{T}; M_{T} (GeV);",nBkgMtBins,bkgMtBins); 
	TH1D *unweight_HT = new TH1D("unweight_HT","HT; HT (GeV);",nBkgHTBins, bkgHTBins); 
	TH1D *unweight_dPhiEleMET = new TH1D("unweight_dPhiEleMET","dPhiEleMET",32,0,3.2); 
// ********** fake lepton tree ************** //
  TChain *fakeEtree = new TChain("fakeLepTree","fakeLepTree");
	if(channelType==1)fakeEtree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_egsignal_DoubleEG_ReMiniAOD_FullEcal.root");
	if(channelType==2)fakeEtree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_mgsignal_MuonEG_FullEcal.root");
  float phoEt(0);
  float phoEta(0);
  float phoPhi(0);
  float lepPt(0);
  float lepEta(0);
  float lepPhi(0);
	float fakeLepMiniIso(0);
	int   fakeLepIsStandardProxy(0);
  float sigMT(0);
  float sigMET(0);
  float sigMETPhi(0);
  float dPhiLepMET(0);
  int   nVertex(0);
  float dRPhoLep(0);
  float HT(0);
  float nJet(0);
	int   nBJet(0); 
 
  fakeEtree->SetBranchAddress("phoEt",     &phoEt);
  fakeEtree->SetBranchAddress("phoEta",    &phoEta);
  fakeEtree->SetBranchAddress("phoPhi",    &phoPhi);
  fakeEtree->SetBranchAddress("lepPt",     &lepPt);
  fakeEtree->SetBranchAddress("lepEta",    &lepEta);
  fakeEtree->SetBranchAddress("lepPhi",    &lepPhi);
  fakeEtree->SetBranchAddress("fakeLepMiniIso", &fakeLepMiniIso);
  fakeEtree->SetBranchAddress("fakeLepIsStandardProxy",&fakeLepIsStandardProxy);
  fakeEtree->SetBranchAddress("sigMT",     &sigMT);
  fakeEtree->SetBranchAddress("sigMET",    &sigMET);
  fakeEtree->SetBranchAddress("sigMETPhi", &sigMETPhi);
  fakeEtree->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
  fakeEtree->SetBranchAddress("nVertex",   &nVertex);
  fakeEtree->SetBranchAddress("dRPhoLep",  &dRPhoLep);
  fakeEtree->SetBranchAddress("HT",        &HT);
  fakeEtree->SetBranchAddress("nJet",      &nJet);
  fakeEtree->SetBranchAddress("nBJet",     &nBJet);

	for(unsigned ievt(0); ievt < fakeEtree->GetEntries(); ievt++){
		fakeEtree->GetEntry(ievt);

		double w_qcd = 0; 
		double w_qcd_up = 0; 
		double w_qcd_unweight = 0;

		if(channelType == 1){
			w_qcd = factorQCD*p_scale->GetBinContent(p_scale->FindBin(lepPt));
			w_qcd_up = factorQCDUP*p_scale->GetBinContent(p_scale->FindBin(lepPt));
			w_qcd_unweight = factorQCD;
		}
		else{
			w_qcd = factorQCD;
			w_qcd_up = factorQCDUP;
			w_qcd_unweight = factorQCD;
		}
	
		p_PU->Fill(nVertex);
		/** cut flow *****/
		if(phoEt < 35 || fabs(phoEta) > 1.4442)continue;
		if(sigMET < lowMET)continue;
		if(highMET > 0 && sigMET > highMET)continue;
		if(sigMT < lowMt)continue;
		if(highMt > 0 && sigMT > highMt)continue;
		if(lepPt < lowPt)continue;
		if(highPt > 0 && lepPt > highPt)continue;

		bool isProxy(false);
		if(channelType==1){if(fakeLepMiniIso < lepIso*0.1)isProxy=true;}
		else if(channelType==2){if((fakeLepMiniIso > 0.2 && fakeLepMiniIso < lepIso*0.1))isProxy=true;}
		if(fakeLepIsStandardProxy == 0)isProxy = false;
		if(!isProxy)continue;

		p_PhoEt->Fill(phoEt, w_qcd);
		p_PhoEta->Fill(phoEta, w_qcd);
		p_LepPt->Fill(lepPt, w_qcd);
		p_LepEta->Fill(lepEta, w_qcd);
		p_MET->Fill(sigMET, w_qcd);
		p_Mt->Fill(sigMT, w_qcd);
		p_HT->Fill(HT, w_qcd);
		p_dPhiEleMET->Fill(fabs(dPhiLepMET), w_qcd);
		p_nJet->Fill(nJet, w_qcd);	
		p_nBJet->Fill(nBJet, w_qcd);

		if(nBJet >= 1){
			p_PhoEt_TT->Fill(phoEt,  w_qcd);
			p_MET_TT->Fill(sigMET,  w_qcd);
			p_Mt_TT->Fill(sigMT,  w_qcd);
			p_HT_TT->Fill(HT,  w_qcd);
		}

		normup_PhoEt->Fill(phoEt, w_qcd_up);
		normup_PhoEta->Fill(phoEta,w_qcd_up);
		normup_LepPt->Fill(lepPt, w_qcd_up);
		normup_LepEta->Fill(lepEta,w_qcd_up);
		normup_MET->Fill(sigMET, w_qcd_up);
		normup_Mt->Fill(sigMT, w_qcd_up);
		normup_HT->Fill(HT, w_qcd_up);
		normup_dPhiEleMET->Fill(fabs(dPhiLepMET), w_qcd_up);

		unweight_PhoEt->Fill(phoEt, w_qcd_unweight);
		unweight_PhoEta->Fill(phoEta,w_qcd_unweight);
		unweight_LepPt->Fill(lepPt, w_qcd_unweight);
		unweight_LepEta->Fill(lepEta,w_qcd_unweight);
		unweight_MET->Fill(sigMET, w_qcd_unweight);
		unweight_Mt->Fill(sigMT, w_qcd_unweight);
		unweight_HT->Fill(HT, w_qcd_unweight);
		unweight_dPhiEleMET->Fill(fabs(dPhiLepMET), w_qcd_unweight);
	}

	for(int ibin(1); ibin < p_PhoEt->GetSize(); ibin++){
		double syserror(0);
		syserror += p_PhoEt->GetBinError(ibin)* p_PhoEt->GetBinError(ibin);
		syserror += pow((normup_PhoEt->GetBinContent(ibin)-p_PhoEt->GetBinContent(ibin)),2);
		syserror += pow((unweight_PhoEt->GetBinContent(ibin)-p_PhoEt->GetBinContent(ibin)),2);
		p_PhoEt->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_LepPt->GetSize(); ibin++){
		double syserror(0);
		syserror += p_LepPt->GetBinError(ibin)* p_LepPt->GetBinError(ibin);
		syserror += pow((normup_LepPt->GetBinContent(ibin)-p_LepPt->GetBinContent(ibin)),2);
		syserror += pow((unweight_LepPt->GetBinContent(ibin)-p_LepPt->GetBinContent(ibin)),2);
		p_LepPt->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_MET->GetSize(); ibin++){
		double syserror(0);
		syserror += p_MET->GetBinError(ibin)* p_MET->GetBinError(ibin);
		syserror += pow((normup_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)),2);
		syserror += pow((unweight_MET->GetBinContent(ibin)-p_MET->GetBinContent(ibin)),2);
		p_MET->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_Mt->GetSize(); ibin++){
		double syserror(0);
		syserror += p_Mt->GetBinError(ibin)* p_Mt->GetBinError(ibin);
		syserror += pow((normup_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)),2);
		syserror += pow((unweight_Mt->GetBinContent(ibin)-p_Mt->GetBinContent(ibin)),2);
		p_Mt->SetBinError(ibin,sqrt(syserror));
	}	
	for(int ibin(1); ibin < p_HT->GetSize(); ibin++){
		double syserror(0);
		syserror += p_HT->GetBinError(ibin)* p_HT->GetBinError(ibin);
		syserror += pow((normup_HT->GetBinContent(ibin)-p_HT->GetBinContent(ibin)),2);
		syserror += pow((unweight_HT->GetBinContent(ibin)-p_HT->GetBinContent(ibin)),2);
		p_HT->SetBinError(ibin,sqrt(syserror));
	}	

	std::ostringstream outputname;
	switch(anatype){
		case 0: outputname << "controlTree_";break;
		case 1: outputname << "bkgTree_";break;	
		case 2: outputname << "validTree_"; break;
		case 3: outputname << "signalTree_"; break;
	}
	if(channelType==1)outputname << "egamma_qcd";
	else if(channelType==2)outputname << "mg_qcd";
	if(anatype ==0)outputname << "_met" << lowMET <<"_" << highMET << "_pt" << lowPt << "_" << highPt << "_iso" << lepIso;
	outputname << ".root";

	TFile *outputfile = TFile::Open(outputname.str().c_str(),"RECREATE");
	outputfile->cd();
	p_PhoEt->Write();
	p_PhoEta->Write();
	p_LepPt->Write();
	p_LepEta->Write();
	p_MET->Write();
	p_Mt->Write();
	p_HT->Write();
	p_dPhiEleMET->Write();
	p_PU->Write();
	p_nJet->Write();
	p_nBJet->Write();
	p_PhoEt_TT->Write();
	p_MET_TT->Write();
	p_Mt_TT->Write();
	p_HT_TT->Write();
	outputfile->Write();
	outputfile->Close();
}


