#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TF3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TChain.h"
#include "TSystem.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TGraphErrors.h"

#include "../../include/analysis_rawData.h"
#include "../../include/analysis_photon.h"
#include "../../include/analysis_muon.h"
#include "../../include/analysis_ele.h"
#include "../../include/analysis_mcData.h"
#include "../../include/analysis_tools.h"
#include "analysis_mcfakes.h"
#include "../../include/analysis_scalefactor.h"
#include "../../include/tdrstyle.C"

#define NTOY 100

bool isElectron(int PID, int momID){
   bool isEle;
   if(fabs(PID) == 11){ 
	   switch(momID){
	     case 1: isEle = true; break;
	     case 2: isEle = true; break;
	     case 3: isEle = true; break;
	     case 4: isEle = true; break;
	     case 5: isEle = true; break;
	     case 6: isEle = true; break;
	     case 21: isEle = true; break;
	     case 23: isEle = true; break;
			 case 24: isEle = true; break;
	     default: isEle = false; break;
	   }
  }
  else isEle = false;

  return isEle;
}
Double_t fakerate_ptDependence(Double_t *x, Double_t *par)
{
  double slope = par[0];
  double constant = par[1];
  double index = par[2];
  double coeff = par[3]; 

  double pt = TMath::Max(x[0],0.000001);

   double arg = 0;
   arg = slope*pt + constant; 
   double fitval = pow(arg, index)*coeff; 
   return fitval;
}

void analysis_rareBkg(int ichannel, double inputmetlow, double inputmetup){//main  

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");
	setTDRStyle();
	gStyle->SetOptStat(0);
	gStyle->SetErrorX(0.5);
	gStyle->SetTitleX(0.5);
	esfScaleFactor  objectESF;

  int channelType = ichannel; // eg = 1; mg =2;
	double met_lowcut = inputmetlow;
	double met_upcut  = inputmetup;

	std::ifstream elefake_file("../eleFakePho/fitMC/EleFakeRate-ByPtVtx.txt");
	double scalefactor(0);
	double ptslope(0);
	double ptconstant(0);
	double ptindex(0);
	double ptcoeff(0);
	double vtxconst(0);
	double vtxslope(0);
	std::string variabletype;
	double variablevalue;
	if(elefake_file.is_open()){
  	for(int i(0); i<7; i++){ 
			elefake_file >> variabletype >> variablevalue; 
			if(variabletype.find("scalefactor")!=std::string::npos)scalefactor = variablevalue;
			else if(variabletype.find("ptslope")!=std::string::npos)ptslope = variablevalue;
			else if(variabletype.find("ptconstant")!=std::string::npos)ptconstant = variablevalue;
			else if(variabletype.find("ptindex")!=std::string::npos)ptindex = variablevalue;
			else if(variabletype.find("ptcoeff")!=std::string::npos)ptcoeff = variablevalue;
			else if(variabletype.find("vtxconst")!=std::string::npos)vtxconst = variablevalue;
			else if(variabletype.find("vtxslope")!=std::string::npos)vtxslope = variablevalue;
	  }
	}
	elefake_file.close();
	TF3 h_nominal_fakerate("h_nominal_fakerate", fakerate_func,10,1000,0,100,0,1.5,7);
	h_nominal_fakerate.SetParameters(scalefactor, ptslope, ptconstant, ptindex, 1, vtxconst, vtxslope);
	TF1 *f1 = new TF1("f1", fakerate_ptDependence,35,1000,4);
	f1->SetParameters(ptslope, ptconstant, ptindex, ptcoeff);

	TF3 *h_toymc_fakerate[NTOY];
	std::ostringstream funcname;
	std::ifstream elefake_toyfile("../eleFakePho/fitMC/ToyFakeRate.txt");
	if(elefake_toyfile.is_open()){
  	for(int i(0); i<NTOY; i++){ 
			elefake_toyfile >> scalefactor >> ptslope >> ptconstant >> ptindex >>  vtxconst >> vtxslope;
			funcname.str("");
			funcname << "h_toymc_fakerate_" << i;
			h_toymc_fakerate[i] = new TF3(funcname.str().c_str(), fakerate_func,10,1000,0,100,0,1.5,7);
			h_toymc_fakerate[i]->SetParameters(scalefactor, ptslope, ptconstant, ptindex, 1, vtxconst, vtxslope); 
	  }
	}
	elefake_toyfile.close();

	//*********** histo list **********************//
	Double_t plotEtBins[]={35,40,50,60,70,80,90,100,115,135,160,200,300};
	Double_t plotMETBins[]={0,10,20,30,40,50,60,70,80,90,100,110,125,140,160,200,400};
	TH1F *p_PhoEt = new TH1F("p_PhoEt","Photon p_{T};;Events",12,plotEtBins);
	TH1F *p_PhoEta = new TH1F("p_PhoEta","#gamma #eta; #eta;",60,-3,3);
	Double_t plotPtBins[]={25,30,35,40,50,60,70,80,90,100,115,135,160,200,300};
	TH1F *p_LepPt = new TH1F("p_LepPt","electron p_{T};;Events",14,plotPtBins);
	TH1F *p_LepEta = new TH1F("p_LepEta","p_LepEta",60,-3,3);
	TH1F *p_MET = new TH1F("p_MET","MET;MET (GeV);Events",16,plotMETBins);
	TH1F *p_Mt = new TH1F("p_Mt","M_{T};;Events",16,plotMETBins); 
	TH1F *p_HT = new TH1F("p_HT","HT; HT (GeV);",100,0,1200); 
	TH1F *p_dPhiEleMET = new TH1F("p_dPhiEleMET","dPhiEleMET",64,-3.2,3.2); 

	TH1F *proxy_PhoEt = new TH1F("proxy_PhoEt","Photon p_{T};;Events",12,plotEtBins);
	TH1F *proxy_PhoEta = new TH1F("proxy_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1F *proxy_LepPt = new TH1F("proxy_LepPt","electron p_{T};;Events",14,plotPtBins);
	TH1F *proxy_LepEta = new TH1F("proxy_LepEta","proxy_LepEta",60,-3,3);
	TH1F *proxy_MET = new TH1F("proxy_MET","MET;MET (GeV);Events",16,plotMETBins);
	TH1F *proxy_Mt = new TH1F("proxy_Mt","M_{T};;Events",16,plotMETBins); 
	TH1F *proxy_HT = new TH1F("proxy_HT","HT; HT (GeV);",100,0,1200); 
	TH1F *proxy_dPhiEleMET = new TH1F("proxy_dPhiEleMET","dPhiEleMET",64,-3.2,3.2); 
	TH1F *proxyEWK_PhoEt = new TH1F("proxyEWK_PhoEt","Photon p_{T};;Events",12,plotEtBins);
	TH1F *proxyEWK_PhoEta = new TH1F("proxyEWK_PhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1F *proxyEWK_LepPt = new TH1F("proxyEWK_LepPt","electron p_{T};;Events",14,plotPtBins);
	TH1F *proxyEWK_LepEta = new TH1F("proxyEWK_LepEta","proxyEWK_LepEta",60,-3,3);
	TH1F *proxyEWK_MET = new TH1F("proxyEWK_MET","MET;MET (GeV);Events",16,plotMETBins);
	TH1F *proxyEWK_Mt = new TH1F("proxyEWK_Mt","M_{T};;Events",16,plotMETBins); 
	TH1F *proxyEWK_HT = new TH1F("proxyEWK_HT","HT; HT (GeV);",100,0,1200); 
	TH1F *proxyEWK_dPhiEleMET = new TH1F("proxyEWK_dPhiEleMET","dPhiEleMET",64,-3.2,3.2); 
	TH1F *toy_PhoEt[NTOY];
	TH1F *toy_LepPt[NTOY];
	TH1F *toy_MET[NTOY];
	TH1F *toy_Mt[NTOY];
	TH1F *toy_dPhiEleMET[NTOY];
	std::stringstream histname;
	for(unsigned ih(0); ih < NTOY; ih++){
		histname.str("");
		histname << "toy_PhoEt_ " << ih;
		toy_PhoEt[ih] = new TH1F(histname.str().c_str(), histname.str().c_str(),12,plotEtBins);
		histname.str("");
		histname << "toy_LepPt_" << ih;
		toy_LepPt[ih] = new TH1F(histname.str().c_str(), histname.str().c_str(),14,plotPtBins);
		histname.str("");
		histname << "toy_MET_ " << ih;
		toy_MET[ih] = new TH1F(histname.str().c_str(), histname.str().c_str(),16,plotMETBins);
		histname.str("");
		histname << "toy_Mt_ " << ih;
		toy_Mt[ih] = new TH1F(histname.str().c_str(), histname.str().c_str(),16,plotMETBins);
		histname.str("");
		histname << "toy_eledPhiEleMET_" << ih;
		toy_dPhiEleMET[ih] = new TH1F(histname.str().c_str(), histname.str().c_str(),64,-3.2,3.2);
	}
	TGraphErrors *proxyError_PhoEt = new TGraphErrors(12);
	TGraphErrors *proxyError_PhoEta = new TGraphErrors(60);
	TGraphErrors *proxyError_LepPt = new TGraphErrors(14);
	TGraphErrors *proxyError_LepEta = new TGraphErrors(60);
	TGraphErrors *proxyError_MET = new TGraphErrors(18);
	TGraphErrors *proxyError_Mt = new TGraphErrors(40);
	TGraphErrors *proxyError_HT = new TGraphErrors(100);
	TGraphErrors *proxyError_dPhiEleMET = new TGraphErrors(64);
	TGraphErrors *ratio_PhoEt = new TGraphErrors(12);
	TGraphErrors *ratio_PhoEta = new TGraphErrors(60);
	TGraphErrors *ratio_LepPt = new TGraphErrors(14);
	TGraphErrors *ratio_LepEta = new TGraphErrors(60);
	TGraphErrors *ratio_MET = new TGraphErrors(18);
	TGraphErrors *ratio_Mt = new TGraphErrors(40);
	TGraphErrors *ratio_HT = new TGraphErrors(100);
	TGraphErrors *ratio_dPhiEleMET = new TGraphErrors(64);
	TGraphErrors *ratioError_PhoEt = new TGraphErrors(12);
	TGraphErrors *ratioError_PhoEta = new TGraphErrors(60);
	TGraphErrors *ratioError_LepPt = new TGraphErrors(14);
	TGraphErrors *ratioError_LepEta = new TGraphErrors(60);
	TGraphErrors *ratioError_MET = new TGraphErrors(18);
	TGraphErrors *ratioError_Mt = new TGraphErrors(40);
	TGraphErrors *ratioError_HT = new TGraphErrors(100);
	TGraphErrors *ratioError_dPhiEleMET = new TGraphErrors(64);
// ********  MC *************************//
  TChain *mctree = new TChain("signalTree","signalTree");
  mctree->Add("/uscms_data/d3/mengleis/closureTree_elefake_DY.root");
  mctree->Add("/uscms_data/d3/mengleis/closureTree_elefake_TT.root");
  mctree->Add("/uscms_data/d3/mengleis/closureTree_elefake_WW.root");
	float crosssection(0);
	float ntotalevent(0);
  float phoEt(0);
  float phoEta(0);
  float phoPhi(0);
  float lepPt(0);
  float lepEta(0);
  float lepPhi(0);
  float sigMT(0);
  float sigMET(0);
  float sigMETPhi(0);
  float dPhiLepMET(0);
  int   nVertex(0);
  float dRPhoLep(0);
  float HT(0);
  float nJet(0);
  std::vector<int> *mcPID=0;
	std::vector<float> *mcEta=0;
	std::vector<float> *mcPhi=0;
  std::vector<float> *mcPt=0;
  std::vector<int> *mcMomPID=0;
  std::vector<int> *mcGMomPID=0;

	mctree->SetBranchAddress("crosssection",&crosssection);
	mctree->SetBranchAddress("ntotalevent", &ntotalevent);
  mctree->SetBranchAddress("phoEt",     &phoEt);
  mctree->SetBranchAddress("phoEta",    &phoEta);
  mctree->SetBranchAddress("phoPhi",    &phoPhi);
  mctree->SetBranchAddress("lepPt",     &lepPt);
  mctree->SetBranchAddress("lepEta",    &lepEta);
  mctree->SetBranchAddress("lepPhi",    &lepPhi);
  mctree->SetBranchAddress("sigMT",     &sigMT);
  mctree->SetBranchAddress("sigMET",    &sigMET);
  mctree->SetBranchAddress("sigMETPhi", &sigMETPhi);
  mctree->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
  mctree->SetBranchAddress("nVertex",   &nVertex);
  mctree->SetBranchAddress("dRPhoLep",  &dRPhoLep);
  mctree->SetBranchAddress("HT",        &HT);
  mctree->SetBranchAddress("nJet",      &nJet);
  mctree->SetBranchAddress("mcPID",     &mcPID);
  mctree->SetBranchAddress("mcEta",     &mcEta);
  mctree->SetBranchAddress("mcPhi",     &mcPhi);
  mctree->SetBranchAddress("mcPt",      &mcPt);
  mctree->SetBranchAddress("mcMomPID",  &mcMomPID);
  mctree->SetBranchAddress("mcGMomPID", &mcGMomPID);

	for(unsigned ievt(0); ievt < mctree->GetEntries(); ievt++){
		mctree->GetEntry(ievt);
		float XS_weight = 35.87*1000*crosssection/ntotalevent;
		float weight = XS_weight; 
		if(sigMET > met_upcut || sigMET < met_lowcut)continue;
		if(phoEt < 35 || lepPt < 25)continue;
		if(fabs(phoEta) > 1.4442 || fabs(lepEta) > 2.5)continue;

		bool isZee(false);
		double mindRprobe(0.3);
		unsigned probeIndex(0);
		for(unsigned iMC(0); iMC<mcPID->size(); iMC++){
			double dR = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], phoEta,  phoPhi);
			double dE = fabs((*mcPt)[iMC] - phoEt)/phoEt; 
			if(dR < mindRprobe && dE < 0.1){mindRprobe=dR; probeIndex=iMC;}
		}
		if(mindRprobe < 0.1){
			if(abs((*mcPID)[probeIndex])== 11 || abs((*mcPID)[probeIndex])== 15)isZee = true;
			//isZee = isElectron(fabs((*mcPID)[probeIndex]), fabs((*mcMomPID)[probeIndex]));
			//if(!isZee)std::cout << "  pho " << fabs((*mcPID)[probeIndex]) << " " << fabs((*mcMomPID)[probeIndex]) << std::endl; 
		}
			
		if(isZee){
			p_PhoEt->Fill(phoEt, weight);
			p_PhoEta->Fill(phoEta,weight);
			p_LepPt->Fill(lepPt, weight);
			p_LepEta->Fill(lepEta,weight);
			p_MET->Fill(sigMET, weight);
			p_Mt->Fill(sigMT, weight);
			p_HT->Fill(HT, weight);
			p_dPhiEleMET->Fill(dPhiLepMET, weight);
		}
	}


	//************ Proxy Tree **********************//
	TChain *proxytree = new TChain("proxyTree");
  proxytree->Add("/uscms_data/d3/mengleis/closureTree_elefake_DY.root");
  proxytree->Add("/uscms_data/d3/mengleis/closureTree_elefake_TT.root");
  proxytree->Add("/uscms_data/d3/mengleis/closureTree_elefake_WW.root");

	float proxycrosssection(0);
	float proxyntotalevent(0);
	float proxyphoEt(0);
	float proxyphoEta(0);
	float proxyphoPhi(0);
	float proxylepPt(0);
	float proxylepEta(0);
	float proxylepPhi(0);
	float proxysigMT(0);
	float proxysigMET(0);
	float proxysigMETPhi(0);
	float proxydPhiLepMET(0);
	int   proxynVertex(0);
	float proxydRPhoLep(0);
	float proxyHT(0);
	float proxynJet(0);
	
	proxytree->SetBranchAddress("crosssection",&proxycrosssection);
	proxytree->SetBranchAddress("ntotalevent", &proxyntotalevent);
	proxytree->SetBranchAddress("phoEt",     &proxyphoEt);
	proxytree->SetBranchAddress("phoEta",    &proxyphoEta);
	proxytree->SetBranchAddress("phoPhi",    &proxyphoPhi);
	proxytree->SetBranchAddress("lepPt",     &proxylepPt);
	proxytree->SetBranchAddress("lepEta",    &proxylepEta);
	proxytree->SetBranchAddress("lepPhi",    &proxylepPhi);
	proxytree->SetBranchAddress("sigMT",     &proxysigMT);
	proxytree->SetBranchAddress("sigMET",    &proxysigMET);
	proxytree->SetBranchAddress("sigMETPhi", &proxysigMETPhi);
	proxytree->SetBranchAddress("dPhiLepMET",&proxydPhiLepMET);
	proxytree->SetBranchAddress("nVertex",   &proxynVertex);
	proxytree->SetBranchAddress("dRPhoLep",  &proxydRPhoLep);
	proxytree->SetBranchAddress("HT",        &proxyHT);
	proxytree->SetBranchAddress("nJet",      &proxynJet);
 
	for (unsigned ievt(0); ievt<proxytree->GetEntries(); ++ievt){//loop on entries
		proxytree->GetEntry(ievt);
		float XS_weight = 35.87*1000*proxycrosssection/proxyntotalevent;
		if(proxysigMET > met_upcut || proxysigMET < met_lowcut)continue;
		if(proxyphoEt < 35 || proxylepPt < 25)continue;
		if(fabs(proxyphoEta) > 1.4442 || fabs(proxylepEta) > 2.5)continue;

		double w_ele = h_nominal_fakerate(proxyphoEt, proxynVertex, fabs(proxyphoEta));
		proxy_PhoEt->Fill(proxyphoEt,w_ele*XS_weight);
		proxy_PhoEta->Fill(proxyphoEta, w_ele*XS_weight);
		proxy_MET->Fill(proxysigMET, w_ele*XS_weight);
		proxy_Mt->Fill(proxysigMT, w_ele*XS_weight);
		proxy_HT->Fill(proxyHT, w_ele*XS_weight);
		proxy_LepPt->Fill(proxylepPt, w_ele*XS_weight);
		proxy_LepEta->Fill(proxylepEta, w_ele*XS_weight);
		if(proxycrosssection < 5000){
			proxyEWK_PhoEt->Fill(proxyphoEt,w_ele*XS_weight);
			proxyEWK_PhoEta->Fill(proxyphoEta, w_ele*XS_weight);
			proxyEWK_MET->Fill(proxysigMET, w_ele*XS_weight);
			proxyEWK_Mt->Fill(proxysigMT, w_ele*XS_weight);
			proxyEWK_HT->Fill(proxyHT, w_ele*XS_weight);
			proxyEWK_LepPt->Fill(proxylepPt, w_ele*XS_weight);
			proxyEWK_LepEta->Fill(proxylepEta, w_ele*XS_weight);
		}

		for(unsigned it(0); it < NTOY; it++){
			double toy_ele = h_toymc_fakerate[it]->Eval(proxyphoEt,proxynVertex,fabs(proxyphoEta));
			toy_PhoEt[it]->Fill(proxyphoEt,toy_ele*XS_weight);
			toy_MET[it]->Fill(proxysigMET, toy_ele*XS_weight);
			toy_Mt[it]->Fill(proxysigMT, toy_ele*XS_weight);
			toy_LepPt[it]->Fill(proxylepPt, toy_ele*XS_weight);
		}
	}

	for(int ibin(1); ibin < toy_PhoEt[0]->GetSize(); ibin++){
		TH1F *p_toydistribution = new TH1F("p_toydistribution","",1000,p_PhoEt->GetBinContent(ibin)/10, p_PhoEt->GetBinContent(ibin)*10);
		for(unsigned it(0); it < NTOY; it++){
			p_toydistribution->Fill(toy_PhoEt[it]->GetBinContent(ibin));
		}
		p_toydistribution->Fit("gaus");
		double bincenter = proxy_PhoEt->GetBinLowEdge(ibin) + (proxy_PhoEt->GetBinLowEdge(ibin+1)-proxy_PhoEt->GetBinLowEdge(ibin))/2;
		double binerror  = (proxy_PhoEt->GetBinLowEdge(ibin+1)-proxy_PhoEt->GetBinLowEdge(ibin))/2;
		double syserror = p_toydistribution->GetFunction("gaus")->GetParameter(2);
		double globalerror = 0.13*proxy_PhoEt->GetBinContent(ibin);
		double staterror = proxy_PhoEt->GetBinError(ibin);
		double totalerror = sqrt(syserror*syserror + globalerror*globalerror + staterror*staterror);	
		proxyError_PhoEt->SetPoint(ibin-1,bincenter , proxy_PhoEt->GetBinContent(ibin));
		proxyError_PhoEt->SetPointError(ibin-1,binerror, totalerror);
		ratio_PhoEt->SetPoint(ibin-1,bincenter, p_PhoEt->GetBinContent(ibin)/proxy_PhoEt->GetBinContent(ibin));
		ratio_PhoEt->SetPointError(ibin-1,binerror, p_PhoEt->GetBinError(ibin)/p_PhoEt->GetBinContent(ibin));
		ratioError_PhoEt->SetPoint(ibin-1, bincenter , 1); 
		ratioError_PhoEt->SetPointError(ibin-1, binerror, totalerror/p_PhoEt->GetBinContent(ibin));
		delete p_toydistribution;
	}
	for(int ibin(1); ibin < toy_LepPt[0]->GetSize(); ibin++){
		TH1F *p_toydistribution = new TH1F("p_toydistribution","",1000,p_LepPt->GetBinContent(ibin)/10, p_LepPt->GetBinContent(ibin)*10);
		for(unsigned it(0); it < NTOY; it++){
			p_toydistribution->Fill(toy_LepPt[it]->GetBinContent(ibin));
		}
		p_toydistribution->Fit("gaus");
		double bincenter = proxy_LepPt->GetBinLowEdge(ibin) + (proxy_LepPt->GetBinLowEdge(ibin+1)-proxy_LepPt->GetBinLowEdge(ibin))/2;
		double binerror  = (proxy_LepPt->GetBinLowEdge(ibin+1)-proxy_LepPt->GetBinLowEdge(ibin))/2;
		double syserror = p_toydistribution->GetFunction("gaus")->GetParameter(2);
		double globalerror = 0.13*proxy_LepPt->GetBinContent(ibin);
		double staterror = proxy_LepPt->GetBinError(ibin);
		double totalerror = sqrt(syserror*syserror + globalerror*globalerror + staterror*staterror);	
		proxyError_LepPt->SetPoint(ibin-1,bincenter , proxy_LepPt->GetBinContent(ibin));
		proxyError_LepPt->SetPointError(ibin-1,binerror, totalerror);
		ratio_LepPt->SetPoint(ibin-1,bincenter, p_LepPt->GetBinContent(ibin)/proxy_LepPt->GetBinContent(ibin));
		ratio_LepPt->SetPointError(ibin-1,binerror, p_LepPt->GetBinError(ibin)/p_LepPt->GetBinContent(ibin));
		ratioError_LepPt->SetPoint(ibin-1, bincenter , 1); 
		ratioError_LepPt->SetPointError(ibin-1, binerror, totalerror/p_LepPt->GetBinContent(ibin));
		delete p_toydistribution;
	}
	for(int ibin(1); ibin < toy_MET[0]->GetSize(); ibin++){
		TH1F *p_toydistribution = new TH1F("p_toydistribution","",1000,p_MET->GetBinContent(ibin)/10, p_MET->GetBinContent(ibin)*10);
		for(unsigned it(0); it < NTOY; it++){
			p_toydistribution->Fill(toy_MET[it]->GetBinContent(ibin));
		}
		p_toydistribution->Fit("gaus");
		double bincenter = proxy_MET->GetBinLowEdge(ibin) + (proxy_MET->GetBinLowEdge(ibin+1)-proxy_MET->GetBinLowEdge(ibin))/2;
		double binerror  = (proxy_MET->GetBinLowEdge(ibin+1)-proxy_MET->GetBinLowEdge(ibin))/2;
		double syserror = p_toydistribution->GetFunction("gaus")->GetParameter(2);
		double globalerror = 0.13*proxy_MET->GetBinContent(ibin);
		double staterror = proxy_MET->GetBinError(ibin);
		double totalerror = sqrt(syserror*syserror + globalerror*globalerror + staterror*staterror);	
		proxyError_MET->SetPoint(ibin-1,bincenter , proxy_MET->GetBinContent(ibin));
		proxyError_MET->SetPointError(ibin-1,binerror, totalerror);
		ratio_MET->SetPoint(ibin-1,bincenter, p_MET->GetBinContent(ibin)/proxy_MET->GetBinContent(ibin));
		ratio_MET->SetPointError(ibin-1,binerror, p_MET->GetBinError(ibin)/p_MET->GetBinContent(ibin));
		ratioError_MET->SetPoint(ibin-1, bincenter , 1); 
		ratioError_MET->SetPointError(ibin-1, binerror, totalerror/p_MET->GetBinContent(ibin));
		delete p_toydistribution;
	}
	for(int ibin(1); ibin < toy_Mt[0]->GetSize(); ibin++){
		TH1F *p_toydistribution = new TH1F("p_toydistribution","",1000,p_Mt->GetBinContent(ibin)/10, p_Mt->GetBinContent(ibin)*10);
		for(unsigned it(0); it < NTOY; it++){
			p_toydistribution->Fill(toy_Mt[it]->GetBinContent(ibin));
		}
		p_toydistribution->Fit("gaus");
		double bincenter = proxy_Mt->GetBinLowEdge(ibin) + (proxy_Mt->GetBinLowEdge(ibin+1)-proxy_Mt->GetBinLowEdge(ibin))/2;
		double binerror  = (proxy_Mt->GetBinLowEdge(ibin+1)-proxy_Mt->GetBinLowEdge(ibin))/2;
		double syserror = p_toydistribution->GetFunction("gaus")->GetParameter(2);
		double globalerror = 0.13*proxy_Mt->GetBinContent(ibin);
		double staterror = proxy_Mt->GetBinError(ibin);
		double totalerror = sqrt(syserror*syserror + globalerror*globalerror + staterror*staterror);	
		proxyError_Mt->SetPoint(ibin-1,bincenter , proxy_Mt->GetBinContent(ibin));
		proxyError_Mt->SetPointError(ibin-1,binerror, totalerror);
		ratio_Mt->SetPoint(ibin-1,bincenter, p_Mt->GetBinContent(ibin)/proxy_Mt->GetBinContent(ibin));
		ratio_Mt->SetPointError(ibin-1,binerror, p_Mt->GetBinError(ibin)/p_Mt->GetBinContent(ibin));
		ratioError_Mt->SetPoint(ibin-1, bincenter , 1); 
		ratioError_Mt->SetPointError(ibin-1, binerror, totalerror/p_Mt->GetBinContent(ibin));
		delete p_toydistribution;
	}


	TLegend *leg = new TLegend(0.5,0.6,0.85,0.85);
	leg->AddEntry(p_PhoEt, "simulation");
	leg->AddEntry(proxy_PhoEt,"DY e#rightarrow#gamma prediction");
	leg->AddEntry(proxyEWK_PhoEt,"tt/WW e#rightarrow#gamma prediction");
	leg->AddEntry(proxyError_PhoEt,"total uncertainty");

	TCanvas *canEt = new TCanvas("canEt","",600,600);
	canEt->cd();
	TPad *canEt_pad1 = new TPad("canEt_pad1", "pad1", 0, 0.3, 1, 1.0);
	canEt_pad1->SetBottomMargin(0);
	canEt_pad1->Draw();          
	canEt_pad1->cd();          
	canEt_pad1->SetLogy(); 
	p_PhoEt->Sumw2();
	p_PhoEt->SetMarkerStyle(8);
	proxy_PhoEt->SetLineColor(kRed);
	proxy_PhoEt->SetFillColor(kRed);
	proxy_PhoEt->SetMinimum(5);
	proxy_PhoEt->Draw("hist");
	proxyEWK_PhoEt->SetLineColor(kYellow);
	proxyEWK_PhoEt->SetFillColor(kYellow);
	proxyEWK_PhoEt->Draw("hist same");
  proxyError_PhoEt->SetFillColor(kBlack);
  proxyError_PhoEt->SetFillStyle(3005);
	proxyError_PhoEt->Draw("E2 same");
	p_PhoEt->Draw("EP same");
	leg->Draw("same");

	canEt->cd();   
	TPad *canEt_pad2 = new TPad("canEt_pad2", "pad2", 0, 0.05, 1, 0.3);
	canEt_pad2->SetTopMargin(0);
	canEt_pad2->SetBottomMargin(0.3);
	canEt_pad2->Draw();
	canEt_pad2->cd(); 	
	TH1F *dummy_PhoEt = new TH1F("dummy_PhoEt",";p_{T} (GeV);obs/pred",12,plotEtBins);
	dummy_PhoEt->SetMaximum(2);
	dummy_PhoEt->SetMinimum(0);
	dummy_PhoEt->Draw();
  TLine *flatratio_PhoEt = new TLine(35,1,300,1);
	flatratio_PhoEt->Draw("same");
	ratioError_PhoEt->SetFillStyle(3005);
	ratioError_PhoEt->Draw("E2 same");		
	ratio_PhoEt->SetMarkerStyle(8);	
	ratio_PhoEt->Draw("EP same");
	canEt->SaveAs("closure_elefakepho_PhotonEt.pdf");

	TCanvas *canPt = new TCanvas("canPt","",600,600);
	canPt->cd();
	TPad *canPt_pad1 = new TPad("canPt_pad1", "pad1", 0, 0.3, 1, 1.0);
	canPt_pad1->SetBottomMargin(0);
	canPt_pad1->Draw();          
	canPt_pad1->cd();          
	canPt_pad1->SetLogy(); 
	p_LepPt->Sumw2();
	p_LepPt->SetMarkerStyle(8);
	proxy_LepPt->SetLineColor(kRed);
	proxy_LepPt->SetFillColor(kRed);
	proxy_LepPt->SetMinimum(5);
	proxy_LepPt->Draw("hist");
	proxyEWK_LepPt->SetLineColor(kYellow);
	proxyEWK_LepPt->SetFillColor(kYellow);
	proxyEWK_LepPt->Draw("hist same");
  proxyError_LepPt->SetFillColor(kBlack);
  proxyError_LepPt->SetFillStyle(3005);
	proxyError_LepPt->Draw("E2 same");
	p_LepPt->Draw("EP same");
	leg->Draw("same");

	canPt->cd();   
	TPad *canPt_pad2 = new TPad("canPt_pad2", "pad2", 0, 0.05, 1, 0.3);
	canPt_pad2->SetTopMargin(0);
	canPt_pad2->SetBottomMargin(0.3);
	canPt_pad2->Draw();
	canPt_pad2->cd(); 	
	TH1F *dummy_LepPt = new TH1F("dummy_LepPt",";p_{T} (GeV);obs/pred",14,plotPtBins);
	dummy_LepPt->SetMaximum(2);
	dummy_LepPt->SetMinimum(0);
	dummy_LepPt->Draw();
  TLine *flatratio_LepPt = new TLine(35,1,300,1);
	flatratio_LepPt->Draw("same");
	ratioError_LepPt->SetFillStyle(3005);
	ratioError_LepPt->Draw("E2 same");		
	ratio_LepPt->SetMarkerStyle(8);	
	ratio_LepPt->Draw("EP same");
	canPt->SaveAs("closure_elefakepho_LepPt.pdf");

	TCanvas *canMET = new TCanvas("canMET","",600,600);
	canMET->cd();
	TPad *canMET_pad1 = new TPad("canMET_pad1", "pad1", 0, 0.3, 1, 1.0);
	canMET_pad1->SetBottomMargin(0);
	canMET_pad1->Draw();          
	canMET_pad1->cd();          
	canMET_pad1->SetLogy(); 
	p_MET->Sumw2();
	p_MET->SetMarkerStyle(8);
	proxy_MET->SetLineColor(kRed);
	proxy_MET->SetFillColor(kRed);
	proxy_MET->SetMinimum(5);
	proxy_MET->Draw("hist");
	proxyEWK_MET->SetLineColor(kYellow);
	proxyEWK_MET->SetFillColor(kYellow);
	proxyEWK_MET->Draw("hist same");
  proxyError_MET->SetFillColor(kBlack);
  proxyError_MET->SetFillStyle(3005);
	proxyError_MET->Draw("E2 same");
	p_MET->Draw("EP same");
	leg->Draw("same");

	canMET->cd();   
	TPad *canMET_pad2 = new TPad("canMET_pad2", "pad2", 0, 0.05, 1, 0.3);
	canMET_pad2->SetTopMargin(0);
	canMET_pad2->SetBottomMargin(0.3);
	canMET_pad2->Draw();
	canMET_pad2->cd(); 	
	TH1F *dummy_MET = new TH1F("dummy_MET",";E_{T}^{miss} (GeV);obs/pred",16,plotMETBins);
	dummy_MET->SetMaximum(2);
	dummy_MET->SetMinimum(0);
	dummy_MET->Draw();
  TLine *flatratio_MET = new TLine(35,1,400,1);
	flatratio_MET->Draw("same");
	ratioError_MET->SetFillStyle(3005);
	ratioError_MET->Draw("E2 same");		
	ratio_MET->SetMarkerStyle(8);	
	ratio_MET->Draw("EP same");
	canMET->SaveAs("closure_elefakepho_MET.pdf");

	TCanvas *canMt = new TCanvas("canMt","",600,600);
	canMt->cd();
	TPad *canMt_pad1 = new TPad("canMt_pad1", "pad1", 0, 0.3, 1, 1.0);
	canMt_pad1->SetBottomMargin(0);
	canMt_pad1->Draw();          
	canMt_pad1->cd();          
	canMt_pad1->SetLogy(); 
	p_Mt->Sumw2();
	p_Mt->SetMarkerStyle(8);
	proxy_Mt->SetLineColor(kRed);
	proxy_Mt->SetFillColor(kRed);
	proxy_Mt->SetMinimum(5);
	proxy_Mt->Draw("hist");
	proxyEWK_Mt->SetLineColor(kYellow);
	proxyEWK_Mt->SetFillColor(kYellow);
	proxyEWK_Mt->Draw("hist same");
  proxyError_Mt->SetFillColor(kBlack);
  proxyError_Mt->SetFillStyle(3005);
	proxyError_Mt->Draw("E2 same");
	p_Mt->Draw("EP same");
	leg->Draw("same");

	canMt->cd();   
	TPad *canMt_pad2 = new TPad("canMt_pad2", "pad2", 0, 0.05, 1, 0.3);
	canMt_pad2->SetTopMargin(0);
	canMt_pad2->SetBottomMargin(0.3);
	canMt_pad2->Draw();
	canMt_pad2->cd(); 	
	TH1F *dummy_Mt = new TH1F("dummy_Mt",";M_{T}(GeV);obs/pred",16,plotMETBins);
	dummy_Mt->SetMaximum(2);
	dummy_Mt->SetMinimum(0);
	dummy_Mt->Draw();
  TLine *flatratio_Mt = new TLine(35,1,400,1);
	flatratio_Mt->Draw("same");
	ratioError_Mt->SetFillStyle(3005);
	ratioError_Mt->Draw("E2 same");		
	ratio_Mt->SetMarkerStyle(8);	
	ratio_Mt->Draw("EP same");
	canMt->SaveAs("closure_elefakepho_Mt.pdf");
}


