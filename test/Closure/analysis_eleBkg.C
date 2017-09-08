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

#include "../../include/analysis_rawData.h"
#include "../../include/analysis_photon.h"
#include "../../include/analysis_muon.h"
#include "../../include/analysis_ele.h"
#include "../../include/analysis_mcData.h"
#include "../../include/analysis_tools.h"
#include "../../include/analysis_fakes.h"

#define NTOY 10

void analysis_eleBkg(int ichannel, double inputmetlow, double inputmetup){//main  

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");
  int channelType = ichannel; // eg = 1; mg =2;
	double met_lowcut = inputmetlow;
	double met_upcut  = inputmetup;

  /**********************************/
	/*	double normfactor = par[0]; 	*/  
  /*	double slope = par[1];				*/
  /*	double constant = par[2];			*/
  /*	double index = par[3];				*/
  /*	double coeff = par[4]; 				*/
	/*	double vtx_constant = par[5];	*/
	/*	double vtx_slope = par[6];		*/
  /**********************************/
	std::ifstream elefake_file("validateresult/EleFakeRate-ByPtVtx.txt");
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
	h_nominal_fakerate.SetParameters(scalefactor, ptslope, ptconstant, ptindex, ptcoeff, vtxconst, vtxslope);

	TF3 *h_toymc_fakerate[NTOY];
	std::ostringstream funcname;
	std::ifstream elefake_toyfile("validateresult/ToyFakeRate.txt");
	if(elefake_toyfile.is_open()){
  	for(int i(0); i<NTOY; i++){ 
			elefake_toyfile >> scalefactor >> ptslope >> ptconstant >> ptindex >> ptcoeff >> vtxconst >> vtxslope;
			funcname.str("");
			funcname << "h_toymc_fakerate_" << i;
			h_toymc_fakerate[i] = new TF3(funcname.str().c_str(), fakerate_func,10,1000,0,100,0,1.5,7);
			h_toymc_fakerate[i]->SetParameters(scalefactor, ptslope, ptconstant, ptindex, ptcoeff, vtxconst, vtxslope); 
	  }
	}
	elefake_toyfile.close();


	//*********** histo list **********************//
	std::ostringstream histname;
	Double_t plotEtBins[]={35,40,45,50,55,60,65,70,75,80, 85,90,95,100,105,110,115,120,125,130, 135,140,146,152,158,164,170,177,184,192, 200,208,216,224,232,240,250,260,275,290, 305,325,345,370,400};
	TH1F *p_PhoEt = new TH1F("p_PhoEt","#gamma E_{T}; E_{T} (GeV)",44,plotEtBins);
	TH1F *p_PhoEta = new TH1F("p_PhoEta","#gamma #eta; #eta;",60,-3,3);
	Double_t plotPtBins[]={25,30,35,40,45,50,55,60,65,70,75,80, 85,90,95,100,105,110,115,120,125,130, 135,140,146,152,158,164,170,177,184,192, 200,208,216,224,232,240,250,260,275,290, 305,325,345,370,400};
	TH1F *p_LepPt = new TH1F("p_LepPt","p_LepPt",46,plotPtBins);
	TH1F *p_LepEta = new TH1F("p_LepEta","p_LepEta",60,-3,3);
	TH1F *p_HT = new TH1F("p_HT","p_HT",100,0,1200);
	TH1F *p_MET = new TH1F("p_MET","MET; MET (GeV);",100,0,100);
	TH1F *p_Mt = new TH1F("p_Mt","M_{T}; M_{T} (GeV);",200,0,400); 
	TH1F *p_dPhiEleMET = new TH1F("p_dPhiEleMET","dPhiEleMET",64,-3.2,3.2); 
	TH1F *p_PU = new TH1F("p_PU","",100,0,100);
	TH1F *toy_PhoEt[NTOY];
	TH1F *toy_LepPt[NTOY];
	TH1F *toy_HT[NTOY];
	TH1F *toy_MET[NTOY];
	TH1F *toy_Mt[NTOY];
	TH1F *toy_dPhiEleMET[NTOY];
	for(unsigned ih(0); ih < NTOY; ih++){
		histname.str("");
		histname << "toy_PhoEt_ " << ih;
		toy_PhoEt[ih] = new TH1F(histname.str().c_str(), histname.str().c_str(),44,plotEtBins);
		histname.str("");
		histname << "toy_LepPt_" << ih;
		toy_LepPt[ih] = new TH1F(histname.str().c_str(), histname.str().c_str(),46,plotPtBins);
		histname.str("");
		histname << "toy_MET_ " << ih;
		toy_MET[ih] = new TH1F(histname.str().c_str(), histname.str().c_str(),100,0,100);
		histname.str("");
		histname << "toy_Mt_ " << ih;
		toy_Mt[ih] = new TH1F(histname.str().c_str(), histname.str().c_str(),200,0,400);
		histname.str("");
		histname << "toy_HT_ " << ih;
		toy_HT[ih] = new TH1F(histname.str().c_str(), histname.str().c_str(),100,0,1200);
		histname.str("");
		histname << "toy_eledPhiEleMET_" << ih;
		toy_dPhiEleMET[ih] = new TH1F(histname.str().c_str(), histname.str().c_str(),64,-3.2,3.2);
	}
	TProfile *pro_PhoEt = new TProfile("pro_PhoEt", "pro_PhoEt" , 44,plotEtBins);
	TProfile *pro_LepPt = new TProfile("pro_LepPt", "pro_LepPt" , 46,plotPtBins);
	TProfile *pro_MET   = new TProfile("pro_MET", "pro_MET",100,0,100);
	TProfile *pro_Mt    = new TProfile("pro_Mt", "pro_Mt",200,0,400);
	TProfile *pro_HT    = new TProfile("pro_HT", "pro_Mt",100,0,1200);
	TProfile *pro_dPhiEleMET = new TProfile("pro_dPhiEleMET","pro_dPhiEleMET",64,-3.2,3.2 );
	pro_PhoEt->SetErrorOption("s");	
	pro_LepPt->SetErrorOption("s");	
	pro_MET->SetErrorOption("s");	
	pro_Mt->SetErrorOption("s");
	pro_HT->SetErrorOption("s");	
	pro_dPhiEleMET->SetErrorOption("s");	
	
	//************ Proxy Tree **********************//
	TChain *proxytree = new TChain("proxyTree");
	if(channelType==1)proxytree->Add("/uscms_data/d3/mengleis/ReMiniAOD/resTree_egsignal_DoubleEG_ReMiniAOD.root");
	if(channelType==2)proxytree->Add("/uscms_data/d3/mengleis/Rereco/resTree_mgsignal_MuonEG_FebReminiAOD.root");

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
	
	proxytree->SetBranchAddress("phoEt",     &phoEt);
	proxytree->SetBranchAddress("phoEta",    &phoEta);
	proxytree->SetBranchAddress("phoPhi",    &phoPhi);
	proxytree->SetBranchAddress("lepPt",     &lepPt);
	proxytree->SetBranchAddress("lepEta",    &lepEta);
	proxytree->SetBranchAddress("lepPhi",    &lepPhi);
	proxytree->SetBranchAddress("sigMT",     &sigMT);
	proxytree->SetBranchAddress("sigMET",    &sigMET);
	proxytree->SetBranchAddress("sigMETPhi", &sigMETPhi);
	proxytree->SetBranchAddress("dPhiLepMET",&dPhiLepMET);
	proxytree->SetBranchAddress("nVertex",   &nVertex);
	proxytree->SetBranchAddress("dRPhoLep",  &dRPhoLep);
	proxytree->SetBranchAddress("HT",        &HT);
	proxytree->SetBranchAddress("nJet",      &nJet);
 
	for (unsigned ievt(0); ievt<proxytree->GetEntries(); ++ievt){//loop on entries
		proxytree->GetEntry(ievt);
		p_PU->Fill(nVertex);
		if(sigMET > met_upcut || sigMET < met_lowcut)continue;
		if(phoEt < 35 || lepPt < 25)continue;
		if(fabs(phoEta) > 1.4442 || fabs(lepEta) > 2.5)continue;

		double w_ele = h_nominal_fakerate(phoEt, nVertex, fabs(phoEta));
		p_PhoEt->Fill(phoEt,w_ele);
		p_PhoEta->Fill(phoEta, w_ele);
		p_MET->Fill(sigMET, w_ele);
		p_Mt->Fill(sigMT, w_ele);
		p_HT->Fill(HT, w_ele);
		p_LepPt->Fill(lepPt, w_ele);
		p_LepEta->Fill(lepEta, w_ele);
		p_dPhiEleMET->Fill(dPhiLepMET, w_ele);

		for(unsigned it(0); it < NTOY; it++){
			double toy_ele = h_toymc_fakerate[it]->Eval(phoEt,nVertex,fabs(phoEta));
			if(phoEt > 200)toy_ele = h_toymc_fakerate[it]->Eval(200,nVertex,fabs(phoEta));
			toy_PhoEt[it]->Fill(phoEt,toy_ele);
			toy_MET[it]->Fill(sigMET, toy_ele);
			toy_Mt[it]->Fill(sigMT, toy_ele);
			toy_HT[it]->Fill(HT, toy_ele);
			toy_LepPt[it]->Fill(lepPt, toy_ele);
			toy_dPhiEleMET[it]->Fill(dPhiLepMET, toy_ele);
		}
	}

	for(unsigned it(0); it < NTOY; it++){
	  for(int ibin(1); ibin < toy_PhoEt[0]->GetSize(); ibin++)	
			pro_PhoEt->Fill(toy_PhoEt[it]->GetBinCenter(ibin), toy_PhoEt[it]->GetBinContent(ibin));
	  for(int ibin(1); ibin < toy_LepPt[0]->GetSize(); ibin++)	
			pro_LepPt->Fill(toy_LepPt[it]->GetBinCenter(ibin), toy_LepPt[it]->GetBinContent(ibin));
	  for(int ibin(1); ibin < toy_MET[0]->GetSize(); ibin++)	
			pro_MET->Fill(toy_MET[it]->GetBinCenter(ibin), toy_MET[it]->GetBinContent(ibin));
	  for(int ibin(1); ibin < toy_Mt[0]->GetSize(); ibin++)	
			pro_Mt->Fill(toy_Mt[it]->GetBinCenter(ibin), toy_Mt[it]->GetBinContent(ibin));
	  for(int ibin(1); ibin < toy_HT[0]->GetSize(); ibin++)	
			pro_HT->Fill(toy_HT[it]->GetBinCenter(ibin), toy_HT[it]->GetBinContent(ibin));
	  for(int ibin(1); ibin < toy_dPhiEleMET[0]->GetSize(); ibin++)	
			pro_dPhiEleMET->Fill(toy_dPhiEleMET[it]->GetBinCenter(ibin), toy_dPhiEleMET[it]->GetBinContent(ibin));
	}
	//for(int ibin(1); ibin < p_PhoEt->GetSize(); ibin++)	
	//	p_PhoEt->SetBinError(pro_PhoEt->GetBinCenter(ibin), pro_PhoEt->GetBinError(ibin));
	//for(int ibin(1); ibin < p_LepPt->GetSize(); ibin++)	
	//	p_LepPt->SetBinError(pro_LepPt->GetBinCenter(ibin), pro_LepPt->GetBinError(ibin));
	//for(int ibin(1); ibin < p_MET->GetSize(); ibin++)	
	//	p_MET->SetBinError(pro_MET->GetBinCenter(ibin), pro_MET->GetBinError(ibin));
	//for(int ibin(1); ibin < p_Mt->GetSize(); ibin++)	
	//	p_Mt->SetBinError(pro_Mt->GetBinCenter(ibin), pro_Mt->GetBinError(ibin));
	//for(int ibin(1); ibin < p_HT->GetSize(); ibin++)	
	//	p_HT->SetBinError(pro_HT->GetBinCenter(ibin), pro_HT->GetBinError(ibin));
	//for(int ibin(1); ibin < p_dPhiEleMET->GetSize(); ibin++)	
	//	p_dPhiEleMET->SetBinError(pro_dPhiEleMET->GetBinCenter(ibin), pro_dPhiEleMET->GetBinError(ibin));
	//	

	std::ostringstream outputname;
	if(channelType==1)outputname << "bkgTree_eg_eleBkg_ReMiniAOD.root";
	else if(channelType==2)outputname << "bkgTree_mg_eleBkg_ReMiniAOD.root";
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
	for(unsigned it(0); it < NTOY; it++){
		toy_PhoEt[it]->Write();
		toy_MET[it]->Write();
		toy_Mt[it]->Write();
		toy_LepPt[it]->Write();
		toy_HT[it]->Write();
		toy_dPhiEleMET[it]->Write();
	}
	outputfile->Write();
	outputfile->Close();
}


