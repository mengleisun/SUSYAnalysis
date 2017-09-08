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
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TPad.h"
#include "TFitResult.h"
#include "TVirtualFitter.h"
#include "TMatrixDSym.h"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooExponential.h"
#include "RooMCStudy.h"
#include "RooChi2MCSModule.h"
#include "RooPlot.h"
#include "TH1.h"
#include "RooFitResult.h"
#include "RooAbsReal.h"
#include "RooMultiVarGaussian.h"

#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_rawData.h"
#include "../../../include/tdrstyle.C"
#include "../../../include/analysis_tools.h"
#include "../../../include/analysis_fakes.h"

#define NTOY 10000
#define NBIN 19
#define REBINSIZE 2

Double_t tmpjetfake_func(Double_t *x, Double_t *par)
{
	double pt_low = x[0] - REBINSIZE/2.0;
	double pt_high = x[0]+ REBINSIZE/2.0;

	double c1 = par[0];
	double c2 = par[1];
	double lamda1 = par[2];
	double lamda2 = par[3];

	double jetfakes_lowedge = c1*exp(lamda1*pt_low)/lamda1 + c2*exp(lamda2*pt_low)/lamda2;
	double jetfakes_highedge =  c1*exp(lamda1*pt_high)/lamda1 + c2*exp(lamda2*pt_high)/lamda2;
	//return (jetfakes_highedge + jetfakes_lowedge)/2.0*REBINSIZE;
	return (jetfakes_highedge - jetfakes_lowedge);
}

//void fitJetFunc(int inpar1, float inpar2){//main   
void fitJetFunc(){
  int channel = 2; // 1 = eg; 2 = mg; 3 = egloose; 4 = mgloose;

	setTDRStyle();
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	gStyle->SetErrorX(0.5);
	gStyle->SetTitleX(0.5);

	 
  TChain *sigtree = new TChain("signalTree");
  //if(channel == 1)sigtree->Add("/uscms_data/d3/mengleis/FullStatusData/resTree_egsignal_DoubleEG_ReMiniAOD.root");
  if(channel == 1)sigtree->Add("/uscms_data/d3/mengleis/resTree_egsignal_DoubleEG_ReMiniAOD_test_2.root");
  else if(channel ==2)sigtree->Add("/uscms_data/d3/mengleis/test/resTree_mgsignal_MuonEG_FebReminiAOD_MiniIso_JES.root");

  TChain *controltree = new TChain("jetTree");
  //if(channel == 1)controltree->Add("/uscms_data/d3/mengleis/FullStatusData/resTree_egsignal_DoubleEG_ReMiniAOD.root");
  if(channel == 1)controltree->Add("/uscms_data/d3/mengleis/resTree_egsignal_DoubleEG_ReMiniAOD_test_2.root");
  else if(channel == 2)controltree->Add("/uscms_data/d3/mengleis/test/resTree_mgsignal_MuonEG_FebReminiAOD_MiniIso_JES.root");

  std::ostringstream outputname;
  if(channel == 1)outputname << "jet-fake-pho_egpt.root";
  else if(channel == 2)outputname << "jet-fake-pho_mgpt.root";
  TFile *outputfile = TFile::Open(outputname.str().c_str(),"RECREATE");
  outputfile->cd();
	TH1F *p_jetbkgPhoEt = new TH1F("p_jetbkgPhoEt","#gamma E_{T}; E_{T} (GeV)",1500,0,1500);
	TH1F *p_jetbkgPhoEta = new TH1F("p_jetbkgPhoEta","#gamma #eta; #eta;",60,-3,3);
	TH1F *p_jetbkgLepPt = new TH1F("p_jetbkgLepPt","p_jetbkgLepPt",1500,0,1500);
	TH1F *p_jetbkgLepEta = new TH1F("p_jetbkgLepEta","p_jetbkgLepEta",60,-3,3);
	TH1F *p_jetbkgdRPhoLep = new TH1F("p_jetbkgdRPhoLep","p_jetbkgdRPhoLep",60,0,6);
	TH1F *p_jetbkgMET = new TH1F("p_jetbkgMET","MET; MET (GeV);",100,0,100);
	TH1F *p_jetbkgMt = new TH1F("p_jetbkgMt","M_{T}; M_{T} (GeV);",200,0,400); 

	std::ifstream jetfake_file("JetFakeRate-ChIso-MuonEG-Aug3.txt");
	//std::ifstream jetfake_file("JetFakeRate-ChIso-DoubleEG-ReMiniAOD-neweleCo.txt");
	double PtBin[NBIN];
	double fracHad[NBIN];
	double fracHadError[NBIN];
	int i(0);
	float pt_lower(0), pt_upper(0);
	float fakerate(0), error(0), systematic(0);
	float truefake; 
	if(jetfake_file.is_open()){
		for(int i(0); i < NBIN; i++){ 
			jetfake_file >> pt_lower >> pt_upper >> fakerate >> error >> systematic >> truefake;
			PtBin[i] = pt_lower;
			fracHad[i]= fakerate;
			fracHadError[i] = error;
		}
		jetfake_file.close(); 
	}

	std::ifstream elefake_file("../../eleFakePho/Aug3/EleFakeRate-ByPtVtx.txt");
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
  /**********************************/
	/*	double normfactor = par[0]; 	*/  
  /*	double slope = par[1];				*/
  /*	double constant = par[2];			*/
  /*	double index = par[3];				*/
  /*	double coeff = par[4]; 				*/
	/*	double vtx_constant = par[5];	*/
	/*	double vtx_slope = par[6];		*/
  /**********************************/
	TF3 f3("f3", fakerate_func,10,1000,0,100,0,1.5,7);
	f3.SetParameters(scalefactor, ptslope, ptconstant, ptindex, 1.0, vtxconst, vtxslope);

  TH1F *p_controlPhoEt = new TH1F("p_controlPhoEt",";p_{T} (GeV);Events",115,35,150);
  TH1F *p_sigPhoEt  = new TH1F("p_sigPhoEt",";p_{T} (GeV);Events",115,35,150);
  TH1F *p_fakesPhoEt = new TH1F("p_fakesPhoEt",";p_{T} (GeV);Events",115,35,150);
  TH1F *p_elebkgPhoEt = new TH1F("p_elebkgPhoEt",";p_{T} (GeV);Events",115,35,150);

  sigtree->Draw("phoEt >> p_sigPhoEt", " phoEt >35 && sigMET < 70 && fabs(phoEta) < 1.5");
  controltree->Draw("phoEt >> p_controlPhoEt", "phoEt > 35 && sigMET < 70 && fabs(phoEta) < 1.5");

  //p_sigPhoEt->SetBinContent(115, p_sigPhoEt->GetBinContent(115) + p_sigPhoEt->GetBinContent(116));
  //p_controlPhoEt->SetBinContent(115, p_controlPhoEt->GetBinContent(115) + p_controlPhoEt->GetBinContent(116));
  
  p_sigPhoEt->Sumw2();
  p_controlPhoEt->Sumw2();
//************ Proxy Tree **********************//
  TChain *proxytree = new TChain("proxyTree");
  //if(channel == 1)proxytree->Add("/uscms_data/d3/mengleis/FullStatusData/resTree_egsignal_DoubleEG_ReMiniAOD.root");
  if(channel == 1)proxytree->Add("/uscms_data/d3/mengleis/resTree_egsignal_DoubleEG_ReMiniAOD_test_2.root");
	else if(channel == 2)proxytree->Add("/uscms_data/d3/mengleis/test/resTree_mgsignal_MuonEG_FebReminiAOD_MiniIso_JES.root");

  float proxyphoEt(0);
  float proxyphoEta(0);
	float proxylepPt(0);
  float proxysigMET(0);
  int   proxynVertex(0);
  
  proxytree->SetBranchAddress("phoEt",     &proxyphoEt);
  proxytree->SetBranchAddress("phoEta",    &proxyphoEta);
  proxytree->SetBranchAddress("lepPt",     &proxylepPt);
  proxytree->SetBranchAddress("sigMET",    &proxysigMET);
  proxytree->SetBranchAddress("nVertex",   &proxynVertex);
 
  for (unsigned ievt(0); ievt<proxytree->GetEntries(); ++ievt){//loop on entries
		proxytree->GetEntry(ievt);
    if(proxysigMET > 70)continue;
		if(proxyphoEt < 35 || fabs(proxyphoEta) > 1.5)continue;
		double w_ele = 1;
		w_ele = f3(proxyphoEt, proxynVertex, fabs(proxyphoEta)); 
		p_elebkgPhoEt->Fill(proxyphoEt,w_ele);
  }
  p_elebkgPhoEt->Sumw2();
  p_sigPhoEt->Add(p_elebkgPhoEt, -1);

  p_sigPhoEt->Sumw2();
  for(unsigned ibin(1); ibin < p_sigPhoEt->GetSize()-1; ibin++){
    double xvalue = p_sigPhoEt->GetBinCenter(ibin);
    double frac(0),fracerror(0);
    for(unsigned i(0); i< NBIN-1; i++) 
      if(xvalue >= PtBin[i] && xvalue < PtBin[i+1]){
        frac = fracHad[i];
        fracerror = fracHadError[i];
      }
    if(xvalue >= PtBin[NBIN-1]){ frac = fracHad[NBIN-1]; fracerror = fracHadError[NBIN-1];}
    
    double binvalue = p_sigPhoEt->GetBinContent(ibin)*frac;
    double binerror = p_sigPhoEt->GetBinError(ibin);  
    if(binvalue == 0)continue;
    double totalerror = sqrt(binvalue*binvalue*fracerror*fracerror + binerror*binerror*frac*frac);
    p_fakesPhoEt->SetBinContent(ibin, binvalue);
    p_fakesPhoEt->SetBinError(ibin, binerror);
  }


Double_t plotPtBins[]={35,40,45,50,55,60,65,70,75,80,85,90,100,110,120,130,140,160,180};
TCanvas *c_pt = new TCanvas("Photon_Pt", "Photon P_{T}",800,800);
c_pt->cd();
TPad *can_pad1 = new TPad("can_pad1", "pad1", 0, 0.05, 1, 1.0);
can_pad1->SetBottomMargin(0.1);
can_pad1->Draw();          
can_pad1->cd();          
gPad->SetLogy();
TH1F *new_controlPhoEt = (TH1F*)p_controlPhoEt->Rebin(REBINSIZE);
TH1F *new_fakesPhoEt = (TH1F*)p_fakesPhoEt->Rebin(REBINSIZE);
new_controlPhoEt->SetMinimum(new_fakesPhoEt->GetBinContent(new_fakesPhoEt->GetSize()-1));
new_fakesPhoEt->SetMinimum(new_fakesPhoEt->GetBinContent(new_fakesPhoEt->GetSize()-1));
new_fakesPhoEt->GetXaxis()->SetTitle("p_{T} (GeV)");
new_controlPhoEt->GetXaxis()->SetTitle("p_{T} (GeV)");
new_controlPhoEt->GetXaxis()->SetTitleOffset(1);
new_controlPhoEt->GetXaxis()->SetTitleSize(20);
new_controlPhoEt->Sumw2();
new_fakesPhoEt->Sumw2();
new_controlPhoEt->GetXaxis()->SetRangeUser(35,145);
new_fakesPhoEt->GetXaxis()->SetRangeUser(35,145);
new_controlPhoEt->Draw("EP");
new_controlPhoEt->SetLineColor(kBlack);
new_controlPhoEt->SetMarkerStyle(20);
new_fakesPhoEt->SetLineColor(kRed);
new_fakesPhoEt->SetMarkerStyle(20);
new_fakesPhoEt->SetMarkerColor(kRed);
new_fakesPhoEt->Draw("EP same");
TLegend *leg =  new TLegend(0.6,0.7,0.85,0.85);
leg->SetFillStyle(0);
leg->SetBorderSize(0);
gStyle->SetLegendFillColor(0);
leg->AddEntry(new_controlPhoEt,"hadron proxies");
leg->AddEntry(new_fakesPhoEt,"fake photons");
leg->Draw("same");

//c_pt->cd();
//TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.25);
//pad2->Draw();
//pad2->cd();
//TH1F *ratio=(TH1F*)new_fakesPhoEt->Clone("transfer factor");
//ratio->Divide(new_controlPhoEt);
//ratio->SetTitle("");
//ratio->GetYaxis()->SetTitle("proxies/fake");
//ratio->GetXaxis()->SetLabelFont(63);
//ratio->GetXaxis()->SetLabelSize(14);
//ratio->GetYaxis()->SetLabelFont(400);
//ratio->GetYaxis()->SetLabelSize(100);
//ratio->Draw();
//

TCanvas* mccan = new TCanvas("mccan","mccan",1200,600) ;
mccan->Divide(2);

//********************   denominator *****************************************************//
new_controlPhoEt->Fit("expo");
TF1 *inifitden = new_controlPhoEt->GetFunction("expo");
double iniLambda_den1 = inifitden->GetParameter(1);
double iniCoeff_den  	= exp(inifitden->GetParameter(0))/2;
double iniLambda_den2 = (log(new_controlPhoEt->GetBinContent(38)- 2*iniCoeff_den*exp(iniLambda_den1*new_controlPhoEt->GetBinCenter(38))) - log(new_controlPhoEt->GetBinContent(56)- 2*iniCoeff_den*exp(iniLambda_den1*new_controlPhoEt->GetBinCenter(56))) )/(new_controlPhoEt->GetBinCenter(38)-new_controlPhoEt->GetBinCenter(56));
TF1 *fitfunc_den= new TF1("fitfunc_den", tmpjetfake_func, 35, 150, 4);
fitfunc_den->SetParameters(iniCoeff_den, iniCoeff_den/10, iniLambda_den1, iniLambda_den2);
new_controlPhoEt->Fit("fitfunc_den","S");
TF1 *fitden = new_controlPhoEt->GetFunction("fitfunc_den");
ofstream myfile;
myfile.open("JetFakeRate-transferfactor-MuonEG-Aug3.txt");

TH1F *ratio=(TH1F*)new_fakesPhoEt->Clone("transfer factor");
ratio->Divide(new_controlPhoEt);
myfile << "den_coeff1 " << fitden->GetParameter(0) << std::endl;
myfile << "den_coeff2 " << fitden->GetParameter(1) << std::endl;
myfile << "den_lambd1 " << fitden->GetParameter(2) << std::endl;
myfile << "den_lambd2 " << fitden->GetParameter(3) << std::endl;

TFitResultPtr rden = new_controlPhoEt->Fit("fitfunc_den","S");
TMatrixDSym covden = rden->GetCovarianceMatrix(); 
rden->Print("V");     
fitfunc_den->Draw("same");
TVectorD muden(4) ;
muden(0) = rden->Parameter(0); 
muden(1) = rden->Parameter(1);
muden(2) = rden->Parameter(2);
muden(3) = rden->Parameter(3);
RooRealVar central_coeff1_den("central_coeff1_den","central_coeff1_den",muden(0)-rden->ParError(0), muden(0)+rden->ParError(0));
RooRealVar central_coeff2_den("central_coeff2_den","central_coeff2_den",muden(1)-rden->ParError(1), muden(1)+rden->ParError(1));
RooRealVar central_lambda1_den("central_lambda1_den","central_lambda1_den",muden(2)-rden->ParError(2),muden(2)+rden->ParError(2));
RooRealVar central_lambda2_den("central_lambda2_den","central_lambda2_den",muden(3)-rden->ParError(3),muden(3)+rden->ParError(3));

RooMultiVarGaussian mvgden("mvgden","mvgden",RooArgList(central_coeff1_den,central_coeff2_den,central_lambda1_den,central_lambda2_den),muden,covden);
RooDataSet* toymcdataden = mvgden.generate(RooArgSet(central_coeff1_den,central_coeff2_den,central_lambda1_den,central_lambda2_den),NTOY);
std::ostringstream modelnameden;
TF1 *gen_den[NTOY];
TH1F *den_upper = new TH1F("den_upper","den_upper",115,35,150);
TH1F *den_lower = new TH1F("den_lower","den_lower",115,35,150);
can_pad1->cd();          
for(unsigned ibin(1); ibin <= 115; ibin++)den_lower->SetBinContent(ibin, fitfunc_den->Eval(den_upper->GetBinCenter(ibin)));
for(int i(0); i<NTOY; i++){
	double data1 = toymcdataden->get(i)->getRealValue("central_coeff1_den");
	double data2 = toymcdataden->get(i)->getRealValue("central_coeff2_den");
	double data3 = toymcdataden->get(i)->getRealValue("central_lambda1_den");
	double data4 = toymcdataden->get(i)->getRealValue("central_lambda2_den");
	modelnameden.str("");
	modelnameden << "gen_den_" << i;
	gen_den[i] = new TF1(modelnameden.str().c_str(), tmpjetfake_func, 35, 150, 4);
	gen_den[i]->SetParameters(data1, data2, data3, data4);
	gen_den[i]->SetLineColorAlpha(kBlue, 0.35);
	//gen_den[i]->Draw("same");
 	for(unsigned ibin(1); ibin <= 115; ibin++){
		double estimated = gen_den[i]->Eval(den_upper->GetBinCenter(ibin));
		if(den_upper->GetBinContent(ibin) < estimated)den_upper->SetBinContent(ibin, estimated);
		if(den_lower->GetBinContent(ibin) > estimated)den_lower->SetBinContent(ibin, estimated);
	}
}

den_upper->Draw("L same");
den_lower->Draw("L same");
// *************************  Numerator ******************************************************************//

new_fakesPhoEt->Fit("expo");
TF1 *inifit = new_fakesPhoEt->GetFunction("expo");
double iniLambda_num1 = inifit->GetParameter(1);
double iniCoeff_num  	= exp(inifit->GetParameter(0))/2;
double iniLambda_num2 = (log(new_fakesPhoEt->GetBinContent(38)- 2*iniCoeff_num*exp(iniLambda_num1*new_fakesPhoEt->GetBinCenter(38))) - log(new_fakesPhoEt->GetBinContent(56)- 2*iniCoeff_num*exp(iniLambda_num1*new_fakesPhoEt->GetBinCenter(56))) )/(new_fakesPhoEt->GetBinCenter(38)-new_fakesPhoEt->GetBinCenter(56));
TF1 *fitfunc_num= new TF1("fitfunc_num", tmpjetfake_func, 35, 150, 4);
//fitfunc_num->SetParameters(iniCoeff_num, iniCoeff_num/10, iniLambda_num1, iniLambda_num2);
if(channel == 1)fitfunc_num->SetParameters(iniCoeff_num, 750, iniLambda_num1, -0.03);
else if(channel == 2)fitfunc_num->SetParameters(iniCoeff_num, 110, iniLambda_num1, -0.82);
TVirtualFitter::SetMaxIterations(1000000);
TFitResultPtr r = new_fakesPhoEt->Fit("fitfunc_num","R S");
TF1 *fit = new_fakesPhoEt->GetFunction("fitfunc_num");
myfile << "num_coeff1 " << fit->GetParameter(0) << std::endl;
myfile << "num_coeff2 " << fit->GetParameter(1) << std::endl;
myfile << "num_lambd1 " << fit->GetParameter(2) << std::endl;
myfile << "num_lambd2 " << fit->GetParameter(3) << std::endl;

//std::ostringstream testname;
//TF1 *test_num[200][200];
//for(unsigned i(0); i < 200; i++){
//	for(unsigned j(0); j < 200; j++){
//		testname.str("");
//		testname << "test_" << i << "_" << j;
//		test_num[i][j] = new TF1(testname.str().c_str(), tmpjetfake_func, 35, 150, 4);
//		test_num[i][j]->SetParameters(iniCoeff_num, 20+i*5, iniLambda_num1, -0.001-0.001*j);
//		//TVirtualFitter::SetMaxIterations(1000000);
//	  int status = new_fakesPhoEt->Fit(testname.str().c_str(),"R");
//		TF1 *tmpf1 = new_fakesPhoEt->GetFunction(testname.str().c_str());
//		float diff = fabs(tmpf1->Eval(new_fakesPhoEt->GetBinCenter(55)) - new_fakesPhoEt->GetBinContent(55))/new_fakesPhoEt->GetBinContent(55);
//		if(tmpf1->GetParError(0)/tmpf1->GetParameter(0) < 1 && tmpf1->GetParError(1)/tmpf1->GetParameter(1) < 1 && tmpf1->GetParError(2)/tmpf1->GetParameter(2) < 1 && tmpf1->GetParError(3)/tmpf1->GetParameter(3) < 1 && diff < 0.2	)std::cout << "good point " << i*5 << " " << -1+0.005*j << " status = " << status << std::endl;
//		delete tmpf1;
//	}
//}
//

can_pad1->cd();          
TMatrixDSym cov = r->GetCovarianceMatrix(); 
r->Print("V");     
fitfunc_num->Draw("same");
float nominalvalue_num[115];
for(unsigned ibin(0); ibin < 115; ibin++){
	nominalvalue_num[ibin] =  fitfunc_num->Eval(35+ibin);
}
TVectorD mu(4) ;
mu(0) = r->Parameter(0); 
mu(1) = r->Parameter(1);
mu(2) = r->Parameter(2);
mu(3) = r->Parameter(3);
RooRealVar central_coeff1_num("central_coeff1_num","central_coeff1_num",mu(0)-r->ParError(0), mu(0)+r->ParError(0));
RooRealVar central_coeff2_num("central_coeff2_num","central_coeff2_num",mu(1)-r->ParError(1), mu(1)+r->ParError(1));
RooRealVar central_lambda1_num("central_lambda1_num","central_lambda1_num",mu(2)-r->ParError(2),mu(2)+r->ParError(2));
RooRealVar central_lambda2_num("central_lambda2_num","central_lambda2_num",mu(3)-r->ParError(3),mu(3)+r->ParError(3));

RooMultiVarGaussian mvg("mvg","mvg",RooArgList(central_coeff1_num,central_coeff2_num,central_lambda1_num,central_lambda2_num),mu,cov);
RooDataSet* toymcdata = mvg.generate(RooArgSet(central_coeff1_num,central_coeff2_num,central_lambda1_num,central_lambda2_num),NTOY);
std::ostringstream modelname;
TF1 *gen_num[NTOY];
TH1F *num_upper = new TH1F("num_upper","num_upper",115,35,150);
TH1F *num_lower = new TH1F("num_lower","num_lower",115,35,150);
float toyptvalue[115][NTOY];
float lowtoyptvalue[115];
float hightoyptvalue[115];
for(unsigned ii(0); ii < 115; ii++){
	lowtoyptvalue[ii] = 100000;
	hightoyptvalue[ii] = 0;
}

for(int i(0); i<NTOY; i++){
	double data1 = toymcdata->get(i)->getRealValue("central_coeff1_num");
	double data2 = toymcdata->get(i)->getRealValue("central_coeff2_num");
	double data3 = toymcdata->get(i)->getRealValue("central_lambda1_num");
	double data4 = toymcdata->get(i)->getRealValue("central_lambda2_num");
	if(data1 < 1e8 && data2 < 1e8 && data3 < 1e8 && data4 < 1e8){
		modelname.str("");
		modelname << "gen_num_" << i;
		gen_num[i] = new TF1(modelname.str().c_str(), tmpjetfake_func, 35, 150, 4);
		gen_num[i]->SetParameters(data1, data2, data3, data4);
		gen_num[i]->SetLineColorAlpha(kBlue, 0.35);
		//gen_num[i]->Draw("same");
		bool exception(false);
		std::cout << std::endl;
		for(unsigned ibin(1); ibin < new_fakesPhoEt->GetSize(); ibin++){
			std::cout << gen_num[i]->Eval(new_fakesPhoEt->GetBinCenter(ibin)) << " " << new_fakesPhoEt->GetBinContent(ibin) << std::endl; 
			if(gen_num[i]->Eval(new_fakesPhoEt->GetBinCenter(ibin)) - new_fakesPhoEt->GetBinContent(ibin) < -1*1.5*new_fakesPhoEt->GetBinError(ibin))exception=true;
		}
		if(exception)continue; 
		for(unsigned ibin(0); ibin < 115; ibin++){
			double estimated = gen_num[i]->Eval(35+ibin);
			toyptvalue[ibin][i] = estimated;
			if(lowtoyptvalue[ibin] > estimated)lowtoyptvalue[ibin] = estimated;
			if(hightoyptvalue[ibin]< estimated)hightoyptvalue[ibin] = estimated;
		}
	}
}

TCanvas *cangaus = new TCanvas("cangaus","",600,600);
cangaus->cd();
for(unsigned ibin(0); ibin < 115; ibin++){
	num_upper->SetBinContent(ibin, hightoyptvalue[ibin]); 
	num_lower->SetBinContent(ibin, lowtoyptvalue[ibin]);
}

can_pad1->cd();          
num_upper->Draw("L same");
num_lower->Draw("L same");
new_controlPhoEt->Draw("EP same");
new_fakesPhoEt->Draw("EP same");


for(unsigned ii(1); ii <= 115; ii++){
	myfile << "numerror " << ii-1 << " " <<  num_upper->GetBinContent(ii) - num_lower->GetBinContent(ii) << std::endl;
} 
for(unsigned ii(1); ii <= 115; ii++){
	myfile << "denerror " << ii-1 << " " <<  den_upper->GetBinContent(ii) - den_lower->GetBinContent(ii) << std::endl;
}

std::cout << std::endl; 
for(unsigned ibin(1); ibin < ratio->GetSize(); ibin++){
	myfile << "ratio " << ibin << " = " << ratio->GetBinContent(ibin) << "  est=" << fitfunc_num->Eval(ratio->GetBinCenter(ibin))/fitfunc_den->Eval(ratio->GetBinCenter(ibin)) << std::endl;
}
myfile.close();

	/************ jet tree **************************/ 
		TChain *jettree = new TChain("jetTree");
		//if(channel==1)jettree->Add("/uscms_data/d3/mengleis/FullStatusData/resTree_egsignal_DoubleEG_ReMiniAOD.root");
		if(channel==1)jettree->Add("/uscms_data/d3/mengleis/resTree_egsignal_DoubleEG_ReMiniAOD_test.root");
		else if(channel==2)jettree->Add("/uscms_data/d3/mengleis/test/resTree_mgsignal_MuonEG_FebReminiAOD_MiniIso_JES.root");

		float jetphoEt(0);
		float jetphoEta(0);
		float jetphoPhi(0);
		float jetlepPt(0);
		float jetlepEta(0);
		float jetlepPhi(0);
		float jetsigMT(0);
		float jetsigMET(0);
		float jetsigMETPhi(0);
		float jetdPhiLepMET(0);
		int   jetnVertex(0);
		float jetdRPhoLep(0);
		float jetHT(0);
		float jetnJet(0);
		
		jettree->SetBranchAddress("phoEt",     &jetphoEt);
		jettree->SetBranchAddress("phoEta",    &jetphoEta);
		jettree->SetBranchAddress("phoPhi",    &jetphoPhi);
		jettree->SetBranchAddress("lepPt",     &jetlepPt);
		jettree->SetBranchAddress("lepEta",    &jetlepEta);
		jettree->SetBranchAddress("lepPhi",    &jetlepPhi);
		jettree->SetBranchAddress("sigMT",     &jetsigMT);
		jettree->SetBranchAddress("sigMET",    &jetsigMET);
		jettree->SetBranchAddress("sigMETPhi", &jetsigMETPhi);
		jettree->SetBranchAddress("dPhiLepMET",&jetdPhiLepMET);
		jettree->SetBranchAddress("nVertex",   &jetnVertex);
		jettree->SetBranchAddress("dRPhoLep",  &jetdRPhoLep);
		jettree->SetBranchAddress("HT",        &jetHT);
		jettree->SetBranchAddress("nJet",      &jetnJet);
	 
		for (unsigned ievt(0); ievt<jettree->GetEntries(); ++ievt){//loop on entries
			jettree->GetEntry(ievt);

			if(jetphoEt <= 35)continue;
			double w_jet(0);
			w_jet = fitfunc_num->Eval(jetphoEt)/fitfunc_den->Eval(jetphoEt);
			if((channel==1) && jetphoEt < 37)w_jet = 0.20;
			else if(channel==2 && jetphoEt < 37)w_jet = 0.238459;
		
			p_jetbkgPhoEt->Fill(jetphoEt,w_jet);
			p_jetbkgPhoEta->Fill(jetphoEta, w_jet);
			p_jetbkgMET->Fill(jetsigMET,w_jet);
			p_jetbkgMt->Fill(jetsigMT, w_jet);
			p_jetbkgLepPt->Fill(jetlepPt, w_jet);
			p_jetbkgLepEta->Fill(jetlepEta, w_jet);
			p_jetbkgdRPhoLep->Fill(jetdRPhoLep, w_jet);
		} 
		
		p_jetbkgPhoEt->Sumw2();
		p_jetbkgPhoEta->Sumw2();
		p_jetbkgMET->Sumw2();
		p_jetbkgMt->Sumw2();

	c_pt->SaveAs("JetFakeRate_transfer_mg_July22.pdf");
	outputfile->Write();
}


