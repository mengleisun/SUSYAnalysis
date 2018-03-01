#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<ctime>
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TF1.h"
#include "TF3.h"
#include "TH2F.h"
#include "TMath.h"
#include "TSystem.h"
#include "TFractionFitter.h"
#include "TLatex.h"
#include "TStyle.h"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooAbsReal.h"
#include "RooConstVar.h"
#include "RooDataHist.h"
#include "RooKeysPdf.h"
#include "RooPlot.h"
#include "RooHistPdf.h"
#include "RooNumIntConfig.h"
#include "RooFFTConvPdf.h"
#include "RooPolynomial.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooFitResult.h"

#include "../../../include/tdrstyle.C"
#include "../../../include/analysis_rawData.h"
#include "../../../include/analysis_photon.h"
#include "../../../include/analysis_muon.h"
#include "../../../include/analysis_ele.h"
#include "../../../include/analysis_mcData.h"
#include "../../../include/analysis_tools.h"
#include "../../../include/analysis_fakes.h"

int  eventType = 1; // 1 = eg, 2 = mg, 3 = MCeg, 4 = MCmg
bool useMC = true;
bool doIterate = true;
float METLOWCUT = 0;
float METHIGHCUT = 70;

int
FitJetFake(float lowercut, float uppercut, int detType, float loweta, float higheta){

	setTDRStyle();   
	time_t now = time(0);
	int nBinTotal = 20;
	if(lowercut > 150)nBinTotal = 10;

	char lowername[3];
	sprintf(lowername, "%d", (int)lowercut);
	char uppername[3];
	if(uppercut < 1000){
	  sprintf(uppername, "%d", (int)uppercut);
	}
	else sprintf(uppername, "Inf");

	std::cout << "start fitting " << std::endl;
	gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");
	ofstream myfile;
	std::string  outputType;
	std::string  outputDet;
	std::ostringstream  myfilename;  myfilename.str("");
	std::ostringstream  Hist1Dname;  Hist1Dname.str("");
	std::ostringstream  Hist2Dname;  Hist2Dname.str("");
	switch(eventType){
		case 1: outputType = "DoubleEG-"; break;
		case 2: outputType = "MuonEG-"; break;
		case 3: outputType = "MCEG-"; break;
		case 4: outputType = "MCMG-"; break;
	}
	if(detType == 1)outputDet = "EB";
	else if(detType == 2)outputDet = "EE";

	//myfilename << "/uscms_data/d3/mengleis/SUSYAnalysis/test/jetFakePho/result/JetFakeRate-" << outputType << outputDet << ".txt";
	myfilename << "./JetFakeRate-" << outputType << outputDet << ".txt";
	//Hist1Dname << "/uscms_data/d3/mengleis/SUSYAnalysis/test/jetFakePho/plot/frac-" << lowername << "-" << uppername << "-" << outputType << outputDet << ".pdf";
	//Hist2Dname << "/uscms_data/d3/mengleis/SUSYAnalysis/test/jetFakePho/plot/can2D-" << lowername << "-" << uppername << "-" << outputType << outputDet << ".pdf";
	Hist1Dname << "./frac-" << lowername << "-" << uppername << "-" << outputType << outputDet << ".pdf";
	Hist2Dname << "./can2D-" << lowername << "-" << uppername << "-" << outputType << outputDet << ".pdf";
	myfile.open(myfilename.str().c_str(), std::ios_base::app | std::ios_base::out);

	ofstream logfile;
	logfile.open("fittingquality.log", std::ios_base::app | std::ios_base::out);
	logfile << outputType << " " << lowername << " " << uppername << std::endl;
	

	TChain *datatree = new TChain("hadronTree");
	switch(eventType){
		case 1: datatree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_egsignal_DoubleEG_ReMiniAOD_FullEcal.root"); break;
		case 2: datatree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_mgsignal_MuonEG_FullEcal.root"); break;
		case 3: datatree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_egsignal_DY.root");
						datatree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_egsignal_WJet.root"); 
						datatree->Add("/uscms_data/d3/mengleis/FullStatusOct/resTree_egsignal_QCD40.root"); 
						break;
	}	

	std::ostringstream mcTreename; 
	mcTreename.str("");
	if(eventType == 1 || eventType == 3)mcTreename << "egTree";
	else if(eventType == 2 || eventType == 4)mcTreename << "mgTree";
	TChain *mctree = new TChain(mcTreename.str().c_str());
	mctree->Add("/uscms_data/d3/mengleis/FullStatusOct/plot_hadron_GJet.root");
	
	/*************************************/
	/*      double normfactor = par[0];  */  
	/*      double slope = par[1];       */
	/*      double constant = par[2];    */
	/*      double index = par[3];       */
	/*      double coeff = par[4];       */
	/*      double vtx_constant = par[5];*/
	/*      double vtx_slope = par[6];   */
	/*************************************/
	std::ostringstream elefake_config;
	elefake_config.str("");
	elefake_config << "/uscms_data/d3/mengleis/SUSYAnalysis/test/eleFakePho/DataResult/EleFakeRate-ByPtVtx-" << outputDet << ".txt";
	std::ifstream elefake_file(elefake_config.str().c_str());
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
	TF3 f_elefake("f_elefake", fakerate_func,10,1000,0,100,0,2.5,7);
	f_elefake.SetParameters(scalefactor, ptslope, ptconstant, ptindex, ptcoeff, vtxconst, vtxslope);

	double   StandardCut = 0;
	double   StandardIso = 0;
	Double_t SigmaCutLower_EB[]={0.0103,0.0104,0.0105,0.0106,0.0107,0.0108,0.0109,0.0110,0.0111,0.0112};
	Double_t SigmaCutUpper_EB[]={0.0140,0.0145,0.0150,0.0155,0.0160,0.0165,0.0170,0.0175,0.0180,0.0185};
	Double_t SigmaCutLower_EE[]={0.03013,0.0302,0.0303,0.0304,0.0305,0.0306,0.0307,0.0308,0.0309,0.031};
	Double_t SigmaCutUpper_EE[]={0.035,0.036,0.037,0.038,0.039,0.04};
	std::vector<double> SigmaCutLower;
	std::vector<double> SigmaCutUpper;
	SigmaCutLower.clear();
	SigmaCutUpper.clear();
	if(detType == 1){
		StandardCut = 0.0103;
		StandardIso = 1.295;
		for(unsigned i(0); i<sizeof(SigmaCutLower_EB)/sizeof(Double_t); i++)SigmaCutLower.push_back(SigmaCutLower_EB[i]);
		for(unsigned i(0); i<sizeof(SigmaCutUpper_EB)/sizeof(Double_t); i++)SigmaCutUpper.push_back(SigmaCutUpper_EB[i]);
	}
	else if(detType == 2){
		StandardCut = 0.03013;
		StandardIso = 1.011;
		for(unsigned i(0); i<sizeof(SigmaCutLower_EE)/sizeof(Double_t); i++)SigmaCutLower.push_back(SigmaCutLower_EE[i]);
		for(unsigned i(0); i<sizeof(SigmaCutUpper_EE)/sizeof(Double_t); i++)SigmaCutUpper.push_back(SigmaCutUpper_EE[i]);
	}
	unsigned nUpper = SigmaCutUpper.size(); 
	unsigned nLower = SigmaCutLower.size();
	std::ostringstream hname;
	hname.str(""); hname <<"ChIsoTar-pt-" << lowername << "-" << uppername;
	TH1D* h_target= new TH1D(hname.str().c_str(),";Iso_{h^{#pm}} (GeV);",nBinTotal, 0.0, 20); 
	hname.str(""); hname <<"mcChIso-pt-" <<  lowername << "-" << uppername;
	TH1D* mc_sig = new TH1D(hname.str().c_str(), ";Iso_{h^#pm} (GeV);",nBinTotal, 0.0, 20);
	hname.str(""); hname <<"mcChIso-fake-pt-" << lowername << "-" << uppername;
	TH1D* mc_fake = new TH1D(hname.str().c_str(),";Iso_{h^#pm} (GeV);",nBinTotal, 0.0, 20);
	TH1D* h_bg[nLower][nUpper];
	
	TH1*  mc_predict[nLower][nUpper];
	TH1D* mc_sbcontamination[nLower][nUpper];
	float templateCorrFactor[nLower][nUpper];
	TObjArray *templateContainer[nLower][nUpper];
	TFractionFitter* fitter[nLower][nUpper];
	TH1D* result[nLower][nUpper];
	TCanvas *can[nLower][nUpper]; 
	hname.str(""); hname <<"can2D" << lowername << "-" << uppername;
	TCanvas *can2D = new TCanvas(hname.str().c_str(), hname.str().c_str(),600,600);
	for(unsigned iUpper(0); iUpper<nUpper; iUpper++){
		for(unsigned iLower(0); iLower<nLower; iLower++){
			templateContainer[iLower][iUpper] = new TObjArray(2);
			hname.str("");
			hname << "can" << lowername << "-" << uppername << "_" << SigmaCutLower[iLower] << "_" << SigmaCutUpper[iUpper] ;
			can[iLower][iUpper]=new TCanvas(hname.str().c_str(),hname.str().c_str(),600,600);
			can[iLower][iUpper]->SetBottomMargin(0.15);
			hname.str(""); hname <<"ChIsoHad-pt-" <<  lowername << "-" << uppername << "_" << SigmaCutLower[iLower] << "_" << SigmaCutUpper[iUpper];
			h_bg[iLower][iUpper] = new TH1D(hname.str().c_str(), hname.str().c_str(), nBinTotal, 0.0, 20);
			hname.str(""); hname <<"mcChIso-sideband-pt-" <<  lowername << "-" << uppername << "_" << SigmaCutLower[iLower] << "_" << SigmaCutUpper[iUpper];
			mc_sbcontamination[iLower][iUpper] = new TH1D(hname.str().c_str(), hname.str().c_str(), nBinTotal, 0.0, 20);   
		}
	}

	//TFile *outputfile = TFile::Open("/uscms_data/d3/mengleis/SUSYAnalysis/test/jetFakePho/result/JetFakeRate-MCeg-FullEcal.root","RECREATE");
	TFile *outputfile = TFile::Open("./JetFakeRate-MCeg-FullEcal.root","RECREATE");
	outputfile->cd();
	hname.str("");
	hname << "fracHad2D_" << lowername << "-" << uppername;
	TH2F* fracHad2D = new TH2F(hname.str().c_str(),"hadron fraction;Lower #sigma_{i#etai#eta} Threshold(GeV);Upper #sigma_{i#etai#eta} Threshold(GeV)",10,SigmaCutLower[0], SigmaCutLower[SigmaCutLower.size()-1],26,0.014,0.04);
	hname.str("");
	hname << "fracHad1D_" << lowername << "-" << uppername;
	TH1D* fracHad1D = new TH1D(hname.str().c_str(),hname.str().c_str(),500,0,1.0);


//************ Signal Tree **********************//
	float mc_phoEt(0);
	float mc_phoEta(0); 
	float mc_phoPhi(0); 
	float mc_phoSigma(0);
	float mc_phoChIso(0);
	std::vector<int> *mc_mcPID=0;
	std::vector<float> *mc_mcEta=0;
	std::vector<float> *mc_mcPhi=0;
	std::vector<float> *mc_mcPt=0;
	std::vector<int> *mc_mcMomPID=0; 

	mctree->SetBranchAddress("phoEt",     &mc_phoEt);
	mctree->SetBranchAddress("phoEta",    &mc_phoEta);
	mctree->SetBranchAddress("phoPhi",    &mc_phoPhi);
	mctree->SetBranchAddress("phoSigma",  &mc_phoSigma);
	mctree->SetBranchAddress("phoChIso",  &mc_phoChIso);
	if(useMC){
		mctree->SetBranchAddress("mcPID",   &mc_mcPID);
		mctree->SetBranchAddress("mcEta",   &mc_mcEta);
		mctree->SetBranchAddress("mcPhi",   &mc_mcPhi);
		mctree->SetBranchAddress("mcPt",    &mc_mcPt);
		mctree->SetBranchAddress("mcMomPID",&mc_mcMomPID);
	}

	float crosssection(1);
	float ntotalevent(1);
	float phoEt(0);
	float phoEta(0);
	float sigMET(0);
	float phoSigma(0);
	float phoChIso(0);
	std::vector<float> *eleproxyEt=0;
	std::vector<float> *eleproxyEta=0;
	std::vector<float> *eleproxyPhi=0;
	std::vector<float> *eleproxySigma=0;
	std::vector<float> *eleproxyChIso=0;
	std::vector<int>   *eleproxynVertex=0;
	if(eventType == 3 || eventType == 4){
		datatree->SetBranchAddress("crosssection",&crosssection);
		datatree->SetBranchAddress("ntotalevent", &ntotalevent);
	}
	datatree->SetBranchAddress("phoEt",     &phoEt);
	datatree->SetBranchAddress("phoEta",    &phoEta);
	datatree->SetBranchAddress("sigMET",    &sigMET);
	datatree->SetBranchAddress("phoSigma",  &phoSigma);
	datatree->SetBranchAddress("phoChIso",  &phoChIso);
	datatree->SetBranchAddress("eleproxyEt",&eleproxyEt);
	datatree->SetBranchAddress("eleproxyEta",&eleproxyEta);
	datatree->SetBranchAddress("eleproxyPhi",&eleproxyPhi);
	datatree->SetBranchAddress("eleproxySigma", &eleproxySigma);
	datatree->SetBranchAddress("eleproxyChIso", &eleproxyChIso);
	datatree->SetBranchAddress("eleproxynVertex",&eleproxynVertex);


	double hadfrac(0);
	double fittingError(0);
	double systematicError(0), lowestfrac(1), highestfrac(0);
	double mcfakerate(0);
	double nden(0), totalden(0), nsignalregion(0);
	double nMCtarget(0),  nMCfaketarget(0);
	double nnum[nLower][nUpper], totalnum[nLower][nUpper];
	double nsideband[nLower][nUpper], nsbPassIso[nLower][nUpper];
	for(unsigned ii(0); ii < nLower; ii++){
		for(unsigned jj(0); jj < nUpper; jj++){
			nnum[ii][jj] = 0;
			totalnum[ii][jj] = 0;
			nsideband[ii][jj] = 0;
			nsbPassIso[ii][jj] = 0;
		}
	}

	for(unsigned ievt(0); ievt < mctree->GetEntries(); ievt++){
		mctree->GetEntry(ievt);

		if(mc_phoEt < lowercut || mc_phoEt > uppercut)continue;
		if(detType == 1 && fabs(mc_phoEta) > 1.4442)continue;
		else if(detType == 2 && (fabs(mc_phoEta) < loweta || fabs(mc_phoEta) > higheta) )continue;

		bool isFake(true);
		unsigned mcIndex(0);
		unsigned mcPhoIndex(0);
		float mindR(0.7);
		float phodR(0.7);
		bool hasMatch(false);
		for(unsigned ii(0); ii < mc_mcPID->size(); ii++){
			float dR = DeltaR(mc_phoEta, mc_phoPhi, (*mc_mcEta)[ii], (*mc_mcPhi)[ii]);
			float dE = fabs(mc_phoEt - (*mc_mcPt)[ii])/mc_phoEt;
			if((*mc_mcPID)[ii] == 22){
				phodR = DeltaR(mc_phoEta, mc_phoPhi, (*mc_mcEta)[ii], (*mc_mcPhi)[ii]);
				mcPhoIndex = ii;
			}
			if(dR < mindR && dE < 0.7){mindR = dR; mcIndex = ii; hasMatch = true;}
		}
		if(phodR < 0.1)mcIndex = mcPhoIndex;
		if(hasMatch)isFake = isHad(fabs((*mc_mcPID)[mcIndex]), fabs((*mc_mcMomPID)[mcIndex]));

		bool isTrueTemplate = (mc_phoSigma <= StandardCut && !isFake);
		bool isInTargetRegion = (mc_phoSigma <= StandardCut && mc_phoChIso < StandardIso);

		if(mc_phoSigma <= StandardCut && isFake)mc_fake->Fill(mc_phoChIso);
		if(isTrueTemplate){
			mc_sig->Fill(mc_phoChIso);
			nsignalregion += 1;
		}
		if(isInTargetRegion){
			nMCtarget += 1;
			if(isFake)nMCfaketarget += 1; 
		}
		for(unsigned iUpper(0); iUpper<nUpper; iUpper++){
			for(unsigned iLower(0); iLower<nLower; iLower++){
				if(mc_phoSigma > SigmaCutLower[iLower] && mc_phoSigma < SigmaCutUpper[iUpper] && !isFake){
					mc_sbcontamination[iLower][iUpper]->Fill(mc_phoChIso);
					nsideband[iLower][iUpper] += 1;
					if(mc_phoChIso < StandardIso)nsbPassIso[iLower][iUpper] += 1;
				}
			}
		}					
	}
	std::cout << "MC tree has been filled" << std::endl;

	if(nMCtarget > 0)mcfakerate = 1.0*nMCfaketarget/nMCtarget;

	for(unsigned ievt(0); ievt < datatree->GetEntries(); ievt++){
		datatree->GetEntry(ievt);
		if(sigMET < METLOWCUT || sigMET > METHIGHCUT)continue;

		double LumiWeight = 1;
		if(eventType == 3 || eventType == 4)LumiWeight = 35.8*crosssection*1000/ntotalevent;

		for(unsigned iele(0); iele < eleproxyEt->size(); iele++){
			if(detType == 1 && fabs( (*eleproxyEta)[iele]) > 1.4442)continue;
			else if(detType == 2 && ( fabs( (*eleproxyEta)[iele]) < loweta || fabs( (*eleproxyEta)[iele]) > higheta) )continue;
			double w_ele = -1.0*f_elefake( (*eleproxyEt)[iele], (*eleproxynVertex)[iele], fabs( (*eleproxyEta)[iele]));
			w_ele = w_ele*LumiWeight;
			if( (*eleproxyEt)[iele] < lowercut || (*eleproxyEt)[iele]  > uppercut)continue;

			if( (*eleproxySigma)[iele] <= StandardCut){
				h_target->Fill( (*eleproxyChIso)[iele], w_ele);
				totalden+= w_ele;
				if( (*eleproxyChIso)[iele] < StandardIso)nden+= w_ele;
			}
			else{
				for(unsigned iUpper(0); iUpper<nUpper; iUpper++){
					for(unsigned iLower(0); iLower<nLower; iLower++){
						if( (*eleproxySigma)[iele] > SigmaCutLower[iLower] && (*eleproxySigma)[iele] < SigmaCutUpper[iUpper]){
							h_bg[iLower][iUpper]->Fill( (*eleproxyChIso)[iele], w_ele);
							totalnum[iLower][iUpper] += w_ele;
							if( (*eleproxyChIso)[iele] < StandardIso)nnum[iLower][iUpper] += w_ele;
						}
					}
				}
			}

		}

		if(phoEt < lowercut || phoEt > uppercut )continue;
		if(sigMET < METLOWCUT || sigMET > METHIGHCUT)continue;
		if(detType == 1 && fabs(phoEta) > 1.4442)continue;
		else if(detType == 2 && ( fabs(phoEta) < loweta || fabs(phoEta) > higheta) )continue;

		if(phoSigma <= StandardCut){
			h_target->Fill(phoChIso,LumiWeight);
			totalden+=1*LumiWeight;
			if(phoChIso < StandardIso)nden+=1*LumiWeight;
		}
		else{
			for(unsigned iUpper(0); iUpper<nUpper; iUpper++){
				for(unsigned iLower(0); iLower<nLower; iLower++){
					if(phoSigma > SigmaCutLower[iLower] && phoSigma < SigmaCutUpper[iUpper]){
						h_bg[iLower][iUpper]->Fill(phoChIso, LumiWeight);
						totalnum[iLower][iUpper] +=1*LumiWeight;
						if(phoChIso < StandardIso)nnum[iLower][iUpper] +=1*LumiWeight;
					}
				}
			}
		}

	}
	h_target->Sumw2();
	for(unsigned iUpper(0); iUpper<nUpper; iUpper++)
		for(unsigned iLower(0); iLower<nLower; iLower++)
			h_bg[iLower][iUpper]->Sumw2();

	std::cout << "data tree has been filled" << std::endl;

	TCanvas *districan = new TCanvas("districan","",1200,600);
	districan->Divide(2);
	districan->cd(1);
	h_target->Draw();
	mc_sig->Scale(h_target->Integral(1,nBinTotal)/mc_sig->Integral(1,nBinTotal));
	mc_sig->SetLineColor(kRed);
	mc_sig->Draw("hist same");	
	districan->cd(2);
	h_bg[0][0]->Draw();

	std::cout << "target " << h_target->Integral(1,nBinTotal) << std::endl;


	for(unsigned iUpper(0); iUpper<nUpper; iUpper++){
		for(unsigned iLower(0); iLower<nLower; iLower++){
			double iteratorUncertainty(1), lastUpdatePurity(1);
			double frac_passSigma = nsbPassIso[iLower][iUpper]/(1.0*mc_sbcontamination[iLower][iUpper]->GetEntries());
			if(mc_sbcontamination[iLower][iUpper]->GetEntries()==0)frac_passSigma = 0;
			double phoPurity(0);
			templateCorrFactor[iLower][iUpper] = 0.0;
	int itertime = 0;
			if(doIterate){
				while(iteratorUncertainty > 0.01){

					itertime += 1; logfile << "iter " << itertime << std::endl;
					TH1D *hist_new=(TH1D*)h_bg[iLower][iUpper]->Clone();
					hist_new->Add(mc_sbcontamination[iLower][iUpper], -1*templateCorrFactor[iLower][iUpper]/mc_sbcontamination[iLower][iUpper]->GetEntries());
					hist_new->Sumw2();
					templateContainer[iLower][iUpper]->Add(mc_sig);
					templateContainer[iLower][iUpper]->Add(hist_new);
					fitter[iLower][iUpper]= new TFractionFitter(h_target, templateContainer[iLower][iUpper], "V");
					fitter[iLower][iUpper]->Constrain(1,0.0,1.0);
					Int_t status = fitter[iLower][iUpper]->Fit();
					double tmpSig, tmpErrorSig;
					if(status == 0){
						fitter[iLower][iUpper]->GetResult(0, tmpSig,  tmpErrorSig);
						phoPurity = tmpSig;
					}
					iteratorUncertainty = fabs(lastUpdatePurity - phoPurity)/phoPurity;
					lastUpdatePurity = phoPurity;
					templateCorrFactor[iLower][iUpper] = totalden*phoPurity*nsideband[iLower][iUpper]/nsignalregion;
					templateContainer[iLower][iUpper]->Clear();
					if( fabs(SigmaCutLower[iLower]- StandardCut) < 0.00005 && iUpper == 0)logfile << "  Standard "  << "iteratorUncertainty " << iteratorUncertainty << " phoPurity " << phoPurity << "  templateCorrFactor " << templateCorrFactor[iLower][iUpper] << std::endl;
					else logfile << "  " << SigmaCutLower[iLower] << " " << SigmaCutLower[iUpper] << " " << "iteratorUncertainty " << iteratorUncertainty << " phoPurity " << phoPurity << "  templateCorrFactor " << templateCorrFactor[iLower][iUpper] << std::endl;
					hist_new->Delete();
				}
			}
			else{
				iteratorUncertainty = 0;
				frac_passSigma = 0;
				templateCorrFactor[iLower][iUpper] = 0;
				templateContainer[iLower][iUpper]->Add(mc_sig);
				templateContainer[iLower][iUpper]->Add(h_bg[iLower][iUpper]);
				fitter[iLower][iUpper]= new TFractionFitter(h_target, templateContainer[iLower][iUpper], "V");
				fitter[iLower][iUpper]->Constrain(1,0.0,1.0);
				Int_t status = fitter[iLower][iUpper]->Fit();
			}
			if(iteratorUncertainty <= 0.01 || !doIterate){
				double den = nden/totalden;  
				double denerror = (1.0/sqrt(nden)+ 1.0/sqrt(totalden))*den;
				double num = (nnum[iLower][iUpper] - templateCorrFactor[iLower][iUpper]*frac_passSigma)/(totalnum[iLower][iUpper] - templateCorrFactor[iLower][iUpper]); 
				double numerror = (1.0/sqrt(nnum[iLower][iUpper]) + 1.0/sqrt(totalnum[iLower][iUpper]))*num;
				double fracBkg, errorBkg;
				double fracSig, errorSig;
				result[iLower][iUpper] = (TH1D*)fitter[iLower][iUpper]->GetPlot();
				result[iLower][iUpper]->SetMinimum(10);
				fitter[iLower][iUpper]->GetResult(0, fracSig, errorSig);
				fitter[iLower][iUpper]->GetResult(1, fracBkg, errorBkg);
				can[iLower][iUpper]->cd();
				gStyle->SetOptStat(0);
				hname.str("");
				hname << lowername << "< pt <" << uppername; 
				h_target->SetTitle(hname.str().c_str());
				h_target->SetMinimum(h_target->GetBinContent(nBinTotal));
				result[iLower][iUpper]->SetMinimum(100);
				gPad->SetLogy();
				h_target->GetXaxis()->SetTitleOffset(0.8);
				h_target->SetMarkerStyle(20);
				h_target->Draw("Ep");
				result[iLower][iUpper]->SetLineColor(kBlue);
				result[iLower][iUpper]->SetFillStyle(1001);
				result[iLower][iUpper]->SetFillColor(kBlue);
				result[iLower][iUpper]->SetMarkerColor(kBlue);
				result[iLower][iUpper]->Draw("same");

				mc_predict[iLower][iUpper] = fitter[iLower][iUpper]->GetMCPrediction(1);
				mc_predict[iLower][iUpper]->SetMinimum(100);
				mc_predict[iLower][iUpper]->Scale(1.0/mc_predict[iLower][iUpper]->Integral(1,21)*result[iLower][iUpper]->Integral(1,21)*fracBkg);
				mc_predict[iLower][iUpper]->SetLineColor(kGreen);
				mc_predict[iLower][iUpper]->SetFillStyle(1001);
				mc_predict[iLower][iUpper]->SetFillColor(kGreen);
				mc_predict[iLower][iUpper]->SetMarkerColor(kGreen);
				mc_predict[iLower][iUpper]->Draw("hist same");

		////		mc_fake->SetLineColor(kRed);
		////		mc_fake->Draw("same");

				h_target->Draw("Ep same");
				TLatex* latex = new TLatex();
		//		hname.str("");
		//		hname << "Bkg% = " << fracBkg << " #pm " << errorBkg << std::endl;
		//		latex->DrawLatex(2,0.9*h_target->GetMaximum(),hname.str().c_str()); 
				hname.str("");
				hname << "#chi^{2}/ndof = " << fitter[iLower][iUpper]->GetChisquare() << "/" << fitter[iLower][iUpper]->GetNDF(); 
				latex->DrawLatex(5,0.15*h_target->GetMaximum(),hname.str().c_str());
		//		hname.str("");
		//		hname << "Hadfrac="<< num*fracBkg/den << " #pm " << (numerror/num + denerror/den + errorBkg/fracBkg)*num*fracBkg/den;
		//		latex->DrawLatex(2,0.15*h_target->GetMaximum(),hname.str().c_str());
		 
				if(SigmaCutLower[iLower] == StandardCut && iUpper == 0){
					hadfrac =  num*fracBkg/den;
					fittingError = (numerror/num + denerror/den + errorBkg/fracBkg)*num*fracBkg/den ;
					TLegend *leg = new TLegend(0.6,0.75,0.87,0.9);
					leg->AddEntry(h_target, "observed");
					leg->AddEntry(result[iLower][iUpper], "true photons");
					leg->AddEntry(mc_predict[iLower][iUpper], "hadrons");
					leg->Draw("same");
					can[iLower][iUpper]->SaveAs(Hist1Dname.str().c_str());
				}
				if(fracBkg >0 && fracBkg < 1)fracHad2D->Fill(SigmaCutLower[iLower]+0.00005, SigmaCutUpper[iUpper]+0.00005, num*fracBkg/den);
				if(fracBkg >0 && fracBkg < 1)fracHad1D->Fill(num*fracBkg/den);
				if(lowestfrac > num*fracBkg/den)lowestfrac = num*fracBkg/den;
				if(highestfrac < num*fracBkg/den)highestfrac = num*fracBkg/den; 
			}
		}
	}

	if(hadfrac < 0.001){
		hadfrac = fracHad1D->GetMean();	
	}

	can2D->cd();
	fracHad2D->Draw("colz");
	can2D->SaveAs(Hist2Dname.str().c_str());
	systematicError = highestfrac - lowestfrac;
	myfile << " " << lowername << " " << uppername << " " << hadfrac << " "; 
	myfile << sqrt(fittingError*fittingError + systematicError*systematicError) << " ";
	myfile << systematicError << "	" << " ";
	myfile << mcfakerate;
	myfile << outputDet <<  std::endl;
	char* dt = ctime(&now);
	if(uppercut > 500){
		myfile<< std::endl;
		myfile<< "ibin pt_lower pt_upper fakerate error errorsystematic" << std::endl;
		myfile << dt << std::endl;
		myfile << "SigmaCutLower = ";
		for(unsigned l(0); l < nLower; l++)myfile << SigmaCutLower[l] << " ";
		myfile << std::endl;
		for(unsigned u(0); u < nUpper; u++)myfile << SigmaCutUpper[u] << " ";
		myfile << std::endl;
	} 
	myfile.close();
	logfile.close();
	outputfile->Write();
	return 1;
}
