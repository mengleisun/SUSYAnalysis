#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<ctime>
#include "TH1F.h"
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

#include "../../include/tdrstyle.C"
#include "../../include/analysis_rawData.h"
#include "../../include/analysis_photon.h"
#include "../../include/analysis_muon.h"
#include "../../include/analysis_ele.h"
#include "../../include/analysis_mcData.h"
#include "../../include/analysis_tools.h"
#include "../../include/analysis_fakes.h"

bool useMC = true;
bool doIterate = true;
float METLOWCUT = 0;
float METHIGHCUT = 1000;

int
FitJetFake(float lowercut, float uppercut){

	setTDRStyle();   
  time_t now = time(0);

	gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");
	ofstream myfile;
	myfile.open("JetFakeRate-ISR.txt", std::ios_base::app | std::ios_base::out);

  char lowername[3];
  sprintf(lowername, "%d", (int)lowercut);
  char uppername[3];
  if(uppercut < 1000){
    sprintf(uppername, "%d", (int)uppercut);
  }
  else sprintf(uppername, "Inf");

	double SigmaCutLower[]={0.0103,0.0104,0.0105,0.0106,0.0107,0.0108,0.0109,0.011};
	double SigmaCutUpper[]={0.0140,0.0145,0.0130,0.0135};
	unsigned nUpper = sizeof(SigmaCutUpper)/sizeof(double);
	unsigned nLower = sizeof(SigmaCutLower)/sizeof(double);
	std::ostringstream hname;
  hname.str(""); hname <<"ChIsoTar-pt-" << lowername << "-" << uppername;
	TH1F* h_target= new TH1F(hname.str().c_str(),";Iso_{h^{#pm}} (GeV);",10, 0.0, 20); 
  hname.str(""); hname <<"mcChIso-pt-" <<  lowername << "-" << uppername;
	TH1F* mc_sig = new TH1F(hname.str().c_str(), ";Iso_{h^#pm} (GeV);",10, 0.0, 20);
  hname.str(""); hname <<"mcChIso-fake-pt-" << lowername << "-" << uppername;
	TH1F* mc_fake = new TH1F(hname.str().c_str(),";Iso_{h^#pm} (GeV);",10, 0.0, 20);
	TH1F* h_bg[nLower][nUpper];
	
	TH1*  mc_predict[nLower][nUpper];
	TH1F* mc_sbcontamination[nLower][nUpper];
  float templateCorrFactor[nLower][nUpper];
	TObjArray *templateContainer[nLower][nUpper];
	TFractionFitter* fitter[nLower][nUpper];
	TH1F* result[nLower][nUpper];
	TCanvas *can[nLower][nUpper]; 
  hname.str(""); hname <<"can2D" << lowername << "-" << uppername;
	TCanvas *can2D = new TCanvas(hname.str().c_str(), hname.str().c_str(),600,600);
	for(unsigned iUpper(0); iUpper<nUpper; iUpper++){
		for(unsigned iLower(0); iLower<nLower; iLower++){
			templateContainer[iLower][iUpper] = new TObjArray(2);
			hname.str("");
			hname << "can" << lowername << "-" << uppername << "_" << SigmaCutLower[iLower] << "_" << SigmaCutUpper[iUpper] ;
			can[iLower][iUpper]=new TCanvas(hname.str().c_str(),hname.str().c_str(),600,600);
			hname.str(""); hname <<"ChIsoHad-pt-" <<  lowername << "-" << uppername << "_" << SigmaCutLower[iLower] << "_" << SigmaCutUpper[iUpper];
			h_bg[iLower][iUpper] = new TH1F(hname.str().c_str(), hname.str().c_str(), 10, 0.0, 20);
			hname.str(""); hname <<"mcChIso-sideband-pt-" <<  lowername << "-" << uppername << "_" << SigmaCutLower[iLower] << "_" << SigmaCutUpper[iUpper];
			mc_sbcontamination[iLower][iUpper] = new TH1F(hname.str().c_str(), hname.str().c_str(),10, 0.0, 20);   
		}
	}
	TCanvas *cantemplate[nLower][nUpper]; 
  hname.str(""); hname <<"template" << lowername << "-" << uppername;
	for(unsigned iUpper(0); iUpper<nUpper; iUpper++){
		for(unsigned iLower(0); iLower<nLower; iLower++){
			hname << "_" << iUpper << "_" << iLower;
			cantemplate[iLower][iUpper]=new TCanvas(hname.str().c_str(),hname.str().c_str(),1200,600);
			cantemplate[iLower][iUpper]->Divide(2);
		}
	}

	TFile *outputfile = TFile::Open("JetFakeRate-ISR.root","RECREATE");
	outputfile->cd();
	hname.str("");
	hname << "fracHad2D_" << lowername << "-" << uppername;
	TH2F* fracHad2D = new TH2F(hname.str().c_str(),"hadron fraction;Lower #sigma_{i#etai#eta} Threshold(GeV);Upper #sigma_{i#etai#eta} Threshold(GeV)",10,0.0103,0.0113,10,0.014,0.019);
	hname.str("");
	hname << "fracHad1D_" << lowername << "-" << uppername;
	TH1F* fracHad1D = new TH1F(hname.str().c_str(),hname.str().c_str(),500,0,1.0);


//************ Signal Tree **********************//
	TChain *mctree = new TChain("egTree");
	mctree->Add("/uscms_data/d3/mengleis/usefuldata/plot_hadron_mgGJet.root");
	float mc_phoEt(0);
	float mc_phoEta(0); 
	float mc_phoPhi(0); 
	float mc_sigMET(0);
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
	mctree->SetBranchAddress("sigMET",    &mc_sigMET);
	mctree->SetBranchAddress("phoSigma",  &mc_phoSigma);
	mctree->SetBranchAddress("phoChIso",  &mc_phoChIso);
	if(useMC){
		mctree->SetBranchAddress("mcPID",   &mc_mcPID);
		mctree->SetBranchAddress("mcEta",   &mc_mcEta);
		mctree->SetBranchAddress("mcPhi",   &mc_mcPhi);
		mctree->SetBranchAddress("mcPt",    &mc_mcPt);
		mctree->SetBranchAddress("mcMomPID",&mc_mcMomPID);
	}

  TChain *datatree = new TChain("signalTree");
  datatree->Add("/uscms_data/d3/mengleis/test/plot_hadron_ISR.root");
	float data_phoEt(0);
	float data_sigMET(0);
	float data_phoSigma(0);
	float data_phoChIso(0);
	float data_dilepMass(0);
	datatree->SetBranchAddress("phoEt",     &data_phoEt);
	datatree->SetBranchAddress("sigMET",    &data_sigMET);
	datatree->SetBranchAddress("phoSigma",  &data_phoSigma);
	datatree->SetBranchAddress("phoChIso",  &data_phoChIso);
	datatree->SetBranchAddress("dilepMass", &data_dilepMass);

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
		if(mc_sigMET < METLOWCUT || mc_sigMET > METHIGHCUT)continue;

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

		bool isTrueTemplate = (mc_phoSigma <= 0.0103 && !isFake);
		bool isInTargetRegion = (mc_phoSigma <= 0.0103 && mc_phoChIso < 1.295);

		if(mc_phoSigma <= 0.0103 && isFake)mc_fake->Fill(mc_phoChIso);
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
					if(mc_phoChIso < 1.295)nsbPassIso[iLower][iUpper] += 1;
				}
			}
		}					
	}

	if(nMCtarget > 0)mcfakerate = 1.0*nMCfaketarget/nMCtarget;

	for(unsigned ievt(0); ievt < datatree->GetEntries(); ievt++){
		datatree->GetEntry(ievt);
		if(data_phoEt < lowercut || data_phoEt > uppercut )continue;
		if(data_sigMET < METLOWCUT || data_sigMET > METHIGHCUT)continue;
		if(data_dilepMass < 80 || data_dilepMass > 100)continue;
		//if(data_dilepMass > 80 || data_dilepMass < 30)continue;

		if(data_phoSigma <= 0.0103){
			h_target->Fill(data_phoChIso);
			totalden+=1;
			if(data_phoChIso < 1.295)nden+=1;
		}
		else{
			for(unsigned iUpper(0); iUpper<nUpper; iUpper++){
				for(unsigned iLower(0); iLower<nLower; iLower++){
					if(data_phoSigma > SigmaCutLower[iLower] && data_phoSigma < SigmaCutUpper[iUpper]){
						h_bg[iLower][iUpper]->Fill(data_phoChIso);
						totalnum[iLower][iUpper] +=1;
						if(data_phoChIso < 1.295)nnum[iLower][iUpper] +=1;
					}
				}
			}
		}

	}
	for(unsigned iUpper(0); iUpper<nUpper; iUpper++){
    for(unsigned iLower(0); iLower<nLower; iLower++){
			cantemplate[iLower][iUpper]->cd(1);
			h_target->SetLineColor(kRed);
			h_target->Draw();
			h_bg[iLower][iUpper]->Draw("same");
			
			cantemplate[iLower][iUpper]->cd(2);
			mc_sig->Draw();
		}
	}
			

	std::cout << "target " << h_target->Integral(1,10) << std::endl;
	h_target->Sumw2();
	for(unsigned iUpper(0); iUpper<nUpper; iUpper++){
		for(unsigned iLower(0); iLower<nLower; iLower++){
			h_bg[iLower][iUpper]->Sumw2();
		}
	}

	for(unsigned iUpper(0); iUpper<nUpper; iUpper++){
		for(unsigned iLower(0); iLower<nLower; iLower++){
			double iteratorUncertainty(1), lastUpdatePurity(1);
			double frac_passSigma = nsbPassIso[iLower][iUpper]/(1.0*mc_sbcontamination[iLower][iUpper]->GetEntries());
			if(mc_sbcontamination[iLower][iUpper]->GetEntries()==0)frac_passSigma = 0;
			double phoPurity(0);
			templateCorrFactor[iLower][iUpper] = 0.0;
			if(doIterate){
				while(iteratorUncertainty > 0.01){

					TH1F *hist_new=(TH1F*)h_bg[iLower][iUpper]->Clone();
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
				result[iLower][iUpper] = (TH1F*)fitter[iLower][iUpper]->GetPlot();
				result[iLower][iUpper]->SetMinimum(10);
				fitter[iLower][iUpper]->GetResult(0, fracSig, errorSig);
				fitter[iLower][iUpper]->GetResult(1, fracBkg, errorBkg);
				can[iLower][iUpper]->cd();
				gStyle->SetOptStat(0);
				hname.str("");
				hname << lowername << "< pt <" << uppername; 
				h_target->SetTitle(hname.str().c_str());
				h_target->SetMinimum(h_target->GetBinContent(10));
				result[iLower][iUpper]->SetMinimum(100);
				gPad->SetLogy();
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
		 
				if(SigmaCutLower[iLower] == 0.0103 && SigmaCutUpper[iUpper]==0.014){
          hadfrac =  num*fracBkg/den;
					fittingError = (numerror/num + denerror/den + errorBkg/fracBkg)*num*fracBkg/den ;
					TLegend *leg = new TLegend(0.6,0.75,0.87,0.9);
					leg->AddEntry(h_target, "observed");
					leg->AddEntry(result[iLower][iUpper], "true photons");
					leg->AddEntry(mc_predict[iLower][iUpper], "hadrons");
					leg->Draw("same");
					hname.str("");
					hname << "frac-" << lowername << "-" << uppername << "_ChIso-MuonEG-July22.pdf";
					can[iLower][iUpper]->SaveAs(hname.str().c_str());
				}
				fracHad2D->Fill(SigmaCutLower[iLower]+0.00005, SigmaCutUpper[iUpper]+0.00005, num*fracBkg/den);
				fracHad1D->Fill(num*fracBkg/den);
				if(lowestfrac > num*fracBkg/den)lowestfrac = num*fracBkg/den;
				if(highestfrac < num*fracBkg/den)highestfrac = num*fracBkg/den; 
			}
    }
  }

	if(hadfrac < 0.001)hadfrac = (highestfrac+lowestfrac)/2;
  systematicError = highestfrac - lowestfrac;
	myfile << " " << lowername << " " << uppername << " " << hadfrac << " "; 
  myfile << sqrt(fittingError*fittingError + systematicError*systematicError) << " ";
  //myfile << systematicError <<"	" << std::endl;
  myfile << systematicError << "	" << " ";
  myfile << mcfakerate  << std::endl;
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
  outputfile->Write();
  return 1;
}
