#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include "TTree.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TF1.h"
#include "TMath.h"
#include "TLegend.h"
#include "../../include/tdrstyle.C"

void analysis_SUSY(){//main  

  //TFile *file = TFile::Open("/uscms_data/d3/mengleis/resTree_TChiWg.root");
  TFile *file = TFile::Open("/uscms_data/d3/mengleis/resTree_T5Wg.root");
  TTree *tree = (TTree*)file->Get("SUSYtree");
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
	float EleRelIso(0);
  float EledEtaIn(0);
  float EledPhiIn(0);
  float EleD0(0);
  float EleDz(0);
  float EleooEmooP(0);
  bool  ElePassID;
  float Mchagino(0);
  //float Mgluino(0);
  float Mneutralino(0);
  float dRPhoEle(0);
  float Invmass(0);
  double MT_(0), ThreeBodyMass_(0);
	float HT(0);
	float sigMET(0);

  tree->SetBranchAddress("MT",&MT_);
  tree->SetBranchAddress("ThreeBodyMass"  ,&ThreeBodyMass_);
  tree->SetBranchAddress("mcPhotonEt"     ,&mcPhotonEt); 
  tree->SetBranchAddress("mcPhotonEta"    ,&mcPhotonEta);
  tree->SetBranchAddress("mcPhotonPhi"    ,&mcPhotonPhi);
  tree->SetBranchAddress("recoPhotonEt"   ,&recoPhotonEt);
  tree->SetBranchAddress("recoPhotonEta"  ,&recoPhotonEta);
  tree->SetBranchAddress("recoPhotonPhi"  ,&recoPhotonPhi);
  tree->SetBranchAddress("PhoR9"          ,&PhoR9);
  tree->SetBranchAddress("PhoHoverE"      ,&PhoHoverE);
  tree->SetBranchAddress("PhoSigma"       ,&PhoSigma);
  tree->SetBranchAddress("PhoChIso"       ,&PhoChIso);
  tree->SetBranchAddress("PhoNeuIso"      ,&PhoNeuIso);
  tree->SetBranchAddress("PhoPhoIso"      ,&PhoPhoIso);
  tree->SetBranchAddress("PhoPassID"      ,&PhoPassID);
  tree->SetBranchAddress("mcElePt"        ,&mcElePt);
  tree->SetBranchAddress("mcEleEta"       ,&mcEleEta);
  tree->SetBranchAddress("mcElePhi"       ,&mcElePhi);
  tree->SetBranchAddress("recoEleEt"      ,&recoEleEt);
  tree->SetBranchAddress("recoEleEta"     ,&recoEleEta);
  tree->SetBranchAddress("recoElePhi"     ,&recoElePhi);
  tree->SetBranchAddress("EleR9"          ,&EleR9);
  tree->SetBranchAddress("EleHoverE"      ,&EleHoverE);
  tree->SetBranchAddress("EleSigma"       ,&EleSigma);
  tree->SetBranchAddress("EleChIso"       ,&EleChIso);
  tree->SetBranchAddress("EleNeuIso"      ,&EleNeuIso);
  tree->SetBranchAddress("ElePhoIso"      ,&ElePhoIso);
  tree->SetBranchAddress("EleMiniIso"     ,&EleMiniIso);
  tree->SetBranchAddress("EleRelIso"      ,&EleRelIso);
  tree->SetBranchAddress("EledEtaIn"      ,&EledEtaIn);
  tree->SetBranchAddress("EledPhiIn"      ,&EledPhiIn);
  tree->SetBranchAddress("EleD0"          ,&EleD0);
  tree->SetBranchAddress("EleDz"          ,&EleDz);
  tree->SetBranchAddress("EleooEmooP"     ,&EleooEmooP);
  tree->SetBranchAddress("ElePassID"      ,&ElePassID);
  tree->SetBranchAddress("Mchagino"       ,&Mchagino);
  //tree->SetBranchAddress("Mgluino"       ,&Mgluino);
  tree->SetBranchAddress("Mneutralino",   &Mneutralino);
  tree->SetBranchAddress("dRPhoEle"       ,&dRPhoEle);
  tree->SetBranchAddress("Invmass"        ,&Invmass);    
  tree->SetBranchAddress("HT"             ,&HT);
 	tree->SetBranchAddress("sigMET"         ,&sigMET);

	setTDRStyle();    
	std::ostringstream histname;

	TFile xSecFile("susyxSec.root");
	//TH1F *p_crosssection = (TH1F*)xSecFile.Get("p_charginoSec"); 
	TH1F *p_crosssection = (TH1F*)xSecFile.Get("p_gluinoxSec"); 

	TFile *outputfile = TFile::Open("plot_T5Wg.root","RECREATE");
	outputfile->cd();

  //********** eventcount ********************//
    double eventcount[3][3][25][25];
      for(unsigned i(0); i < 3; i++)
          for(unsigned j(0); j < 3; j++)
							for(unsigned m(0); m < 25; m++)
								for(unsigned n(0); n < 25; n++)
                eventcount[i][j][m][n] = 0;
  
  

	TH2D *p_SUSYMass = new TH2D("p_SUSYMass","",26,800,2100,26,800,2100);
	TH2D *p_SUSYMt = new TH2D("p_SUSYMt","",26,800,2100,26,800,2100);

	long gridcount[25][25];
	TH1F *p_MET[25][25];
	for(unsigned i(0); i<25; i++){
		for(unsigned j(0); j<25; j++){
			gridcount[i][j] = 0;
			histname.str("");
			histname << "p_MET_";
			if(i==0)histname << "Glu800-";
			else if(i==1)histname << "Glu900-";
			else if(i > 1)histname << "Glu" << 1000 + (i-1)*50 << "-";
			if(j==0)histname << "Chi800";
			else if(j==1)histname << "Chi900";
			else if(j > 1)histname << "Chi" << 1000 + (j-1)*50; 
			p_MET[i][j] = new TH1F(histname.str().c_str(), histname.str().c_str(),200, 0, 2000);
			p_MET[i][j]->SetLineColor(j);
		}
	}
	
	TH1F *p_eventcount[25][25];
	for(unsigned i(0); i<25; i++){
		for(unsigned j(0); j<25; j++){
			histname.str("");
			histname << "p_eventcount_";
			if(i==0)histname << "Glu800_";
			else if(i==1)histname << "Glu900_";
			else if(i > 1)histname << "Glu" << 1000 + (i-1)*50 << "_";
			if(j==0)histname << "Chi800";
			else if(j==1)histname << "Chi900";
			else if(j > 1)histname << "Chi" << 1000 + (j-1)*50; 
			p_eventcount[i][j] = new TH1F(histname.str().c_str(), histname.str().c_str(),9,0,9);
			p_eventcount[i][j]->SetLineColor(j);
		}
	}
	

	float massSUSY[25];
	massSUSY[0] = 100, massSUSY[1] = 900;
	for(unsigned ibin(2); ibin < 25; ibin++)massSUSY[ibin] = 900 + ibin*50;

	int nEvt = tree->GetEntries();
	for(unsigned ievt(0); ievt < nEvt; ievt++){
		tree->GetEntry(ievt);
    if (ievt%1000000==0) std::cout << " -- Processing event " << ievt << std::endl;
		if(Mchagino >0 && Mgluino >0){

			p_SUSYMass->Fill(Mchagino, Mgluino);
			double xsweight = p_crosssection->GetBinContent(p_crosssection->FindBin(Mgluino))*35.8*1000;
			if(sigMET > 100 && mcElePt > 25 && mcPhotonEt > 35 && MT_ > 1 && MT_ < 100)p_SUSYMt->Fill(Mchagino, Mgluino,xsweight);
			for(unsigned i(0); i < 24; i++){
				if(Mgluino > massSUSY[i] && Mgluino <= massSUSY[i+1]){
					for(unsigned j(0); j < 24; j++){
						if(Mchagino > massSUSY[j] && Mchagino <= massSUSY[j+1]){
			
								gridcount[i][j]+=1;					
								if(mcElePt < 25 || mcPhotonEt < 35 || fabs(mcPhotonEta) < 1.444)continue;
							  if(MT_ > 1 && MT_ < 100 && sigMET > 0)p_MET[i][j]->Fill(sigMET,xsweight); 	

								if(sigMET > 120 && sigMET <= 200){
								  if(HT  < 100)eventcount[0][0][i][j] += xsweight;
									else if(HT  > 100 && HT < 400)eventcount[0][1][i][j] += xsweight;
									else if(HT >= 400)eventcount[0][2][i][j] += xsweight;
								}
								else if(sigMET > 200 && sigMET <= 300){
									if(HT  < 100)eventcount[1][0][i][j] += xsweight;
									else if(HT  > 100 && HT < 400)eventcount[1][1][i][j] += xsweight;
									else if(HT >= 400)eventcount[1][2][i][j] += xsweight;
								}
								else if(sigMET > 300){
									if(HT  < 100)eventcount[2][0][i][j] += xsweight;
									else if(HT  > 100 && HT < 400)eventcount[2][1][i][j] += xsweight;
									else if(HT >= 400)eventcount[2][2][i][j] += xsweight;
								}
			

						}
					}
				}
			}
		}

	}


	for(int i(1); i<= 26; i++){
		for(int j(1); j<=26; j++){
			if(p_SUSYMt->GetBinContent(i,j) > 0)p_SUSYMt->SetBinContent(i, j, p_SUSYMt->GetBinContent(i,j)/p_SUSYMass->GetBinContent(i,j));
		}
	}

	for(int i(0); i < 25; i++){
		for(int j(0); j < 25; j++){
			if(p_MET[i][j]->Integral(1,100) > 0)p_MET[i][j]->Scale(1.0/gridcount[i][j]);
		}
	}

	for(unsigned mass1(0); mass1 < 24; mass1++){
			for(unsigned mass2(0); mass2 < 24; mass2++){
			if(p_SUSYMass->GetBinContent(mass2+1,mass1+1) > 1){
					 if(massSUSY[mass2] > massSUSY[mass1] )continue;
		        std::cout << massSUSY[mass1]  << " "  << massSUSY[mass2]  << std::endl;
					  for(unsigned i(0); i < 3; i++){
    					for(unsigned j(0); j < 3; j++){
      						std::cout << " " << eventcount[i][j][mass1][mass2]/p_SUSYMass->GetBinContent(mass2+1,mass1+1);
									p_eventcount[mass1][mass2]->Fill(i*3+j, eventcount[i][j][mass1][mass2]/p_SUSYMass->GetBinContent(mass2+1,mass1+1));
							}
						}
					std::cout << std::endl;
				}
			}
	}

  xSecFile.Close();

	TH1F *p_databkg = new TH1F("p_databkg","",9,0,9);
	p_databkg->SetBinContent(1,277.58);
	p_databkg->SetBinContent(2,392.8);
	p_databkg->SetBinContent(3,83.42);
	p_databkg->SetBinContent(4,32.9);
	p_databkg->SetBinContent(5,55.16);
	p_databkg->SetBinContent(6,24.74);
	p_databkg->SetBinContent(7,8.77);
	p_databkg->SetBinContent(8,6.58);
	p_databkg->SetBinContent(9,9.31);

	p_SUSYMt->Write();
	p_SUSYMass->Write();
	for(unsigned mass1(0); mass1 < 24; mass1++){
			for(unsigned mass2(0); mass2 < 24; mass2++){	
				p_eventcount[mass1][mass2]->Write();
				p_MET[mass1][mass2]->Write();
			}
	}
	p_databkg->Write();
	outputfile->Write();
	outputfile->Close();

}


