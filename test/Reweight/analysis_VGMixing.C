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
#include "../../include/analysis_scalefactor.h"

void analysis_VGMixing(){
	
  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

  TFile *outputfile = TFile::Open("/uscms_data/d3/mengleis/test/resTree_ZGamma_ZG_inclusive.root","RECREATE");

  TChain *mgtree;
  mgtree = new TChain("ZTree");
  mgtree->Add("/uscms_data/d3/mengleis/test/resTree_ZGamma_ZG.root");
  mgtree->Add("/uscms_data/d3/mengleis/test/resTree_ZGamma_ZG130.root");

	float MCweight(0);
	float phoEt(0);
	mgtree->SetBranchAddress("MCweight",    &MCweight);
  mgtree->SetBranchAddress("phoEt",       &phoEt);

  TTree *tree_out = mgtree->CloneTree(0);
  for(unsigned ievt(0); ievt < mgtree->GetEntries(); ievt++){
    mgtree->GetEntry(ievt); 
    if( MCweight > 0.01 && phoEt >= 140)continue;
		else if(MCweight < 0.01 && phoEt < 140)continue;
    tree_out->Fill();
  }

	
  outputfile->Write();
  outputfile->Close();
}


