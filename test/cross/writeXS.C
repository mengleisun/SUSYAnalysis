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

void writeXS(){
	std::ifstream T6WG_file("CrossSectionT6WG.txt");
	std::ofstream newfile("newCrossSectionT6WG.txt");

	double susymass(0);
	double xsvalue(0);
	double xserror(0);
	if(T6WG_file.is_open()){
  	for(int i(0); i<500; i++){ 
			T6WG_file >> susymass >> xsvalue >> xserror;
			newfile << susymass << " " << xsvalue*1000.0 << " " << xsvalue*1000.0*xserror*0.01 << std::endl; 
	  }
	}
	T6WG_file.close();

}

