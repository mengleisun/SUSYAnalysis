#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include "TH1F.h"
#include "TCanvas.h"

#define NTOY 9000

void plotsys(){
	TH1F *p_slope = new TH1F("p_slope","slope",100,0,1000);
	TH1F *p_const = new TH1F("p_const","const",100,-4000,0);
	TH1F *p_index = new TH1F("p_index","index",100,-1,0);
	double scalefactor(0), ptslope(0), ptconstant(0), ptindex(0), vtxconst(0), vtxslope(0);
	std::ifstream elefake_toyfile("ToyFakeRate.txt");
	if(elefake_toyfile.is_open()){
  	for(int i(0); i<NTOY; i++){ 
			elefake_toyfile >> scalefactor >> ptslope >> ptconstant >> ptindex >>  vtxconst >> vtxslope;
			p_slope->Fill(ptslope);
			p_const->Fill(ptconstant);
			p_index->Fill(ptindex);
	  }
	}
	elefake_toyfile.close();
	TCanvas *can = new TCanvas("can","",800,800);
	can->Divide(2,2);
	can->cd(1);
	p_slope->Draw();
	can->cd(2);
	p_const->Draw();
	can->cd(3);
	p_index->Draw();
}
