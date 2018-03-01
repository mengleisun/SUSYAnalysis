#include "TChain.h"
#include<iostream>
#include "TTree.h"
#include "TFile.h"
#include "TFileCollection.h"
#include "TMath.h"

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
	     default: isEle = false; break;
	   }
  }
  else isEle = false;

  return isEle;
}

float DeltaR(float eta1,float phi1,float eta2,float phi2)
{
	float deltaPhi = TMath::Abs(phi1-phi2);
	float deltaEta = eta1-eta2;
	if(deltaPhi > TMath::Pi())
	deltaPhi = TMath::TwoPi() - deltaPhi;
		return TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
}

void skimMC(){
  
  TChain *es= new TChain("eeTree");
  es->Add("/uscms_data/d3/mengleis/FullStatusOct/plot_R9_DY.root");

  int nevt = es->GetEntries();

	float tagPt(0);
	float tagEta(0);
	float tagPhi(0);
	float tagR9(0);
	float probeEt(0);
	float probeEta(0);
	float probePhi(0);
	float probeR9(0);
	bool  probeMatchLeading;
	bool  probeMatchTrailing;
	float invmass;
	std::vector<int>   *mcPID=0;
	std::vector<float> *mcEta=0;
	std::vector<float> *mcPhi=0;
	std::vector<float> *mcPt=0;
	std::vector<int>   *mcMomPID=0;
	std::vector<int>   *mcGMomPID=0;

	es->SetBranchAddress("tagPt",                 &tagPt);
	es->SetBranchAddress("tagEta",                &tagEta);
	es->SetBranchAddress("tagPhi",                &tagPhi);
	es->SetBranchAddress("tagR9",                 &tagR9);
	es->SetBranchAddress("probeEt",            		&probeEt);
	es->SetBranchAddress("probeEta",           		&probeEta);
	es->SetBranchAddress("probePhi",           		&probePhi);
	es->SetBranchAddress("probeR9",            		&probeR9);
	es->SetBranchAddress("probeMatchLeading",  		&probeMatchLeading);
	es->SetBranchAddress("probeMatchTrailing", 		&probeMatchTrailing);
	es->SetBranchAddress("invmass",               &invmass);
	es->SetBranchAddress("mcPID",			   					&mcPID);
	es->SetBranchAddress("mcEta",			   					&mcEta);
	es->SetBranchAddress("mcPhi",			   					&mcPhi);
	es->SetBranchAddress("mcPt",				   				&mcPt);
	es->SetBranchAddress("mcMomPID",			   			&mcMomPID);
	es->SetBranchAddress("mcGMomPID",		   				&mcGMomPID);

  TFile *outputfile = TFile::Open("/uscms_data/d3/mengleis/FullStatusOct/plot_R9_DY_mcTrue.root","RECREATE");
  TTree *tree_out = es->CloneTree(0);
	bool mcTrue;
	tree_out->Branch("mcTrue", &mcTrue);

  std::cout << "total " << nevt << std::endl;
  for(unsigned ievt(0); ievt < nevt; ievt++){
    es->GetEntry(ievt); 

    if (ievt%100000==0) std::cout << " -- Processing event " << ievt << std::endl;
	  bool isZee(false);
	  double mindRtag(0.3), mindRprobe(0.3);
	  unsigned tagIndex(0), probeIndex(0);
	  for(unsigned iMC(0); iMC<mcPID->size(); iMC++){
			double dR1 = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], tagEta, tagPhi);
			double dR2 = DeltaR((*mcEta)[iMC], (*mcPhi)[iMC], probeEta,probePhi);
			double dE1 = fabs((*mcPt)[iMC] - tagPt)/tagPt;
			double dE2 = fabs((*mcPt)[iMC] - probeEt)/probeEt;
			if(dR1 < mindRtag && dE1 < 0.1){mindRtag=dR1; tagIndex=iMC;}
			if(dR2 < mindRprobe && dE2 < 0.1){mindRprobe=dR2; probeIndex=iMC;}
	  }
	  if(mindRtag < 0.1 && mindRprobe < 0.1){
			bool isZe(false),isZg(false);
			isZe = isElectron(fabs((*mcPID)[tagIndex]), fabs((*mcMomPID)[tagIndex]));
			isZg = isElectron(fabs((*mcPID)[probeIndex]), fabs((*mcMomPID)[probeIndex]));
  		if(isZe && isZg)isZee=true; 
    }
		if(isZee)mcTrue = true;
		else mcTrue = false;

    tree_out->Fill();
  }

  outputfile->Write();
  outputfile->Close();
}

