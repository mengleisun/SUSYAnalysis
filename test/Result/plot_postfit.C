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
#include "../../include/tdrstyle.C"

int channel = 2; // 1 = eg, 2 = mg

void plot_postfit(){//main  

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");
	setTDRStyle();

	std::ostringstream signame;  signame.str("");
	std::ostringstream elename;  elename.str("");
	std::ostringstream jetname;  jetname.str("");
	std::ostringstream qcdname;  qcdname.str("");
	std::ostringstream VGname;   VGname.str("");
	std::ostringstream rarename; rarename.str(""); 

	if(channel == 1){
			signame <<  "vetoDiPhoData/signalTree_egamma_signal.root";
			elename <<  "vetoDiPhoData/signalTree_egamma_eleBkg.root";  
			jetname <<  "vetoDiPhoData/signalTree_egamma_jetbkg.root";
			qcdname <<  "vetoDiPhoData/signalTree_egamma_qcd.root";
			VGname  <<  "vetoDiPhoData/signalTree_egamma_VGBkg.root";
			rarename << "vetoDiPhoData/signalTree_egamma_rareBkg.root";
	}
	else if(channel == 2){
			signame << "vetoDiPhoData/signalTree_mg_signal.root";
			elename << "vetoDiPhoData/signalTree_mg_eleBkg.root";  
			jetname << "vetoDiPhoData/signalTree_mg_jetbkg.root";
			qcdname << "vetoDiPhoData/signalTree_mg_qcd.root";
			VGname  << "vetoDiPhoData/signalTree_mg_VGBkg.root";
			rarename <<"vetoDiPhoData/signalTree_mg_rareBkg.root";
	}

	TFile *file_sig = TFile::Open(signame.str().c_str());
	TFile *file_ele = TFile::Open(elename.str().c_str());
	TFile *file_jet = TFile::Open(jetname.str().c_str());
	TFile *file_qcd = TFile::Open(qcdname.str().c_str());
	TFile *file_VG  = TFile::Open(VGname.str().c_str());
	TFile *file_rare= TFile::Open(rarename.str().c_str()); 

	TFile *file_t5 = TFile::Open("signalTree_T5WG.root");
	TFile *file_tchi=TFile::Open("signalTree_TChiWG.root");
	TH1D *p_t5wg_MET_signal_1700_1000= (TH1D*)file_t5->Get("p_t5wg_MET_signal_1700_1000_mg");
	TH1D *p_tchiwg_MET_signal_800  = (TH1D*)file_tchi->Get("p_tchiwg_MET_signal_mg");
	TH1D *p_t5wg_HT_signal_1700_1000= (TH1D*)file_t5->Get("p_t5wg_HT_signal_1700_1000_mg");
	TH1D *p_tchiwg_HT_signal_800  = (TH1D*)file_tchi->Get("p_tchiwg_HT_signal_mg");
	TH1D *p_t5wg_PhoEt_signal_1700_1000= (TH1D*)file_t5->Get("p_t5wg_PhoEt_signal_1700_1000_mg");
	TH1D *p_tchiwg_PhoEt_signal_800  = (TH1D*)file_tchi->Get("p_tchiwg_PhoEt_signal_mg");

	Double_t bkgEtBins[]={35,40,45,50,55,60,65,70,75,80, 85,90,95,100,105,110,115,120,125,130, 135,140,146,152,158,164,170,177,184,192, 200,208,216,224,232,240,250,260,275,290, 305,325,345,370,400,500,800};
	int nBkgEtBins= sizeof(bkgEtBins)/sizeof(bkgEtBins[0]) -1;
	Double_t bkgPtBins[]={25,30,35,40,45,50,55,60,65,70,75,80, 85,90,95,100,105,110,115,120,125,130, 135,140,146,152,158,164,170,177,184,192, 200,208,216,224,232,240,250,260,275,290, 305,325,345,370,400,500,800};
	int nBkgPtBins= sizeof(bkgPtBins)/sizeof(bkgPtBins[0])-1;
	Double_t bkgMETBins[]={0,40,60,80,100,120,140,160,180,210,240,280,320,400,600,1000};
	int nBkgMETBins= sizeof(bkgMETBins)/sizeof(bkgMETBins[0]) -1;
	Double_t bkgMtBins[]={0,20,40,60,80,100,120,140,160,180,200,300,400,500,600,800,1000};
	int nBkgMtBins= sizeof(bkgMtBins)/sizeof(bkgMtBins[0]) -1;
	Double_t bkgHTBins[]={0,40,60,80,100,120,140,160,180,200,225,250,275,300,340,380,420,500,600,1000};
	int nBkgHTBins= sizeof(bkgHTBins)/sizeof(bkgHTBins[0]) -1;

	TGraphErrors *error_PhoEt = new TGraphErrors(nBkgEtBins);
	TGraphErrors *error_LepPt = new TGraphErrors(nBkgPtBins);
	TGraphErrors *error_MET = new TGraphErrors(nBkgMETBins); 
	TGraphErrors *error_Mt = new TGraphErrors(nBkgMtBins); 
	TGraphErrors *error_HT = new TGraphErrors(nBkgHTBins);
	TGraphErrors *ratioerror_PhoEt = new TGraphErrors(nBkgEtBins);
	TGraphErrors *ratioerror_LepPt = new TGraphErrors(nBkgPtBins);
	TGraphErrors *ratioerror_MET = new TGraphErrors(nBkgMETBins); 
	TGraphErrors *ratioerror_Mt = new TGraphErrors(nBkgMtBins); 
	TGraphErrors *ratioerror_HT = new TGraphErrors(nBkgHTBins);

	TH1F *p_allPhoEt = (TH1F*)file_sig->Get("p_PhoEt");
	TH1F *p_allLepPt = (TH1F*)file_sig->Get("p_LepPt");
	TH1F *p_allMET   = (TH1F*)file_sig->Get("p_MET");
	TH1F *p_allMt    = (TH1F*)file_sig->Get("p_Mt");
	TH1F *p_allHT  = (TH1F*)file_sig->Get("p_HT");
	TH1F *p_allPU    = (TH1F*)file_sig->Get("p_PU");

	p_allPhoEt->SetMarkerColor(1);
	p_allLepPt->SetMarkerColor(1);
	p_allMET->SetMarkerColor(1);
	p_allMt->SetMarkerColor(1);  
	p_allHT->SetMarkerColor(1);
	p_allPU->SetMarkerColor(1);
	p_allPhoEt->SetTitle("");
	p_allLepPt->SetTitle("");
	p_allMET->SetTitle("");
	p_allMt->SetTitle("");
	p_allHT->SetTitle("");
	p_allPU->SetTitle("");
	p_allPhoEt->GetYaxis()->SetTitle("Events / bin");
	p_allLepPt->GetYaxis()->SetTitle("Events / bin");
	p_allMET->GetYaxis()->SetTitle("Events / bin");
	p_allMt->GetYaxis()->SetTitle("Events / bin");
	p_allHT->GetYaxis()->SetTitle("Events / bin");
	p_allPU->GetYaxis()->SetTitle("Events / bin");
	
	TH1F *p_elePhoEt = (TH1F*)file_ele->Get("p_PhoEt");
	TH1F *p_eleLepPt = (TH1F*)file_ele->Get("p_LepPt");
	TH1F *p_eleMET   = (TH1F*)file_ele->Get("p_MET");
	TH1F *p_eleMt    = (TH1F*)file_ele->Get("p_Mt");
	TH1F *p_eleHT  = (TH1F*)file_ele->Get("p_HT");
	TH1F *p_elePU    = (TH1F*)file_ele->Get("p_PU");

	TH1F *p_jetPhoEt = (TH1F*)file_jet->Get("p_PhoEt");
	TH1F *p_jetLepPt = (TH1F*)file_jet->Get("p_LepPt");
	TH1F *p_jetMET   = (TH1F*)file_jet->Get("p_MET");
	TH1F *p_jetMt    = (TH1F*)file_jet->Get("p_Mt");
	TH1F *p_jetHT  = (TH1F*)file_jet->Get("p_HT");
	TH1F *p_jetPU    = (TH1F*)file_jet->Get("p_PU");

	TH1F *p_qcdPhoEt = (TH1F*)file_qcd->Get("p_PhoEt");
	TH1F *p_qcdLepPt = (TH1F*)file_qcd->Get("p_LepPt");
	TH1F *p_qcdMET   = (TH1F*)file_qcd->Get("p_MET");
	TH1F *p_qcdMt    = (TH1F*)file_qcd->Get("p_Mt");
	TH1F *p_qcdHT  = (TH1F*)file_qcd->Get("p_HT");
	TH1F *p_qcdPU    = (TH1F*)file_qcd->Get("p_PU");

	TH1F *p_VGPhoEt = (TH1F*)file_VG->Get("p_PhoEt");
	TH1F *p_VGLepPt = (TH1F*)file_VG->Get("p_LepPt");
	TH1F *p_VGMET   = (TH1F*)file_VG->Get("p_MET");
	TH1F *p_VGMt    = (TH1F*)file_VG->Get("p_Mt");
	TH1F *p_VGHT  = (TH1F*)file_VG->Get("p_HT");
	TH1F *p_VGPU    = (TH1F*)file_VG->Get("p_PU");

	TH1F *p_rarePhoEt = (TH1F*)file_rare->Get("p_PhoEt");
	TH1F *p_rareLepPt = (TH1F*)file_rare->Get("p_LepPt");
	TH1F *p_rareMET   = (TH1F*)file_rare->Get("p_MET");
	TH1F *p_rareMt    = (TH1F*)file_rare->Get("p_Mt");
	TH1F *p_rareHT  = (TH1F*)file_rare->Get("p_HT");
	TH1F *p_rarePU    = (TH1F*)file_rare->Get("p_PU");

	int binnumber;
 binnumber = p_allPhoEt->GetSize()-2;  p_allPhoEt->SetBinContent(binnumber,  p_allPhoEt->GetBinContent(binnumber) +   p_allPhoEt->GetBinContent(binnumber+1) );   
 binnumber = p_allLepPt->GetSize()-2;  p_allLepPt->SetBinContent(binnumber,  p_allLepPt->GetBinContent(binnumber) +   p_allLepPt->GetBinContent(binnumber+1) ); 
 binnumber = p_allMET->GetSize()-2;    p_allMET->SetBinContent(binnumber,    p_allMET->GetBinContent(binnumber) +     p_allMET->GetBinContent(binnumber+1) );   
 binnumber = p_allMt->GetSize()-2;     p_allMt->SetBinContent(binnumber,     p_allMt->GetBinContent(binnumber) +      p_allMt->GetBinContent(binnumber+1) );    
 binnumber = p_allHT->GetSize()-2;     p_allHT->SetBinContent(binnumber,     p_allHT->GetBinContent(binnumber) +      p_allHT->GetBinContent(binnumber+1) );    
                                                                                                                          
 binnumber = p_elePhoEt->GetSize()-2;  p_elePhoEt->SetBinContent(binnumber,  p_elePhoEt->GetBinContent(binnumber) +   p_elePhoEt->GetBinContent(binnumber+1) ); 
 binnumber = p_eleLepPt->GetSize()-2;  p_eleLepPt->SetBinContent(binnumber,  p_eleLepPt->GetBinContent(binnumber) +   p_eleLepPt->GetBinContent(binnumber+1) ); 
 binnumber = p_eleMET->GetSize()-2;    p_eleMET->SetBinContent(binnumber,    p_eleMET->GetBinContent(binnumber) +     p_eleMET->GetBinContent(binnumber+1) );   
 binnumber = p_eleMt->GetSize()-2;     p_eleMt->SetBinContent(binnumber,     p_eleMt->GetBinContent(binnumber) +      p_eleMt->GetBinContent(binnumber+1) );    
 binnumber = p_eleHT->GetSize()-2;     p_eleHT->SetBinContent(binnumber,     p_eleHT->GetBinContent(binnumber) +      p_eleHT->GetBinContent(binnumber+1) );    
                                                                                                                                                                           
 binnumber = p_jetPhoEt->GetSize()-2;  p_jetPhoEt->SetBinContent(binnumber,  p_jetPhoEt->GetBinContent(binnumber) +   p_jetPhoEt->GetBinContent(binnumber+1) );             
 binnumber = p_jetLepPt->GetSize()-2;  p_jetLepPt->SetBinContent(binnumber,  p_jetLepPt->GetBinContent(binnumber) +   p_jetLepPt->GetBinContent(binnumber+1) );             
 binnumber = p_jetMET->GetSize()-2;    p_jetMET->SetBinContent(binnumber,    p_jetMET->GetBinContent(binnumber) +     p_jetMET->GetBinContent(binnumber+1) );               
 binnumber = p_jetMt->GetSize()-2;     p_jetMt->SetBinContent(binnumber,     p_jetMt->GetBinContent(binnumber) +      p_jetMt->GetBinContent(binnumber+1) );                
 binnumber = p_jetHT->GetSize()-2;     p_jetHT->SetBinContent(binnumber,     p_jetHT->GetBinContent(binnumber) +      p_jetHT->GetBinContent(binnumber+1) );  
                                                                                                                                                                           
 binnumber = p_qcdPhoEt->GetSize()-2;  p_qcdPhoEt->SetBinContent(binnumber,  p_qcdPhoEt->GetBinContent(binnumber) +   p_qcdPhoEt->GetBinContent(binnumber+1) );  
 binnumber = p_qcdLepPt->GetSize()-2;  p_qcdLepPt->SetBinContent(binnumber,  p_qcdLepPt->GetBinContent(binnumber) +   p_qcdLepPt->GetBinContent(binnumber+1) );  
 binnumber = p_qcdMET->GetSize()-2;    p_qcdMET->SetBinContent(binnumber,    p_qcdMET->GetBinContent(binnumber) +     p_qcdMET->GetBinContent(binnumber+1) );    
 binnumber = p_qcdMt->GetSize()-2;     p_qcdMt->SetBinContent(binnumber,     p_qcdMt->GetBinContent(binnumber) +      p_qcdMt->GetBinContent(binnumber+1) );     
 binnumber = p_qcdHT->GetSize()-2;     p_qcdHT->SetBinContent(binnumber,     p_qcdHT->GetBinContent(binnumber) +      p_qcdHT->GetBinContent(binnumber+1) );  
                                                                                                                                                                           
 binnumber = p_VGPhoEt->GetSize()-2;   p_VGPhoEt->SetBinContent(binnumber,   p_VGPhoEt->GetBinContent(binnumber) +    p_VGPhoEt->GetBinContent(binnumber+1) );   
 binnumber = p_VGLepPt->GetSize()-2;   p_VGLepPt->SetBinContent(binnumber,   p_VGLepPt->GetBinContent(binnumber) +    p_VGLepPt->GetBinContent(binnumber+1) );   
 binnumber = p_VGMET->GetSize()-2;     p_VGMET->SetBinContent(binnumber,     p_VGMET->GetBinContent(binnumber) +      p_VGMET->GetBinContent(binnumber+1) );     
 binnumber = p_VGMt->GetSize()-2;      p_VGMt->SetBinContent(binnumber,      p_VGMt->GetBinContent(binnumber) +       p_VGMt->GetBinContent(binnumber+1) );      
 binnumber = p_VGHT->GetSize()-2;      p_VGHT->SetBinContent(binnumber,      p_VGHT->GetBinContent(binnumber) +       p_VGHT->GetBinContent(binnumber+1) );   
                                                                                                                                                                           
 binnumber = p_rarePhoEt->GetSize()-2; p_rarePhoEt->SetBinContent(binnumber, p_rarePhoEt->GetBinContent(binnumber) +  p_rarePhoEt->GetBinContent(binnumber+1) );            
 binnumber = p_rareLepPt->GetSize()-2; p_rareLepPt->SetBinContent(binnumber, p_rareLepPt->GetBinContent(binnumber) +  p_rareLepPt->GetBinContent(binnumber+1) );            
 binnumber = p_rareMET->GetSize()-2;   p_rareMET->SetBinContent(binnumber,   p_rareMET->GetBinContent(binnumber) +    p_rareMET->GetBinContent(binnumber+1) );              
 binnumber = p_rareMt->GetSize()-2;    p_rareMt->SetBinContent(binnumber,    p_rareMt->GetBinContent(binnumber) +     p_rareMt->GetBinContent(binnumber+1) );               
 binnumber = p_rareHT->GetSize()-2;    p_rareHT->SetBinContent(binnumber,    p_rareHT->GetBinContent(binnumber) +     p_rareHT->GetBinContent(binnumber+1) );                  
// Errors
 binnumber = p_allPhoEt->GetSize()-2;  p_allPhoEt->SetBinError(binnumber,  p_allPhoEt->GetBinError(binnumber) +   p_allPhoEt->GetBinError(binnumber+1) );   
 binnumber = p_allLepPt->GetSize()-2;  p_allLepPt->SetBinError(binnumber,  p_allLepPt->GetBinError(binnumber) +   p_allLepPt->GetBinError(binnumber+1) ); 
 binnumber = p_allMET->GetSize()-2;    p_allMET->SetBinError(binnumber,    p_allMET->GetBinError(binnumber) +     p_allMET->GetBinError(binnumber+1) );   
 binnumber = p_allMt->GetSize()-2;     p_allMt->SetBinError(binnumber,     p_allMt->GetBinError(binnumber) +      p_allMt->GetBinError(binnumber+1) );    
 binnumber = p_allHT->GetSize()-2;     p_allHT->SetBinError(binnumber,     p_allHT->GetBinError(binnumber) +      p_allHT->GetBinError(binnumber+1) );    
                                                                                                                          
 binnumber = p_elePhoEt->GetSize()-2;  p_elePhoEt->SetBinError(binnumber,  p_elePhoEt->GetBinError(binnumber) +   p_elePhoEt->GetBinError(binnumber+1) ); 
 binnumber = p_eleLepPt->GetSize()-2;  p_eleLepPt->SetBinError(binnumber,  p_eleLepPt->GetBinError(binnumber) +   p_eleLepPt->GetBinError(binnumber+1) ); 
 binnumber = p_eleMET->GetSize()-2;    p_eleMET->SetBinError(binnumber,    p_eleMET->GetBinError(binnumber) +     p_eleMET->GetBinError(binnumber+1) );   
 binnumber = p_eleMt->GetSize()-2;     p_eleMt->SetBinError(binnumber,     p_eleMt->GetBinError(binnumber) +      p_eleMt->GetBinError(binnumber+1) );    
 binnumber = p_eleHT->GetSize()-2;     p_eleHT->SetBinError(binnumber,     p_eleHT->GetBinError(binnumber) +      p_eleHT->GetBinError(binnumber+1) );    
                                                                                                                                                                           
 binnumber = p_jetPhoEt->GetSize()-2;  p_jetPhoEt->SetBinError(binnumber,  p_jetPhoEt->GetBinError(binnumber) +   p_jetPhoEt->GetBinError(binnumber+1) );             
 binnumber = p_jetLepPt->GetSize()-2;  p_jetLepPt->SetBinError(binnumber,  p_jetLepPt->GetBinError(binnumber) +   p_jetLepPt->GetBinError(binnumber+1) );             
 binnumber = p_jetMET->GetSize()-2;    p_jetMET->SetBinError(binnumber,    p_jetMET->GetBinError(binnumber) +     p_jetMET->GetBinError(binnumber+1) );               
 binnumber = p_jetMt->GetSize()-2;     p_jetMt->SetBinError(binnumber,     p_jetMt->GetBinError(binnumber) +      p_jetMt->GetBinError(binnumber+1) );                
 binnumber = p_jetHT->GetSize()-2;     p_jetHT->SetBinError(binnumber,     p_jetHT->GetBinError(binnumber) +      p_jetHT->GetBinError(binnumber+1) );  
                                                                                                                                                                           
 binnumber = p_qcdPhoEt->GetSize()-2;  p_qcdPhoEt->SetBinError(binnumber,  p_qcdPhoEt->GetBinError(binnumber) +   p_qcdPhoEt->GetBinError(binnumber+1) );  
 binnumber = p_qcdLepPt->GetSize()-2;  p_qcdLepPt->SetBinError(binnumber,  p_qcdLepPt->GetBinError(binnumber) +   p_qcdLepPt->GetBinError(binnumber+1) );  
 binnumber = p_qcdMET->GetSize()-2;    p_qcdMET->SetBinError(binnumber,    p_qcdMET->GetBinError(binnumber) +     p_qcdMET->GetBinError(binnumber+1) );    
 binnumber = p_qcdMt->GetSize()-2;     p_qcdMt->SetBinError(binnumber,     p_qcdMt->GetBinError(binnumber) +      p_qcdMt->GetBinError(binnumber+1) );     
 binnumber = p_qcdHT->GetSize()-2;     p_qcdHT->SetBinError(binnumber,     p_qcdHT->GetBinError(binnumber) +      p_qcdHT->GetBinError(binnumber+1) );  
                                                                                                                                                                           
 binnumber = p_VGPhoEt->GetSize()-2;   p_VGPhoEt->SetBinError(binnumber,   p_VGPhoEt->GetBinError(binnumber) +    p_VGPhoEt->GetBinError(binnumber+1) );   
 binnumber = p_VGLepPt->GetSize()-2;   p_VGLepPt->SetBinError(binnumber,   p_VGLepPt->GetBinError(binnumber) +    p_VGLepPt->GetBinError(binnumber+1) );   
 binnumber = p_VGMET->GetSize()-2;     p_VGMET->SetBinError(binnumber,     p_VGMET->GetBinError(binnumber) +      p_VGMET->GetBinError(binnumber+1) );     
 binnumber = p_VGMt->GetSize()-2;      p_VGMt->SetBinError(binnumber,      p_VGMt->GetBinError(binnumber) +       p_VGMt->GetBinError(binnumber+1) );      
 binnumber = p_VGHT->GetSize()-2;      p_VGHT->SetBinError(binnumber,      p_VGHT->GetBinError(binnumber) +       p_VGHT->GetBinError(binnumber+1) );   
                                                                                                                                                                           
 binnumber = p_rarePhoEt->GetSize()-2; p_rarePhoEt->SetBinError(binnumber, p_rarePhoEt->GetBinError(binnumber) +  p_rarePhoEt->GetBinError(binnumber+1) );            
 binnumber = p_rareLepPt->GetSize()-2; p_rareLepPt->SetBinError(binnumber, p_rareLepPt->GetBinError(binnumber) +  p_rareLepPt->GetBinError(binnumber+1) );            
 binnumber = p_rareMET->GetSize()-2;   p_rareMET->SetBinError(binnumber,   p_rareMET->GetBinError(binnumber) +    p_rareMET->GetBinError(binnumber+1) );              
 binnumber = p_rareMt->GetSize()-2;    p_rareMt->SetBinError(binnumber,    p_rareMt->GetBinError(binnumber) +     p_rareMt->GetBinError(binnumber+1) );               
 binnumber = p_rareHT->GetSize()-2;    p_rareHT->SetBinError(binnumber,    p_rareHT->GetBinError(binnumber) +     p_rareHT->GetBinError(binnumber+1) );                  


	binnumber = p_t5wg_MET_signal_1700_1000->GetSize()-2; p_t5wg_MET_signal_1700_1000->SetBinContent(binnumber, p_t5wg_MET_signal_1700_1000->GetBinContent(binnumber) + p_t5wg_MET_signal_1700_1000->GetBinContent(binnumber+1) );
	binnumber = p_tchiwg_MET_signal_800->GetSize()-2; p_tchiwg_MET_signal_800->SetBinContent(binnumber, p_tchiwg_MET_signal_800->GetBinContent(binnumber) + p_tchiwg_MET_signal_800->GetBinContent(binnumber+1) );
	binnumber = p_t5wg_HT_signal_1700_1000->GetSize()-2; p_t5wg_HT_signal_1700_1000->SetBinContent(binnumber, p_t5wg_HT_signal_1700_1000->GetBinContent(binnumber) + p_t5wg_HT_signal_1700_1000->GetBinContent(binnumber+1) );
	binnumber = p_tchiwg_HT_signal_800->GetSize()-2; p_tchiwg_HT_signal_800->SetBinContent(binnumber, p_tchiwg_HT_signal_800->GetBinContent(binnumber) + p_tchiwg_HT_signal_800->GetBinContent(binnumber+1) );

	std::ostringstream etplot;  etplot.str("");
	std::ostringstream ptplot;  ptplot.str("");
	std::ostringstream metplot; metplot.str("");
	std::ostringstream mtplot;  mtplot.str("");
	std::ostringstream htplot;  htplot.str("");
	
	if(channel == 1){
			etplot << "SIGNAL_egamma_2016ReMiniAOD_pt.pdf";     
			ptplot << "SIGNAL_egamma_2016ReMiniAOD_leppt.pdf";	 
			metplot << "SIGNAL_egamma_2016ReMiniAOD_met.pdf";   
			mtplot << "SIGNAL_egamma_2016ReMiniAOD_mt.pdf";     
			htplot << "SIGNAL_egamma_2016ReMiniAOD_ht.pdf";	   
	}                                                      
	else if(channel == 2){
			etplot << "SIGNAL_mg_2016ReMiniAOD_pt.pdf";
			ptplot << "SIGNAL_mg_2016ReMiniAOD_leppt.pdf";		
			metplot << "SIGNAL_mg_2016ReMiniAOD_met.pdf";
			mtplot << "SIGNAL_mg_2016ReMiniAOD_mt.pdf";
			htplot << "SIGNAL_mg_2016ReMiniAOD_ht.pdf";	
	}
	// ******** Mt ************************//
	gStyle->SetOptStat(0);
	TCanvas *c_mt = new TCanvas("Mt", "Mt",600,600);
	setCanvas(c_mt); 
	c_mt->cd();
	TPad *mt_pad1 = new TPad("mt_pad1", "mt_pad1", 0, 0.3, 1, 1.0);
	setTopPad(mt_pad1); 
	mt_pad1->Draw();  
	mt_pad1->cd();  
	gPad->SetLogy();
	p_allMt->SetMinimum(0.05);
  p_allMt->SetMaximum(50*p_allMt->GetBinContent(p_allMt->GetMaximumBin()));
	p_allMt->GetXaxis()->SetRangeUser(0,1000);
	p_allMt->SetLineColor(1);
	p_allMt->SetMarkerStyle(20);
	p_allMt->Draw("P");
	p_VGMt->SetFillStyle(1001);
	p_VGMt->SetLineColor(kGreen-7);
	p_VGMt->SetFillColor(kGreen-7);
	p_rareMt->SetFillStyle(1001);
	p_rareMt->SetLineColor(kAzure-9);
	p_rareMt->SetFillColor(kAzure-9);
	p_qcdMt->SetFillStyle(1001);
	p_qcdMt->SetLineColor(kYellow-9);
	p_qcdMt->SetFillColor(kYellow-9);
	p_eleMt->SetFillStyle(1001);
	p_eleMt->SetLineColor(kMagenta-9);
	p_eleMt->SetFillColor(kMagenta-9);
	p_jetMt->SetFillStyle(1001);
	p_jetMt->SetLineColor(kOrange-9);
	p_jetMt->SetFillColor(kOrange-9);
	p_eleMt->Add(p_rareMt); // ele 2nd
	p_jetMt->Add(p_eleMt);  // jet 3rd
	p_qcdMt->Add(p_jetMt);  // qcd 4th
	p_VGMt->Add(p_qcdMt);   // VG  5th
	for(int ibin(1); ibin < p_VGMt->GetSize(); ibin++){
		error_Mt->SetPoint(ibin-1,p_VGMt->GetBinCenter(ibin), p_VGMt->GetBinContent(ibin));
		float prederror = p_VGMt->GetBinError(ibin);
//		prederror += p_eleMt->GetBinError(ibin);
//		prederror += p_jetMt->GetBinError(ibin);
//		prederror += p_qcdMt->GetBinError(ibin);
//		prederror += p_rareMt->GetBinError(ibin);
		error_Mt->SetPointError(ibin-1,(p_VGMt->GetBinLowEdge(ibin+1)-p_VGMt->GetBinLowEdge(ibin))/2,prederror);
		ratioerror_Mt->SetPoint(ibin-1,p_VGMt->GetBinCenter(ibin), 1); 
		ratioerror_Mt->SetPointError(ibin-1,(p_VGMt->GetBinLowEdge(ibin+1)-p_VGMt->GetBinLowEdge(ibin))/2, prederror/p_VGMt->GetBinContent(ibin)); 
	}
	p_VGMt->Draw("hist same");
	p_qcdMt->Draw("hist same");
	p_jetMt->Draw("hist same");
	p_eleMt->Draw("hist same");
	p_rareMt->Draw("hist same");
  error_Mt->SetFillColor(12);
  error_Mt->SetFillStyle(3345);
	//error_Mt->Draw("E2 same");
	TLegend *leg_mt =  new TLegend(0.35,0.55,0.9,0.83);
	leg_mt->SetNColumns(2);
	leg_mt->SetFillStyle(0);
	leg_mt->SetBorderSize(0);
	leg_mt->SetFillColor(0);
	p_rareMt->SetMarkerSize(0);
	p_eleMt->SetMarkerSize(0);
	p_jetMt->SetMarkerSize(0);
	p_qcdMt->SetMarkerSize(0);
	p_VGMt->SetMarkerSize(0);
	ratioerror_Mt->SetMarkerSize(0);
	ratioerror_Mt->SetLineWidth(0);
  p_t5wg_MET_signal_1700_1000->SetLineColor(9);
  p_t5wg_MET_signal_1700_1000->SetLineWidth(4);
	p_tchiwg_MET_signal_800->SetLineColor(28);
	p_tchiwg_MET_signal_800->SetLineStyle(2);
	p_tchiwg_MET_signal_800->SetLineWidth(4);
	leg_mt->AddEntry(p_allMt,"observed","epl");
	leg_mt->AddEntry(p_rareMt,"t#bar{t}#gamma / WW#gamma / WZ#gamma");
	leg_mt->AddEntry(p_eleMt,"e#rightarrow#gamma fakes");
	leg_mt->AddEntry(p_jetMt,"j#rightarrow#gamma fakes");
	leg_mt->AddEntry(p_qcdMt,"fake lepton");
	leg_mt->AddEntry(p_VGMt, "W#gamma / Z#gamma");
  leg_mt->AddEntry(p_t5wg_MET_signal_1700_1000, "T5Wg");
	leg_mt->AddEntry(p_tchiwg_MET_signal_800, "TChiWg");
	leg_mt->AddEntry(ratioerror_Mt, "Unc");
	leg_mt->Draw("same");
	p_allMt->Draw("E same");
 	gPad->RedrawAxis();
  CMS_lumi( mt_pad1, 11 );

	c_mt->cd();
	TPad *mt_pad2 = new TPad("mt_pad2", "mt_pad2", 0, 0, 1, 0.3);
	mt_pad2->SetBottomMargin(0.4);
	mt_pad2->Draw();
	mt_pad2->cd();
  TLine *flatratio_mt = new TLine(0,1,800,1);
	TH1F *ratio_mt=(TH1F*)p_allMt->Clone("transfer factor");
	ratio_mt->SetMarkerStyle(20);
	ratio_mt->SetLineColor(kBlack);
	ratio_mt->GetXaxis()->SetRangeUser(0,800);
	ratio_mt->GetYaxis()->SetRangeUser(0,2);
	ratio_mt->GetYaxis()->SetNdivisions(504);
	ratio_mt->SetMinimum(0);
	ratio_mt->SetMaximum(2);
	ratio_mt->Divide(p_VGMt);
	ratio_mt->SetTitle("");
	ratio_mt->GetYaxis()->SetTitle("#frac{Obs.}{Bkg.}");
	ratio_mt->Draw();
	ratioerror_Mt->SetFillColor(12);
	ratioerror_Mt->SetFillStyle(3345);
	ratioerror_Mt->Draw("E2 same");
	ratio_mt->Draw("same");
	flatratio_mt->Draw("same");
	c_mt->SaveAs(mtplot.str().c_str());
 
  gStyle->SetOptStat(0);
	TCanvas *c_pt = new TCanvas("Photon_Pt", "Photon P_{T}",600,600);
	setCanvas(c_pt); 
	c_pt->cd();
	TPad *pt_pad1 = new TPad("pt_pad1", "pt_pad1", 0, 0.3, 1, 1.0);
	setTopPad(pt_pad1); 
	pt_pad1->Draw();  
	pt_pad1->cd();  
	gPad->SetLogy();
	p_allPhoEt->SetMaximum(50*p_allPhoEt->GetBinContent(p_allPhoEt->GetMaximumBin()));
	p_allPhoEt->SetMinimum(0.5);
	p_allPhoEt->GetXaxis()->SetRangeUser(35,800);
	p_allPhoEt->SetLineColor(1);
	p_allPhoEt->SetMarkerStyle(20);
	p_allPhoEt->Draw("P");
	p_VGPhoEt->SetFillStyle(1001);
	p_VGPhoEt->SetLineColor(kGreen-7);
	p_VGPhoEt->SetFillColor(kGreen-7);
	p_rarePhoEt->SetFillStyle(1001);
	p_rarePhoEt->SetLineColor(kAzure-9);
	p_rarePhoEt->SetFillColor(kAzure-9);
	p_qcdPhoEt->SetFillStyle(1001);
	p_qcdPhoEt->SetLineColor(kYellow-9);
	p_qcdPhoEt->SetFillColor(kYellow-9);
	p_elePhoEt->SetFillStyle(1001);
	p_elePhoEt->SetLineColor(kMagenta-9);
	p_elePhoEt->SetFillColor(kMagenta-9);
	p_jetPhoEt->SetFillStyle(1001);
	p_jetPhoEt->SetLineColor(kOrange-9);
	p_jetPhoEt->SetFillColor(kOrange-9);
	p_elePhoEt->Add(p_rarePhoEt); // ele 2nd
	p_jetPhoEt->Add(p_elePhoEt);  // jet 3rd
	p_qcdPhoEt->Add(p_jetPhoEt);  // qcd 4th
	p_VGPhoEt->Add(p_qcdPhoEt);   // VG  5th
	p_VGPhoEt->Sumw2();
	for(int ibin(1); ibin < p_VGPhoEt->GetSize(); ibin++){
		error_PhoEt->SetPoint(ibin-1,p_VGPhoEt->GetBinCenter(ibin), p_VGPhoEt->GetBinContent(ibin));
		float prederror = p_VGPhoEt->GetBinError(ibin);
		//prederror += p_elePhoEt->GetBinError(ibin);
		//prederror += p_jetPhoEt->GetBinError(ibin);
		//prederror += p_qcdPhoEt->GetBinError(ibin);
		//prederror += p_rarePhoEt->GetBinError(ibin)*0.6;
		error_PhoEt->SetPointError(ibin-1,(p_VGPhoEt->GetBinLowEdge(ibin+1)-p_VGPhoEt->GetBinLowEdge(ibin))/2,prederror);
		std::cout << p_elePhoEt->GetBinError(ibin) << " " << p_jetPhoEt->GetBinError(ibin) << " " << p_qcdPhoEt->GetBinError(ibin) << " " << p_rarePhoEt->GetBinError(ibin)  << std::endl;
		ratioerror_PhoEt->SetPoint(ibin-1,p_VGPhoEt->GetBinCenter(ibin), 1); 
		ratioerror_PhoEt->SetPointError(ibin-1,(p_VGPhoEt->GetBinLowEdge(ibin+1)-p_VGPhoEt->GetBinLowEdge(ibin))/2, prederror/p_VGPhoEt->GetBinContent(ibin)); 
	}
	p_VGPhoEt->Draw("hist same");
	p_qcdPhoEt->Draw("hist same");
	p_jetPhoEt->Draw("hist same");
	p_elePhoEt->Draw("hist same");
	p_rarePhoEt->Draw("hist same");
  error_PhoEt->SetFillColor(12);
  error_PhoEt->SetFillStyle(3345);
	error_PhoEt->Draw("E2 same");
	leg_mt->Draw("same");
	p_allPhoEt->Draw("E same");
  p_t5wg_PhoEt_signal_1700_1000->SetLineColor(9);
  p_t5wg_PhoEt_signal_1700_1000->SetLineWidth(4);
  p_t5wg_PhoEt_signal_1700_1000->Draw("same");
	p_tchiwg_PhoEt_signal_800->SetLineColor(28);
	p_tchiwg_PhoEt_signal_800->SetLineStyle(2);
	p_tchiwg_PhoEt_signal_800->SetLineWidth(4);
	p_tchiwg_PhoEt_signal_800->Draw("same");
  TLatex chantex;
  chantex.SetNDC();
  chantex.SetTextFont(42);
  chantex.SetTextSize(0.05);    
  chantex.DrawLatex(0.4,0.87,"M_{T} > 100 GeV, p_{T}^{miss} > 120 GeV");
  chantex.SetTextSize(0.07);    
	chantex.DrawLatex(0.8,0.5,"(d)");
 	gPad->RedrawAxis();
  CMS_lumi( pt_pad1, 11 );

	c_pt->cd();
	TPad *pt_pad2 = new TPad("pt_pad2", "pt_pad2", 0, 0, 1, 0.3);
	pt_pad2->SetBottomMargin(0.4);
	pt_pad2->Draw();
	pt_pad2->cd();
  TLine *flatratio = new TLine(35,1,800,1);
	TH1F *ratio=(TH1F*)p_allPhoEt->Clone("transfer factor");
	ratio->SetMinimum(0);
	ratio->SetMaximum(2);
	ratio->SetMarkerStyle(20);
	ratio->SetLineColor(kBlack);
	ratio->Divide(p_VGPhoEt);
	ratio->SetTitle("");
	ratio->GetYaxis()->SetTitle("#frac{Obs.}{Bkg.}");
	ratio->GetYaxis()->SetNdivisions(504);
	ratio->Draw();
	ratioerror_PhoEt->SetFillColor(12);
	ratioerror_PhoEt->SetFillStyle(3345);
	ratioerror_PhoEt->Draw("E2 same");
	ratio->Draw("same");
	flatratio->Draw("same");
	c_pt->SaveAs(etplot.str().c_str());

// ******** MET ************************//
	gStyle->SetOptStat(0);
	TCanvas *c_met = new TCanvas("MET", "MET",600,600);
	setCanvas(c_met); 
	c_met->cd();
	TPad *met_pad1 = new TPad("met_pad1", "met_pad1", 0, 0.3, 1, 1.0);
	setTopPad(met_pad1); 
	met_pad1->Draw();  
	met_pad1->cd();  
	gPad->SetLogy();
	p_allMET->GetYaxis()->SetRangeUser(0.05, 50*p_allMET->GetBinContent(p_allMET->GetMaximumBin()));
	p_allMET->SetMinimum(0.5);
	p_allMET->GetXaxis()->SetRangeUser(0,600);
	p_allMET->SetLineColor(1);
	p_allMET->SetMarkerStyle(20);
	p_allMET->Draw("P");
	p_VGMET->SetFillStyle(1001);
	p_VGMET->SetLineColor(kGreen-7);
	p_VGMET->SetFillColor(kGreen-7);
	p_rareMET->SetFillStyle(1001);
	p_rareMET->SetLineColor(kAzure-9);
	p_rareMET->SetFillColor(kAzure-9);
	p_qcdMET->SetFillStyle(1001);
	p_qcdMET->SetLineColor(kYellow-9);
	p_qcdMET->SetFillColor(kYellow-9);
	p_eleMET->SetFillStyle(1001);
	p_eleMET->SetLineColor(kMagenta-9);
	p_eleMET->SetFillColor(kMagenta-9);
	p_jetMET->SetFillStyle(1001);
	p_jetMET->SetLineColor(kOrange-9);
	p_jetMET->SetFillColor(kOrange-9);
	p_eleMET->Add(p_rareMET); // ele 2nd
	p_jetMET->Add(p_eleMET);  // jet 3rd
	p_qcdMET->Add(p_jetMET);  // qcd 4th
	p_VGMET->Add(p_qcdMET);   // VG  5th
	p_VGMET->Sumw2();
	for(int ibin(1); ibin < p_VGMET->GetSize(); ibin++){
		float prederror = p_VGMET->GetBinError(ibin);
	//	prederror += p_eleMET->GetBinError(ibin);
	//	prederror += p_jetMET->GetBinError(ibin);
	//	prederror += p_qcdMET->GetBinError(ibin);
	//	prederror += p_rareMET->GetBinError(ibin);
		error_MET->SetPoint(ibin-1,p_VGMET->GetBinCenter(ibin), p_VGMET->GetBinContent(ibin));
		error_MET->SetPointError(ibin-1,(p_VGMET->GetBinLowEdge(ibin+1)-p_VGMET->GetBinLowEdge(ibin))/2,prederror);
		ratioerror_MET->SetPoint(ibin-1,p_VGMET->GetBinCenter(ibin), 1); 
		ratioerror_MET->SetPointError(ibin-1,(p_VGMET->GetBinLowEdge(ibin+1)-p_VGMET->GetBinLowEdge(ibin))/2, prederror/p_VGMET->GetBinContent(ibin)); 
	}
	p_VGMET->Draw("hist same");
	p_qcdMET->Draw("hist same");
	p_jetMET->Draw("hist same");
	p_eleMET->Draw("hist same");
	p_rareMET->Draw("hist same");
  error_MET->SetFillColor(12);
  error_MET->SetFillStyle(3345);
	error_MET->Draw("E2 same");
//	TLegend *leg_met =  new TLegend(0.5,0.65,0.9,0.9);
//	leg_met->SetFillStyle(0);
//	gStyle->SetLegendBorderSize(1);
//	gStyle->SetLegendFillColor(0);
//	leg_met->AddEntry(p_allMET,"observed (MT < 100 GeV)");
//	leg_met->AddEntry(p_rareMET,"t#bar{t}#gamma/WW#gamma/WZ#gamma");
//	leg_met->AddEntry(p_eleMET,"e->#gamma fake");
//	leg_met->AddEntry(p_jetMET,"j->#gamma fake");
//	leg_met->AddEntry(p_qcdMET,"j->e fake");
//	leg_met->AddEntry(p_VGMET, "WG/ZG");
	leg_mt->Draw("same");
	p_allMET->Draw("E same");
  chantex.SetTextSize(0.05);    
  chantex.DrawLatex(0.4,0.87,"M_{T} > 100 GeV");
  chantex.SetTextSize(0.07);    
	chantex.DrawLatex(0.8,0.5,"(b)");
 	gPad->RedrawAxis();
  CMS_lumi( met_pad1, 11 );

  p_t5wg_MET_signal_1700_1000->SetLineColor(9);
  p_t5wg_MET_signal_1700_1000->SetLineWidth(4);
  p_t5wg_MET_signal_1700_1000->Draw("same");
	p_tchiwg_MET_signal_800->SetLineColor(28);
	p_tchiwg_MET_signal_800->SetLineStyle(2);
	p_tchiwg_MET_signal_800->SetLineWidth(4);
	p_tchiwg_MET_signal_800->Draw("same");
//	TLegend *leg_susy =  new TLegend(0.1,0.75,0.5,0.9);
//	leg_susy->SetFillStyle(0);
//	gStyle->SetLegendBorderSize(1);
//	gStyle->SetLegendFillColor(0);
//  leg_susy->AddEntry(p_t5wg_MET_signal_1800_1000,"T5WG,1800,1000");
//  leg_susy->AddEntry(p_t5wg_MET_signal_1800_200, "T5WG,1800,200");
//	leg_susy->AddEntry(p_tchiwg_MET_signal_800, "TChiWG, 800");
//	leg_susy->AddEntry(p_tchiwg_MET_signal_1000,"TChiWG, 1000");
//	leg_susy->Draw("same");

	c_met->cd();
	TPad *met_pad2 = new TPad("met_pad2", "met_pad2", 0, 0, 1, 0.3);
	met_pad2->SetBottomMargin(0.4);
	met_pad2->Draw();
	met_pad2->cd();
  TLine *flatratio_met = new TLine(0,1,600,1);
	TH1F *ratio_met=(TH1F*)p_allMET->Clone("transfer factor");
	ratio_met->GetXaxis()->SetRangeUser(0,600);
	ratio_met->GetYaxis()->SetNdivisions(504);
	ratio_met->SetLineColor(kBlack);
	ratio_met->SetMarkerStyle(20);
	ratio_met->Divide(p_VGMET);
	ratio_met->SetTitle("");
	ratio_met->GetYaxis()->SetTitle("#frac{Obs.}{Bkg.}");
	ratio_met->GetYaxis()->SetRangeUser(0,2.3);
	ratio_met->Draw();
	ratioerror_MET->SetFillColor(12);
	ratioerror_MET->SetFillStyle(3345);
	ratioerror_MET->Draw("E2 same");
	ratio_met->Draw("same");
	flatratio_met->Draw("same");
	c_met->SaveAs(metplot.str().c_str());


// ******** LepPt ************************//
	gStyle->SetOptStat(0);
	TCanvas *c_leppt = new TCanvas("LepPt", "LepPt",600,600);
	setCanvas(c_leppt); 
	c_leppt->cd();
	TPad *leppt_pad1 = new TPad("leppt_pad1", "leppt_pad1", 0, 0.3, 1, 1.0);
	setTopPad(leppt_pad1); 
	leppt_pad1->Draw();  
	leppt_pad1->cd();  
	gPad->SetLogy();
	p_allLepPt->SetMinimum(p_rareLepPt->GetBinContent(p_rareLepPt->GetSize()));
	p_allLepPt->GetXaxis()->SetRangeUser(35,800);
	p_allLepPt->SetLineColor(1);
	p_allLepPt->SetMarkerStyle(20);
	p_allLepPt->Draw("P");
	p_VGLepPt->SetFillStyle(1001);
	p_VGLepPt->SetLineColor(kGreen-7);
	p_VGLepPt->SetFillColor(kGreen-7);
	p_rareLepPt->SetFillStyle(1001);
	p_rareLepPt->SetLineColor(kAzure-9);
	p_rareLepPt->SetFillColor(kAzure-9);
	p_qcdLepPt->SetFillStyle(1001);
	p_qcdLepPt->SetLineColor(kYellow-9);
	p_qcdLepPt->SetFillColor(kYellow-9);
	p_eleLepPt->SetFillStyle(1001);
	p_eleLepPt->SetLineColor(kMagenta-9);
	p_eleLepPt->SetFillColor(kMagenta-9);
	p_jetLepPt->SetFillStyle(1001);
	p_jetLepPt->SetLineColor(kOrange-9);
	p_jetLepPt->SetFillColor(kOrange-9);
	p_eleLepPt->Add(p_rareLepPt); // ele 2nd
	p_jetLepPt->Add(p_eleLepPt);  // jet 3rd
	p_qcdLepPt->Add(p_jetLepPt);  // qcd 4th
	p_VGLepPt->Add(p_qcdLepPt);   // VG  5th
	for(int ibin(1); ibin < p_VGLepPt->GetSize(); ibin++){
		error_LepPt->SetPoint(ibin-1,p_VGLepPt->GetBinCenter(ibin), p_VGLepPt->GetBinContent(ibin));
		float prederror = p_VGLepPt->GetBinError(ibin);
//		prederror += p_eleLepPt->GetBinError(ibin);
//		prederror += p_jetLepPt->GetBinError(ibin);
//		prederror += p_qcdLepPt->GetBinError(ibin);
//		prederror += p_rareLepPt->GetBinError(ibin);
		error_LepPt->SetPointError(ibin-1,(p_VGLepPt->GetBinLowEdge(ibin+1)-p_VGLepPt->GetBinLowEdge(ibin))/2,prederror);
		ratioerror_LepPt->SetPoint(ibin-1,p_VGLepPt->GetBinCenter(ibin), 1); 
		ratioerror_LepPt->SetPointError(ibin-1,(p_VGLepPt->GetBinLowEdge(ibin+1)-p_VGLepPt->GetBinLowEdge(ibin))/2, prederror/p_VGLepPt->GetBinContent(ibin)); 
	}
	p_VGLepPt->Draw("hist same");
	p_qcdLepPt->Draw("hist same");
	p_jetLepPt->Draw("hist same");
	p_eleLepPt->Draw("hist same");
	p_rareLepPt->Draw("hist same");
  error_LepPt->SetFillColor(12);
  error_LepPt->SetFillStyle(3345);
	//error_LepPt->Draw("E2 same");
	TLegend *leg_leppt =  new TLegend(0.6,0.75,0.9,0.9);
	leg_leppt->SetFillStyle(0);
	gStyle->SetLegendBorderSize(1);
	gStyle->SetLegendFillColor(0);
	leg_leppt->AddEntry(p_allLepPt,"observed (LepPt < 100 GeV)");
	leg_leppt->AddEntry(p_rareLepPt,"t#bar{t}#gamma/WW#gamma/WZ#gamma");
	leg_leppt->AddEntry(p_eleLepPt,"e->#gamma fake");
	leg_leppt->AddEntry(p_jetLepPt,"j->#gamma fake");
	leg_leppt->AddEntry(p_qcdLepPt,"l fakes");
	leg_leppt->AddEntry(p_VGLepPt, "W#gamma/Z#gamma");
	leg_leppt->Draw("same");
	p_allLepPt->Draw("E same");
 	gPad->RedrawAxis();
  CMS_lumi( leppt_pad1, 11 );

	c_leppt->cd();
	TPad *leppt_pad2 = new TPad("leppt_pad2", "leppt_pad2", 0, 0, 1, 0.3);
	leppt_pad2->SetBottomMargin(0.4);
	leppt_pad2->Draw();
	leppt_pad2->cd();
  TLine *flatratio_leppt = new TLine(25,1,800,1);
	TH1F *ratio_leppt=(TH1F*)p_allLepPt->Clone("transfer factor");
	ratio_leppt->SetMarkerStyle(20);
	ratio_leppt->SetLineColor(kBlack);
	ratio_leppt->Divide(p_VGLepPt);
	ratio_leppt->SetTitle("");
	ratio_leppt->GetYaxis()->SetTitle("#frac{Obs.}{Bkg.}");
	ratio_leppt->GetYaxis()->SetRangeUser(0,2);
	ratio_leppt->Draw();
	ratioerror_LepPt->SetFillColor(12);
	ratioerror_LepPt->SetFillStyle(3345);
	ratioerror_LepPt->Draw("E2 same");
	ratio_leppt->Draw("same");
	flatratio_leppt->Draw("same");
	c_leppt->SaveAs(ptplot.str().c_str());

// ******** HT ************************//
	gStyle->SetOptStat(0);
	TCanvas *c_HT = new TCanvas("HT", "HT",600,600);
	setCanvas(c_HT); 
	c_HT->cd();
	TPad *HT_pad1 = new TPad("HT_pad1", "HT_pad1", 0, 0.3, 1, 1.0);
	setTopPad(HT_pad1); 
	HT_pad1->Draw();  
	HT_pad1->cd();  
	gPad->SetLogy();
	p_allHT->SetMinimum(0.5);
	p_allHT->SetMaximum(50*p_allHT->GetBinContent(p_allHT->GetMaximumBin()));
	p_allHT->SetLineColor(1);
	p_allHT->SetMarkerStyle(20);
	p_allHT->Draw("P");
	p_VGHT->SetFillStyle(1001);
	p_VGHT->SetLineColor(kGreen-7);
	p_VGHT->SetFillColor(kGreen-7);
	p_rareHT->SetFillStyle(1001);
	p_rareHT->SetLineColor(kAzure-9);
	p_rareHT->SetFillColor(kAzure-9);
	p_qcdHT->SetFillStyle(1001);
	p_qcdHT->SetLineColor(kYellow-9);
	p_qcdHT->SetFillColor(kYellow-9);
	p_eleHT->SetFillStyle(1001);
	p_eleHT->SetLineColor(kMagenta-9);
	p_eleHT->SetFillColor(kMagenta-9);
	p_jetHT->SetFillStyle(1001);
	p_jetHT->SetLineColor(kOrange-9);
	p_jetHT->SetFillColor(kOrange-9);
	p_eleHT->Add(p_rareHT); // ele 2nd
	p_jetHT->Add(p_eleHT);  // jet 3rd
	p_qcdHT->Add(p_jetHT);  // qcd 4th
	p_VGHT->Add(p_qcdHT);   // VG  5th
	for(int ibin(1); ibin < p_VGHT->GetSize(); ibin++){
		error_HT->SetPoint(ibin-1,p_VGHT->GetBinCenter(ibin), p_VGHT->GetBinContent(ibin));
		float prederror = p_VGHT->GetBinError(ibin);
		error_HT->SetPointError(ibin-1,(p_VGHT->GetBinLowEdge(ibin+1)-p_VGHT->GetBinLowEdge(ibin))/2,prederror);
		ratioerror_HT->SetPoint(ibin-1,p_VGHT->GetBinCenter(ibin), 1); 
		ratioerror_HT->SetPointError(ibin-1,(p_VGHT->GetBinLowEdge(ibin+1)-p_VGHT->GetBinLowEdge(ibin))/2, prederror/p_VGHT->GetBinContent(ibin)); 
	}
	p_VGHT->Draw("hist same");
	p_qcdHT->Draw("hist same");
	p_jetHT->Draw("hist same");
	p_eleHT->Draw("hist same");
	p_rareHT->Draw("hist same");
	error_HT->SetFillColor(12);
	error_HT->SetFillStyle(3345);
	error_HT->Draw("E2 same");
	leg_mt->Draw("same");
	p_allHT->Draw("E same");
  chantex.SetTextSize(0.05);    
  chantex.DrawLatex(0.4,0.87,"M_{T} > 100 GeV, p_{T}^{miss} > 120 GeV");
  chantex.SetTextSize(0.07);    
	chantex.DrawLatex(0.8,0.5,"(f)");
 	gPad->RedrawAxis();
  CMS_lumi( HT_pad1, 11 );

  p_t5wg_HT_signal_1700_1000->SetLineColor(9);
  p_t5wg_HT_signal_1700_1000->SetLineWidth(4);
  p_t5wg_HT_signal_1700_1000->Draw("same");
	p_tchiwg_HT_signal_800->SetLineColor(28);
	p_tchiwg_HT_signal_800->SetLineStyle(2);
	p_tchiwg_HT_signal_800->SetLineWidth(4);
	p_tchiwg_HT_signal_800->Draw("same");

	c_HT->cd();
	TPad *HT_pad2 = new TPad("HT_pad2", "HT_pad2", 0, 0, 1, 0.3);
	HT_pad2->SetBottomMargin(0.4);
	HT_pad2->Draw();
	HT_pad2->cd();
  TLine *flatratio_HT = new TLine(0,1,2000,1);
	TH1F *ratio_HT=(TH1F*)p_allHT->Clone("transfer factor");
	ratio_HT->SetMarkerStyle(20);
	ratio_HT->SetLineColor(kBlack);
	ratio_HT->GetXaxis()->SetRangeUser(0,2000);
	ratio_HT->GetYaxis()->SetRangeUser(0,2);
	ratio_HT->GetYaxis()->SetNdivisions(504);
	ratio_HT->SetMinimum(0);
	ratio_HT->SetMaximum(2);
	ratio_HT->Divide(p_VGHT);
	ratio_HT->SetTitle("");
	ratio_HT->GetYaxis()->SetTitle("#frac{Obs.}{Bkg.}");
	ratio_HT->Draw();
	ratioerror_HT->SetFillColor(12);
	ratioerror_HT->SetFillStyle(3345);
	ratioerror_HT->Draw("E2 same");
	flatratio_HT->Draw("same");
	c_HT->SaveAs(htplot.str().c_str());
}


