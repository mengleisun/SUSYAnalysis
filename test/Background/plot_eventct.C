#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TF3.h"
#include "TH1D.h"
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
void plot_eventct(){//main  

  gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");
	
	TFile *egfile_ele = TFile::Open("signalTree_egamma_eleBkg.root");
	TFile *egfile_jet = TFile::Open("signalTree_egamma_jetbkg.root");
	TFile *egfile_qcd = TFile::Open("signalTree_egamma_qcd.root");
	TFile *egfile_VG  = TFile::Open("signalTree_egamma_VGBkg.root");
	TFile *egfile_rare= TFile::Open("signalTree_egamma_rareBkg.root");
	TFile *mgfile_ele = TFile::Open("signalTree_mg_eleBkg.root");
	TFile *mgfile_jet = TFile::Open("signalTree_mg_jetbkg.root");
	TFile *mgfile_qcd = TFile::Open("signalTree_mg_qcd.root");
	TFile *mgfile_VG  = TFile::Open("signalTree_mg_VGBkg.root");
	TFile *mgfile_rare= TFile::Open("signalTree_mg_rareBkg.root");

	TH1D *h_elefakepho_norm           = new TH1D("h_elefakepho_norm",            "h_elefakepho_norm",           36,0,36); 
	TH1D *h_elefakepho_syserr_jes     = new TH1D("h_elefakepho_syserr_jes",      "h_elefakepho_syserr_jes",     36,0,36); 
	TH1D *h_elefakepho_syserr_jer     = new TH1D("h_elefakepho_syserr_jer",      "h_elefakepho_syserr_jer",     36,0,36);   
	TH1D *h_elefakepho_syserr_esf     = new TH1D("h_elefakepho_syserr_esf",      "h_elefakepho_syserr_esf",     36,0,36); 
	TH1D *h_elefakepho_syserr_scale   = new TH1D("h_elefakepho_syserr_scale",    "h_elefakepho_syserr_scale",   36,0,36); 
	TH1D *h_elefakepho_syserr_eleshape= new TH1D("h_elefakepho_syserr_eleshape", "h_elefakepho_syserr_eleshape",36,0,36); 
	TH1D *h_elefakepho_syserr_jetshape= new TH1D("h_elefakepho_syserr_jetshape", "h_elefakepho_syserr_jetshape",36,0,36); 
	TH1D *h_elefakepho_syserr_xs      = new TH1D("h_elefakepho_syserr_xs",       "h_elefakepho_syserr_xs",      36,0,36);
	TH1D *h_elefakepho_syserr_lumi    = new TH1D("h_elefakepho_syserr_lumi",     "h_elefakepho_syserr_lumi",    36,0,36); 
	                                                                                                          
	TH1D *h_jetfakepho_norm           = new TH1D("h_jetfakepho_norm",            "h_jetfakepho_norm",           36,0,36); 
	TH1D *h_jetfakepho_syserr_jes     = new TH1D("h_jetfakepho_syserr_jes",      "h_jetfakepho_syserr_jes",     36,0,36); 
	TH1D *h_jetfakepho_syserr_jer     = new TH1D("h_jetfakepho_syserr_jer",      "h_jetfakepho_syserr_jer",     36,0,36); 
	TH1D *h_jetfakepho_syserr_esf     = new TH1D("h_jetfakepho_syserr_esf",      "h_jetfakepho_syserr_esf",     36,0,36); 
	TH1D *h_jetfakepho_syserr_scale   = new TH1D("h_jetfakepho_syserr_scale",    "h_jetfakepho_syserr_scale",   36,0,36); 
	TH1D *h_jetfakepho_syserr_eleshape= new TH1D("h_jetfakepho_syserr_eleshape", "h_jetfakepho_syserr_eleshape",36,0,36); 
	TH1D *h_jetfakepho_syserr_jetshape= new TH1D("h_jetfakepho_syserr_jetshape", "h_jetfakepho_syserr_jetshape",36,0,36); 
	TH1D *h_jetfakepho_syserr_xs      = new TH1D("h_jetfakepho_syserr_xs",       "h_jetfakepho_syserr_xs",      36,0,36);
	TH1D *h_jetfakepho_syserr_lumi    = new TH1D("h_jetfakepho_syserr_lumi",     "h_jetfakepho_syserr_lumi",    36,0,36); 
	                                                                             
	TH1D *h_qcdfakelep_norm           = new TH1D("h_qcdfakelep_norm",            "h_qcdfakelep_norm",           36,0,36); 
	TH1D *h_qcdfakelep_syserr_jes     = new TH1D("h_qcdfakelep_syserr_jes",      "h_qcdfakelep_syserr_jes",     36,0,36); 
	TH1D *h_qcdfakelep_syserr_jer     = new TH1D("h_qcdfakelep_syserr_jer",      "h_qcdfakelep_syserr_jer",     36,0,36); 
	TH1D *h_qcdfakelep_syserr_esf     = new TH1D("h_qcdfakelep_syserr_esf",      "h_qcdfakelep_syserr_esf",     36,0,36); 
	TH1D *h_qcdfakelep_syserr_scale   = new TH1D("h_qcdfakelep_syserr_scale",    "h_qcdfakelep_syserr_scale",   36,0,36); 
	TH1D *h_qcdfakelep_syserr_eleshape= new TH1D("h_qcdfakelep_syserr_eleshape", "h_qcdfakelep_syserr_eleshape",36,0,36); 
	TH1D *h_qcdfakelep_syserr_jetshape= new TH1D("h_qcdfakelep_syserr_jetshape", "h_qcdfakelep_syserr_jetshape",36,0,36); 
	TH1D *h_qcdfakelep_syserr_xs      = new TH1D("h_qcdfakelep_syserr_xs",       "h_qcdfakelep_syserr_xs",      36,0,36);
	TH1D *h_qcdfakelep_syserr_lumi    = new TH1D("h_qcdfakelep_syserr_lumi",     "h_qcdfakelep_syserr_lumi",    36,0,36); 
	
	TH1D *h_VGamma_norm           = new TH1D("h_VGamma_norm",             "h_VGamma_norm",           36,0,36);   
	TH1D *h_VGamma_syserr_jes     = new TH1D("h_VGamma_syserr_jes",       "h_VGamma_syserr_jes",     36,0,36);   
	TH1D *h_VGamma_syserr_jer     = new TH1D("h_VGamma_syserr_jer",       "h_VGamma_syserr_jer",     36,0,36);   
	TH1D *h_VGamma_syserr_esf     = new TH1D("h_VGamma_syserr_esf",       "h_VGamma_syserr_esf",     36,0,36);   
	TH1D *h_VGamma_syserr_scale   = new TH1D("h_VGamma_syserr_scale",     "h_VGamma_syserr_scale",   36,0,36);   
	TH1D *h_VGamma_syserr_eleshape= new TH1D("h_VGamma_syserr_eleshape",  "h_VGamma_syserr_eleshape",36,0,36);   
	TH1D *h_VGamma_syserr_jetshape= new TH1D("h_VGamma_syserr_jetshape",  "h_VGamma_syserr_jetshape",36,0,36);   
	TH1D *h_VGamma_syserr_xs      = new TH1D("h_VGamma_syserr_xs",        "h_VGamma_syserr_xs",      36,0,36);   
	TH1D *h_VGamma_syserr_lumi    = new TH1D("h_VGamma_syserr_lumi",      "h_VGamma_syserr_lumi",    36,0,36);   
	
	TH1D *h_rare_norm           = new TH1D("h_rare_norm",              "h_rare_norm",           36,0,36);   
	TH1D *h_rare_syserr_jes     = new TH1D("h_rare_syserr_jes",        "h_rare_syserr_jes",     36,0,36);   
	TH1D *h_rare_syserr_jer     = new TH1D("h_rare_syserr_jer",        "h_rare_syserr_jer",     36,0,36);   
	TH1D *h_rare_syserr_esf     = new TH1D("h_rare_syserr_esf",        "h_rare_syserr_esf",     36,0,36);   
	TH1D *h_rare_syserr_scale   = new TH1D("h_rare_syserr_scale",      "h_rare_syserr_scale",   36,0,36);   
	TH1D *h_rare_syserr_eleshape= new TH1D("h_rare_syserr_eleshape",   "h_rare_syserr_eleshape",36,0,36);   
	TH1D *h_rare_syserr_jetshape= new TH1D("h_rare_syserr_jetshape",   "h_rare_syserr_jetshape",36,0,36);   
	TH1D *h_rare_syserr_xs      = new TH1D("h_rare_syserr_xs",         "h_rare_syserr_xs",      36,0,36);   
	TH1D *h_rare_syserr_lumi    = new TH1D("h_rare_syserr_lumi",       "h_rare_syserr_lumi",    36,0,36);   



	TH1D *eg_elefakepho_norm           = (TH1D*)egfile_ele->Get("eg_elefakepho_norm");         
	TH1D *eg_elefakepho_syserr_jes     = (TH1D*)egfile_ele->Get("eg_elefakepho_syserr_jes");     
	TH1D *eg_elefakepho_syserr_jer     = (TH1D*)egfile_ele->Get("eg_elefakepho_syserr_jer");  
	TH1D *eg_elefakepho_syserr_esf     = (TH1D*)egfile_ele->Get("eg_elefakepho_syserr_esf");  
	TH1D *eg_elefakepho_syserr_scale   = (TH1D*)egfile_ele->Get("eg_elefakepho_syserr_scale");  
	TH1D *eg_elefakepho_syserr_eleshape= (TH1D*)egfile_ele->Get("eg_elefakepho_syserr_eleshape");
	TH1D *eg_elefakepho_syserr_jetshape= (TH1D*)egfile_ele->Get("eg_elefakepho_syserr_jetshape");
	TH1D *eg_elefakepho_syserr_xs      = (TH1D*)egfile_ele->Get("eg_elefakepho_syserr_xs");
	TH1D *eg_elefakepho_syserr_lumi    = (TH1D*)egfile_ele->Get("eg_elefakepho_syserr_lumi");    
	
	TH1D *eg_jetfakepho_norm           = (TH1D*)egfile_jet->Get("eg_jetfakepho_norm");         
	TH1D *eg_jetfakepho_syserr_jes     = (TH1D*)egfile_jet->Get("eg_jetfakepho_syserr_jes");     
	TH1D *eg_jetfakepho_syserr_jer     = (TH1D*)egfile_jet->Get("eg_jetfakepho_syserr_jer");  
	TH1D *eg_jetfakepho_syserr_esf     = (TH1D*)egfile_jet->Get("eg_jetfakepho_syserr_esf");  
	TH1D *eg_jetfakepho_syserr_scale   = (TH1D*)egfile_jet->Get("eg_jetfakepho_syserr_scale");  
	TH1D *eg_jetfakepho_syserr_eleshape= (TH1D*)egfile_jet->Get("eg_jetfakepho_syserr_eleshape");
	TH1D *eg_jetfakepho_syserr_jetshape= (TH1D*)egfile_jet->Get("eg_jetfakepho_syserr_jetshape");
	TH1D *eg_jetfakepho_syserr_xs      = (TH1D*)egfile_jet->Get("eg_jetfakepho_syserr_xs");
	TH1D *eg_jetfakepho_syserr_lumi    = (TH1D*)egfile_jet->Get("eg_jetfakepho_syserr_lumi");    
	
	TH1D *eg_qcdfakelep_norm           = (TH1D*)egfile_qcd->Get("eg_qcdfakelep_norm");         
	TH1D *eg_qcdfakelep_syserr_jes     = (TH1D*)egfile_qcd->Get("eg_qcdfakelep_syserr_jes");     
	TH1D *eg_qcdfakelep_syserr_jer     = (TH1D*)egfile_qcd->Get("eg_qcdfakelep_syserr_jer");  
	TH1D *eg_qcdfakelep_syserr_esf     = (TH1D*)egfile_qcd->Get("eg_qcdfakelep_syserr_esf");  
	TH1D *eg_qcdfakelep_syserr_scale   = (TH1D*)egfile_qcd->Get("eg_qcdfakelep_syserr_scale");  
	TH1D *eg_qcdfakelep_syserr_eleshape= (TH1D*)egfile_qcd->Get("eg_qcdfakelep_syserr_eleshape");
	TH1D *eg_qcdfakelep_syserr_jetshape= (TH1D*)egfile_qcd->Get("eg_qcdfakelep_syserr_jetshape");
	TH1D *eg_qcdfakelep_syserr_xs      = (TH1D*)egfile_qcd->Get("eg_qcdfakelep_syserr_xs");
	TH1D *eg_qcdfakelep_syserr_lumi    = (TH1D*)egfile_qcd->Get("eg_qcdfakelep_syserr_lumi");    
	
	TH1D *eg_VGamma_norm           = (TH1D*)egfile_VG->Get("eg_VGamma_norm");         
	TH1D *eg_VGamma_syserr_jes     = (TH1D*)egfile_VG->Get("eg_VGamma_syserr_jes");     
	TH1D *eg_VGamma_syserr_jer     = (TH1D*)egfile_VG->Get("eg_VGamma_syserr_jer");  
	TH1D *eg_VGamma_syserr_esf     = (TH1D*)egfile_VG->Get("eg_VGamma_syserr_esf");  
	TH1D *eg_VGamma_syserr_scale   = (TH1D*)egfile_VG->Get("eg_VGamma_syserr_scale");  
	TH1D *eg_VGamma_syserr_eleshape= (TH1D*)egfile_VG->Get("eg_VGamma_syserr_eleshape");
	TH1D *eg_VGamma_syserr_jetshape= (TH1D*)egfile_VG->Get("eg_VGamma_syserr_jetshape");
	TH1D *eg_VGamma_syserr_xs      = (TH1D*)egfile_VG->Get("eg_VGamma_syserr_xs");
	TH1D *eg_VGamma_syserr_lumi    = (TH1D*)egfile_VG->Get("eg_VGamma_syserr_lumi");    
	
	TH1D *eg_rare_norm           = (TH1D*)egfile_rare->Get("eg_rare_norm");         
	TH1D *eg_rare_syserr_jes     = (TH1D*)egfile_rare->Get("eg_rare_syserr_jes");     
	TH1D *eg_rare_syserr_jer     = (TH1D*)egfile_rare->Get("eg_rare_syserr_jer");  
	TH1D *eg_rare_syserr_esf     = (TH1D*)egfile_rare->Get("eg_rare_syserr_esf");  
	TH1D *eg_rare_syserr_scale   = (TH1D*)egfile_rare->Get("eg_rare_syserr_scale");  
	TH1D *eg_rare_syserr_eleshape= (TH1D*)egfile_rare->Get("eg_rare_syserr_eleshape");
	TH1D *eg_rare_syserr_jetshape= (TH1D*)egfile_rare->Get("eg_rare_syserr_jetshape");
	TH1D *eg_rare_syserr_xs      = (TH1D*)egfile_rare->Get("eg_rare_syserr_xs");
	TH1D *eg_rare_syserr_lumi    = (TH1D*)egfile_rare->Get("eg_rare_syserr_lumi");    

	TH1D *mg_elefakepho_norm           = (TH1D*)mgfile_ele->Get("mg_elefakepho_norm");         
	TH1D *mg_elefakepho_syserr_jes     = (TH1D*)mgfile_ele->Get("mg_elefakepho_syserr_jes");     
	TH1D *mg_elefakepho_syserr_jer     = (TH1D*)mgfile_ele->Get("mg_elefakepho_syserr_jer");  
	TH1D *mg_elefakepho_syserr_esf     = (TH1D*)mgfile_ele->Get("mg_elefakepho_syserr_esf");  
	TH1D *mg_elefakepho_syserr_scale   = (TH1D*)mgfile_ele->Get("mg_elefakepho_syserr_scale");  
	TH1D *mg_elefakepho_syserr_eleshape= (TH1D*)mgfile_ele->Get("mg_elefakepho_syserr_eleshape");
	TH1D *mg_elefakepho_syserr_jetshape= (TH1D*)mgfile_ele->Get("mg_elefakepho_syserr_jetshape");
	TH1D *mg_elefakepho_syserr_xs      = (TH1D*)mgfile_ele->Get("mg_elefakepho_syserr_xs");
	TH1D *mg_elefakepho_syserr_lumi    = (TH1D*)mgfile_ele->Get("mg_elefakepho_syserr_lumi");    
	
	TH1D *mg_jetfakepho_norm           = (TH1D*)mgfile_jet->Get("mg_jetfakepho_norm");         
	TH1D *mg_jetfakepho_syserr_jes     = (TH1D*)mgfile_jet->Get("mg_jetfakepho_syserr_jes");     
	TH1D *mg_jetfakepho_syserr_jer     = (TH1D*)mgfile_jet->Get("mg_jetfakepho_syserr_jer");  
	TH1D *mg_jetfakepho_syserr_esf     = (TH1D*)mgfile_jet->Get("mg_jetfakepho_syserr_esf");  
	TH1D *mg_jetfakepho_syserr_scale   = (TH1D*)mgfile_jet->Get("mg_jetfakepho_syserr_scale");  
	TH1D *mg_jetfakepho_syserr_eleshape= (TH1D*)mgfile_jet->Get("mg_jetfakepho_syserr_eleshape");
	TH1D *mg_jetfakepho_syserr_jetshape= (TH1D*)mgfile_jet->Get("mg_jetfakepho_syserr_jetshape");
	TH1D *mg_jetfakepho_syserr_xs      = (TH1D*)mgfile_jet->Get("mg_jetfakepho_syserr_xs");
	TH1D *mg_jetfakepho_syserr_lumi    = (TH1D*)mgfile_jet->Get("mg_jetfakepho_syserr_lumi");    
	
	TH1D *mg_qcdfakelep_norm           = (TH1D*)mgfile_qcd->Get("mg_qcdfakelep_norm");         
	TH1D *mg_qcdfakelep_syserr_jes     = (TH1D*)mgfile_qcd->Get("mg_qcdfakelep_syserr_jes");     
	TH1D *mg_qcdfakelep_syserr_jer     = (TH1D*)mgfile_qcd->Get("mg_qcdfakelep_syserr_jer");  
	TH1D *mg_qcdfakelep_syserr_esf     = (TH1D*)mgfile_qcd->Get("mg_qcdfakelep_syserr_esf");  
	TH1D *mg_qcdfakelep_syserr_scale   = (TH1D*)mgfile_qcd->Get("mg_qcdfakelep_syserr_scale");  
	TH1D *mg_qcdfakelep_syserr_eleshape= (TH1D*)mgfile_qcd->Get("mg_qcdfakelep_syserr_eleshape");
	TH1D *mg_qcdfakelep_syserr_jetshape= (TH1D*)mgfile_qcd->Get("mg_qcdfakelep_syserr_jetshape");
	TH1D *mg_qcdfakelep_syserr_xs      = (TH1D*)mgfile_qcd->Get("mg_qcdfakelep_syserr_xs");
	TH1D *mg_qcdfakelep_syserr_lumi    = (TH1D*)mgfile_qcd->Get("mg_qcdfakelep_syserr_lumi");    
	
	TH1D *mg_VGamma_norm           = (TH1D*)mgfile_VG->Get("mg_VGamma_norm");         
	TH1D *mg_VGamma_syserr_jes     = (TH1D*)mgfile_VG->Get("mg_VGamma_syserr_jes");     
	TH1D *mg_VGamma_syserr_jer     = (TH1D*)mgfile_VG->Get("mg_VGamma_syserr_jer");  
	TH1D *mg_VGamma_syserr_esf     = (TH1D*)mgfile_VG->Get("mg_VGamma_syserr_esf");  
	TH1D *mg_VGamma_syserr_scale   = (TH1D*)mgfile_VG->Get("mg_VGamma_syserr_scale");  
	TH1D *mg_VGamma_syserr_eleshape= (TH1D*)mgfile_VG->Get("mg_VGamma_syserr_eleshape");
	TH1D *mg_VGamma_syserr_jetshape= (TH1D*)mgfile_VG->Get("mg_VGamma_syserr_jetshape");
	TH1D *mg_VGamma_syserr_xs      = (TH1D*)mgfile_VG->Get("mg_VGamma_syserr_xs");
	TH1D *mg_VGamma_syserr_lumi    = (TH1D*)mgfile_VG->Get("mg_VGamma_syserr_lumi");    
	
	TH1D *mg_rare_norm           = (TH1D*)mgfile_rare->Get("mg_rare_norm");         
	for(unsigned ibin(1); ibin <= 18; ibin++){
		std::cout << "rare_norm_mg " << ibin << " " << mg_rare_norm->GetBinContent(ibin) << std::endl;
	}
	TH1D *mg_rare_syserr_jes     = (TH1D*)mgfile_rare->Get("mg_rare_syserr_jes");     
	TH1D *mg_rare_syserr_jer     = (TH1D*)mgfile_rare->Get("mg_rare_syserr_jer");  
	TH1D *mg_rare_syserr_esf     = (TH1D*)mgfile_rare->Get("mg_rare_syserr_esf");  
	TH1D *mg_rare_syserr_scale   = (TH1D*)mgfile_rare->Get("mg_rare_syserr_scale");  
	TH1D *mg_rare_syserr_eleshape= (TH1D*)mgfile_rare->Get("mg_rare_syserr_eleshape");
	TH1D *mg_rare_syserr_jetshape= (TH1D*)mgfile_rare->Get("mg_rare_syserr_jetshape");
	TH1D *mg_rare_syserr_xs      = (TH1D*)mgfile_rare->Get("mg_rare_syserr_xs");
	TH1D *mg_rare_syserr_lumi    = (TH1D*)mgfile_rare->Get("mg_rare_syserr_lumi");    


	for(unsigned ibin(1); ibin <= 18; ibin++){
		h_elefakepho_norm->SetBinContent(ibin, mg_elefakepho_norm->GetBinContent(ibin));
		h_elefakepho_syserr_jes->SetBinContent(ibin,      mg_elefakepho_syserr_jes->GetBinContent(ibin));      
		h_elefakepho_syserr_jer->SetBinContent(ibin,      mg_elefakepho_syserr_jer->GetBinContent(ibin)); 
		h_elefakepho_syserr_esf->SetBinContent(ibin,      mg_elefakepho_syserr_esf->GetBinContent(ibin));      
		h_elefakepho_syserr_scale->SetBinContent(ibin,    mg_elefakepho_syserr_scale->GetBinContent(ibin));    
		h_elefakepho_syserr_eleshape->SetBinContent(ibin, mg_elefakepho_syserr_eleshape->GetBinContent(ibin));
		h_elefakepho_syserr_jetshape->SetBinContent(ibin, mg_elefakepho_syserr_jetshape->GetBinContent(ibin));
		h_elefakepho_syserr_xs->SetBinContent(ibin,       mg_elefakepho_syserr_xs->GetBinContent(ibin));       
		h_elefakepho_syserr_lumi->SetBinContent(ibin,     mg_elefakepho_syserr_lumi->GetBinContent(ibin));     
																																																																																																							
		h_jetfakepho_norm->SetBinContent(ibin,            mg_jetfakepho_norm->GetBinContent(ibin));            
		h_jetfakepho_syserr_jes->SetBinContent(ibin,      mg_jetfakepho_syserr_jes->GetBinContent(ibin));      
		h_jetfakepho_syserr_jer->SetBinContent(ibin,      mg_jetfakepho_syserr_jer->GetBinContent(ibin));      
		h_jetfakepho_syserr_esf->SetBinContent(ibin,      mg_jetfakepho_syserr_esf->GetBinContent(ibin));      
		h_jetfakepho_syserr_scale->SetBinContent(ibin,    mg_jetfakepho_syserr_scale->GetBinContent(ibin));    
		h_jetfakepho_syserr_eleshape->SetBinContent(ibin, mg_jetfakepho_syserr_eleshape->GetBinContent(ibin));
		h_jetfakepho_syserr_jetshape->SetBinContent(ibin, mg_jetfakepho_syserr_jetshape->GetBinContent(ibin));
		h_jetfakepho_syserr_xs->SetBinContent(ibin,       mg_jetfakepho_syserr_xs->GetBinContent(ibin));       
		h_jetfakepho_syserr_lumi->SetBinContent(ibin,     mg_jetfakepho_syserr_lumi->GetBinContent(ibin));     
																														
		h_qcdfakelep_norm->SetBinContent(ibin,            mg_qcdfakelep_norm->GetBinContent(ibin));            
		h_qcdfakelep_syserr_jes->SetBinContent(ibin,      mg_qcdfakelep_syserr_jes->GetBinContent(ibin));      
		h_qcdfakelep_syserr_jer->SetBinContent(ibin,      mg_qcdfakelep_syserr_jer->GetBinContent(ibin));      
		h_qcdfakelep_syserr_esf->SetBinContent(ibin,      mg_qcdfakelep_syserr_esf->GetBinContent(ibin));      
		h_qcdfakelep_syserr_scale->SetBinContent(ibin,    mg_qcdfakelep_syserr_scale->GetBinContent(ibin));    
		h_qcdfakelep_syserr_eleshape->SetBinContent(ibin, mg_qcdfakelep_syserr_eleshape->GetBinContent(ibin));
		h_qcdfakelep_syserr_jetshape->SetBinContent(ibin, mg_qcdfakelep_syserr_jetshape->GetBinContent(ibin));
		h_qcdfakelep_syserr_xs->SetBinContent(ibin,       mg_qcdfakelep_syserr_xs->GetBinContent(ibin));       
		h_qcdfakelep_syserr_lumi->SetBinContent(ibin,     mg_qcdfakelep_syserr_lumi->GetBinContent(ibin));     
																
		h_VGamma_norm->SetBinContent(ibin,                mg_VGamma_norm->GetBinContent(ibin));            
		h_VGamma_syserr_jes->SetBinContent(ibin,          mg_VGamma_syserr_jes->GetBinContent(ibin));      
		h_VGamma_syserr_jer->SetBinContent(ibin,          mg_VGamma_syserr_jer->GetBinContent(ibin));      
		h_VGamma_syserr_esf->SetBinContent(ibin,          mg_VGamma_syserr_esf->GetBinContent(ibin));      
		h_VGamma_syserr_scale->SetBinContent(ibin,        mg_VGamma_syserr_scale->GetBinContent(ibin));    
		h_VGamma_syserr_eleshape->SetBinContent(ibin,     mg_VGamma_syserr_eleshape->GetBinContent(ibin));
		h_VGamma_syserr_jetshape->SetBinContent(ibin,     mg_VGamma_syserr_jetshape->GetBinContent(ibin));
		h_VGamma_syserr_xs->SetBinContent(ibin,           mg_VGamma_syserr_xs->GetBinContent(ibin));       
		h_VGamma_syserr_lumi->SetBinContent(ibin,         mg_VGamma_syserr_lumi->GetBinContent(ibin));     
																
		h_rare_norm->SetBinContent(ibin,                  mg_rare_norm->GetBinContent(ibin));            
		h_rare_syserr_jes->SetBinContent(ibin,            mg_rare_syserr_jes->GetBinContent(ibin));      
		h_rare_syserr_jer->SetBinContent(ibin,            mg_rare_syserr_jer->GetBinContent(ibin));      
		h_rare_syserr_esf->SetBinContent(ibin,            mg_rare_syserr_esf->GetBinContent(ibin));      
		h_rare_syserr_scale->SetBinContent(ibin,          mg_rare_syserr_scale->GetBinContent(ibin));   
		h_rare_syserr_eleshape->SetBinContent(ibin,       mg_rare_syserr_eleshape->GetBinContent(ibin));   
		h_rare_syserr_jetshape->SetBinContent(ibin,       mg_rare_syserr_jetshape->GetBinContent(ibin)); 
		h_rare_syserr_xs->SetBinContent(ibin,             mg_rare_syserr_xs->GetBinContent(ibin));       
		h_rare_syserr_lumi->SetBinContent(ibin,           mg_rare_syserr_lumi->GetBinContent(ibin));     


		h_elefakepho_norm->SetBinContent(ibin+18, eg_elefakepho_norm->GetBinContent(ibin));
		h_elefakepho_syserr_jes->SetBinContent(ibin+18,      eg_elefakepho_syserr_jes->GetBinContent(ibin));      
		h_elefakepho_syserr_jer->SetBinContent(ibin+18,      eg_elefakepho_syserr_jer->GetBinContent(ibin)); 
		h_elefakepho_syserr_esf->SetBinContent(ibin+18,      eg_elefakepho_syserr_esf->GetBinContent(ibin));      
		h_elefakepho_syserr_scale->SetBinContent(ibin+18,    eg_elefakepho_syserr_scale->GetBinContent(ibin));    
		h_elefakepho_syserr_eleshape->SetBinContent(ibin+18, eg_elefakepho_syserr_eleshape->GetBinContent(ibin));
		h_elefakepho_syserr_jetshape->SetBinContent(ibin+18, eg_elefakepho_syserr_jetshape->GetBinContent(ibin));
		h_elefakepho_syserr_xs->SetBinContent(ibin+18,       eg_elefakepho_syserr_xs->GetBinContent(ibin));       
		h_elefakepho_syserr_lumi->SetBinContent(ibin+18,     eg_elefakepho_syserr_lumi->GetBinContent(ibin));     
																																																																																																							
		h_jetfakepho_norm->SetBinContent(ibin+18,            eg_jetfakepho_norm->GetBinContent(ibin));            
		h_jetfakepho_syserr_jes->SetBinContent(ibin+18,      eg_jetfakepho_syserr_jes->GetBinContent(ibin));      
		h_jetfakepho_syserr_jer->SetBinContent(ibin+18,      eg_jetfakepho_syserr_jer->GetBinContent(ibin));      
		h_jetfakepho_syserr_esf->SetBinContent(ibin+18,      eg_jetfakepho_syserr_esf->GetBinContent(ibin));      
		h_jetfakepho_syserr_scale->SetBinContent(ibin+18,    eg_jetfakepho_syserr_scale->GetBinContent(ibin));    
		h_jetfakepho_syserr_eleshape->SetBinContent(ibin+18, eg_jetfakepho_syserr_eleshape->GetBinContent(ibin));
		h_jetfakepho_syserr_jetshape->SetBinContent(ibin+18, eg_jetfakepho_syserr_jetshape->GetBinContent(ibin));
		h_jetfakepho_syserr_xs->SetBinContent(ibin+18,       eg_jetfakepho_syserr_xs->GetBinContent(ibin));       
		h_jetfakepho_syserr_lumi->SetBinContent(ibin+18,     eg_jetfakepho_syserr_lumi->GetBinContent(ibin));     
																														
		h_qcdfakelep_norm->SetBinContent(ibin+18,            eg_qcdfakelep_norm->GetBinContent(ibin));            
		h_qcdfakelep_syserr_jes->SetBinContent(ibin+18,      eg_qcdfakelep_syserr_jes->GetBinContent(ibin));      
		h_qcdfakelep_syserr_jer->SetBinContent(ibin+18,      eg_qcdfakelep_syserr_jer->GetBinContent(ibin));      
		h_qcdfakelep_syserr_esf->SetBinContent(ibin+18,      eg_qcdfakelep_syserr_esf->GetBinContent(ibin));      
		h_qcdfakelep_syserr_scale->SetBinContent(ibin+18,    eg_qcdfakelep_syserr_scale->GetBinContent(ibin));    
		h_qcdfakelep_syserr_eleshape->SetBinContent(ibin+18, eg_qcdfakelep_syserr_eleshape->GetBinContent(ibin));
		h_qcdfakelep_syserr_jetshape->SetBinContent(ibin+18, eg_qcdfakelep_syserr_jetshape->GetBinContent(ibin));
		h_qcdfakelep_syserr_xs->SetBinContent(ibin+18,       eg_qcdfakelep_syserr_xs->GetBinContent(ibin));       
		h_qcdfakelep_syserr_lumi->SetBinContent(ibin+18,     eg_qcdfakelep_syserr_lumi->GetBinContent(ibin));     
																
		h_VGamma_norm->SetBinContent(ibin+18,                eg_VGamma_norm->GetBinContent(ibin));            
		h_VGamma_syserr_jes->SetBinContent(ibin+18,          eg_VGamma_syserr_jes->GetBinContent(ibin));      
		h_VGamma_syserr_jer->SetBinContent(ibin+18,          eg_VGamma_syserr_jer->GetBinContent(ibin));      
		h_VGamma_syserr_esf->SetBinContent(ibin+18,          eg_VGamma_syserr_esf->GetBinContent(ibin));      
		h_VGamma_syserr_scale->SetBinContent(ibin+18,        eg_VGamma_syserr_scale->GetBinContent(ibin));    
		h_VGamma_syserr_eleshape->SetBinContent(ibin+18,     eg_VGamma_syserr_eleshape->GetBinContent(ibin));
		h_VGamma_syserr_jetshape->SetBinContent(ibin+18,     eg_VGamma_syserr_jetshape->GetBinContent(ibin));
		h_VGamma_syserr_xs->SetBinContent(ibin+18,           eg_VGamma_syserr_xs->GetBinContent(ibin));       
		h_VGamma_syserr_lumi->SetBinContent(ibin+18,         eg_VGamma_syserr_lumi->GetBinContent(ibin));     
																
		h_rare_norm->SetBinContent(ibin+18,                  eg_rare_norm->GetBinContent(ibin));            
		h_rare_syserr_jes->SetBinContent(ibin+18,            eg_rare_syserr_jes->GetBinContent(ibin));      
		h_rare_syserr_jer->SetBinContent(ibin+18,            eg_rare_syserr_jer->GetBinContent(ibin));      
		h_rare_syserr_esf->SetBinContent(ibin+18,            eg_rare_syserr_esf->GetBinContent(ibin));      
		h_rare_syserr_scale->SetBinContent(ibin+18,          eg_rare_syserr_scale->GetBinContent(ibin));   
		h_rare_syserr_eleshape->SetBinContent(ibin+18,       eg_rare_syserr_eleshape->GetBinContent(ibin));   
		h_rare_syserr_jetshape->SetBinContent(ibin+18,       eg_rare_syserr_jetshape->GetBinContent(ibin)); 
		h_rare_syserr_xs->SetBinContent(ibin+18,             eg_rare_syserr_xs->GetBinContent(ibin));       
		h_rare_syserr_lumi->SetBinContent(ibin+18,           eg_rare_syserr_lumi->GetBinContent(ibin));     


		h_elefakepho_norm->SetBinError(ibin,            mg_elefakepho_norm->GetBinError(ibin));
		h_jetfakepho_norm->SetBinError(ibin,            mg_jetfakepho_norm->GetBinError(ibin));            
		h_qcdfakelep_norm->SetBinError(ibin,            mg_qcdfakelep_norm->GetBinError(ibin));            
		h_VGamma_norm->SetBinError(ibin,                mg_VGamma_norm->GetBinError(ibin));            
		h_rare_norm->SetBinError(ibin,                  mg_rare_norm->GetBinError(ibin));            
		h_elefakepho_norm->SetBinError(ibin+18,          eg_elefakepho_norm->GetBinError(ibin));
		h_jetfakepho_norm->SetBinError(ibin+18,          eg_jetfakepho_norm->GetBinError(ibin));            
		h_qcdfakelep_norm->SetBinError(ibin+18,          eg_qcdfakelep_norm->GetBinError(ibin));            
		h_VGamma_norm->SetBinError(ibin+18,              eg_VGamma_norm->GetBinError(ibin));            
		h_rare_norm->SetBinError(ibin+18,                eg_rare_norm->GetBinError(ibin));            
	}

	TFile *outputfile = TFile::Open("SignalSystematic_egmg.root","RECREATE");
	outputfile->cd();
	h_elefakepho_norm->Write();       
	h_elefakepho_syserr_jes->Write();       
	h_elefakepho_syserr_jer->Write();       
	h_elefakepho_syserr_esf->Write();       
	h_elefakepho_syserr_scale->Write();     
	h_elefakepho_syserr_eleshape->Write();  
	h_elefakepho_syserr_jetshape->Write();  
	h_elefakepho_syserr_xs->Write();        
	h_elefakepho_syserr_lumi->Write();      
	h_jetfakepho_norm->Write();             
	h_jetfakepho_syserr_jes->Write();       
	h_jetfakepho_syserr_jer->Write();       
	h_jetfakepho_syserr_esf->Write();       
	h_jetfakepho_syserr_scale->Write();     
	h_jetfakepho_syserr_eleshape->Write();  
	h_jetfakepho_syserr_jetshape->Write();  
	h_jetfakepho_syserr_xs->Write();        
	h_jetfakepho_syserr_lumi->Write();      
	h_qcdfakelep_norm->Write();             
	h_qcdfakelep_syserr_jes->Write();       
	h_qcdfakelep_syserr_jer->Write();       
	h_qcdfakelep_syserr_esf->Write();       
	h_qcdfakelep_syserr_scale->Write();     
	h_qcdfakelep_syserr_eleshape->Write();  
	h_qcdfakelep_syserr_jetshape->Write();  
	h_qcdfakelep_syserr_xs->Write();        
	h_qcdfakelep_syserr_lumi->Write();      
	h_VGamma_norm->Write();             
	h_VGamma_syserr_jes->Write();       
	h_VGamma_syserr_jer->Write();       
	h_VGamma_syserr_esf->Write();       
	h_VGamma_syserr_scale->Write();     
	h_VGamma_syserr_eleshape->Write();  
	h_VGamma_syserr_jetshape->Write();  
	h_VGamma_syserr_xs->Write();        
	h_VGamma_syserr_lumi->Write();      
	h_rare_norm->Write();             
	h_rare_syserr_jes->Write();       
	h_rare_syserr_jer->Write();       
	h_rare_syserr_esf->Write();       
	h_rare_syserr_scale->Write();     
	h_rare_syserr_eleshape->Write();  
	h_rare_syserr_jetshape->Write();  
	h_rare_syserr_xs->Write();        
	h_rare_syserr_lumi->Write();      
	outputfile->Close();

}
