#include "../../include/analysis_commoncode.h"
#include "TProfile2D.h"
#define NBIN 18

int DIM;
int NBIN_2DX;
int NBIN_2DY;
double XBIN_low;
double XBIN_high;
double YBIN_low;
double YBIN_high;

TFile *outputfile_susy;
TH2D *p_SUSYMass;
TH2D *p_SUSYselect;
TH2D *p_lowPU_susy_pass;
TH2D *p_lowPU_susy_all;
TH2D *p_highPU_susy_pass;
TH2D *p_highPU_susy_all;

TH2D *susy_wg_chan_rate_nom[NBIN*2]; 
TH2D *susy_wg_chan_rate_jesUp[NBIN*2]; 
TH2D *susy_wg_chan_rate_jesDown[NBIN*2];    
TH2D *susy_wg_chan_rate_jerUp[NBIN*2];    
TH2D *susy_wg_chan_rate_jerDown[NBIN*2];    
TH2D *susy_wg_chan_rate_xsUp[NBIN*2];       
TH2D *susy_wg_chan_rate_esfUp[NBIN*2];       
                                
TH2D *susy_wg_syserr_PU;
TH2D *susy_wg_chan_syserr_jes[NBIN*2];      
TH2D *susy_wg_chan_syserr_jer[NBIN*2];     
TH2D *susy_wg_chan_syserr_esf[NBIN*2];     
TH2D *susy_wg_chan_syserr_scale[NBIN*2];   
TH2D *susy_wg_chan_syserr_eleshape[NBIN*2]; 
TH2D *susy_wg_chan_syserr_jetshape[NBIN*2];
TH2D *susy_wg_chan_syserr_qcdshape[NBIN*2];
TH2D *susy_wg_chan_syserr_xs[NBIN*2];
TH2D *susy_wg_chan_syserr_lumi[NBIN*2];     
TH2D *susy_wg_chan_syserr_isr[NBIN*2];     

TH2D *susy_gg_chan_rate_nom[NBIN*2]; 
TH2D *susy_gg_chan_rate_jesUp[NBIN*2]; 
TH2D *susy_gg_chan_rate_jesDown[NBIN*2];    
TH2D *susy_gg_chan_rate_jerUp[NBIN*2];    
TH2D *susy_gg_chan_rate_jerDown[NBIN*2];    
TH2D *susy_gg_chan_rate_xsUp[NBIN*2];       
TH2D *susy_gg_chan_rate_esfUp[NBIN*2];       
                                
TH2D *susy_gg_syserr_PU;
TH2D *susy_gg_chan_syserr_jes[NBIN*2];      
TH2D *susy_gg_chan_syserr_jer[NBIN*2];     
TH2D *susy_gg_chan_syserr_esf[NBIN*2];     
TH2D *susy_gg_chan_syserr_scale[NBIN*2];   
TH2D *susy_gg_chan_syserr_eleshape[NBIN*2]; 
TH2D *susy_gg_chan_syserr_jetshape[NBIN*2];
TH2D *susy_gg_chan_syserr_qcdshape[NBIN*2];
TH2D *susy_gg_chan_syserr_xs[NBIN*2];
TH2D *susy_gg_chan_syserr_lumi[NBIN*2];     
TH2D *susy_gg_chan_syserr_isr[NBIN*2];     

TH1D *p_susy_MET_signal_mg;
TH1D *p_susy_HT_signal_mg;
TH1D *p_susy_PhoEt_signal_mg;
TH1D *p_susy_MET_signal_eg;
TH1D *p_susy_HT_signal_eg;
TH1D *p_susy_PhoEt_signal_eg;


void init_histo(TString modelname, TString filename, int dim, int nxbin, double xbinlow, double xbinhigh, int nybin, double ybinlow, double ybinhigh){

	DIM = dim;	
	NBIN_2DX = nxbin;
	NBIN_2DY = nybin;
	XBIN_low = xbinlow;
	XBIN_high = xbinhigh;
	YBIN_low = ybinlow;
	YBIN_high = ybinhigh;

	std::ostringstream outputfilename;
	outputfilename.str("");
	outputfilename << filename;
	outputfile_susy = TFile::Open(outputfilename.str().c_str(),"RECREATE");
	outputfile_susy->cd();

	susy_wg_syserr_PU 	 = new TH2D("susy_wg_syserr_PU","susy_wg_syserr_PU", NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
	susy_gg_syserr_PU 	 = new TH2D("susy_gg_syserr_PU","susy_gg_syserr_PU", NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
	p_SUSYMass         = new TH2D("SUSYMass","",NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
	p_SUSYselect       = new TH2D("p_SUSYselect","",NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high); 
	p_lowPU_susy_pass  = new TH2D("p_lowPU_susy_pass","", NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
	p_lowPU_susy_all   = new TH2D("p_lowPU_susy_all", "", NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
	p_highPU_susy_pass = new TH2D("p_highPU_susy_pass","",NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
	p_highPU_susy_all  = new TH2D("p_highPU_susy_all","", NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);

	std::ostringstream histname;
	for(unsigned i(0); i < NBIN*2; i++){
		histname.str("");
		histname << "wg_chan" << i+1 << "_rate_nom";
		susy_wg_chan_rate_nom[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
	}

	for(unsigned i(0); i < NBIN*2; i++){
		histname.str("");
		histname << "wg_chan" << i+1 << "_rate_jesUp";
		susy_wg_chan_rate_jesUp[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
		histname.str("");
		histname << "wg_chan" << i+1 << "_rate_jesDown";
		susy_wg_chan_rate_jesDown[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
		histname.str("");
		histname << "wg_chan" << i+1 << "_rate_jerUp";
		susy_wg_chan_rate_jerUp[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
		histname.str("");
		histname << "wg_chan" << i+1 << "_rate_jerDown";
		susy_wg_chan_rate_jerDown[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
		histname.str("");
		histname << "wg_chan" << i+1 << "_rate_xsUp";
		susy_wg_chan_rate_xsUp[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
		histname.str("");
		histname << "wg_chan" << i+1 << "_rate_esfUp";
		susy_wg_chan_rate_esfUp[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
															
		histname.str("");
		histname << "wg_chan" << i+1 << "_syserr_jes";
		susy_wg_chan_syserr_jes[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
		histname.str("");
		histname << "wg_chan" << i+1 << "_syserr_jer";
		susy_wg_chan_syserr_jer[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
		histname.str("");
		histname << "wg_chan" << i+1 << "_syserr_esf";
		susy_wg_chan_syserr_esf[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
		histname.str("");
		histname << "wg_chan" << i+1 << "_syserr_scale";
		susy_wg_chan_syserr_scale[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
		histname.str("");
		histname << "wg_chan" << i+1 << "_syserr_eleshape";
		susy_wg_chan_syserr_eleshape[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
		histname.str("");
		histname << "wg_chan" << i+1 << "_syserr_jetshape";
		susy_wg_chan_syserr_jetshape[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
		histname.str("");
		histname << "wg_chan" << i+1 << "_syserr_qcdshape";
		susy_wg_chan_syserr_qcdshape[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
		histname.str("");
		histname << "wg_chan" << i+1 << "_syserr_xs";
		susy_wg_chan_syserr_xs[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
		histname.str("");
		histname << "wg_chan" << i+1 << "_syserr_lumi";
		susy_wg_chan_syserr_lumi[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
		histname.str("");
		histname << "wg_chan" << i+1 << "_syserr_isr";
		susy_wg_chan_syserr_isr[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
	}

	for(unsigned i(0); i < NBIN*2; i++){
		histname.str("");
		histname << "gg_chan" << i+1 << "_rate_nom";
		susy_gg_chan_rate_nom[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
	}

	for(unsigned i(0); i < NBIN*2; i++){
		histname.str("");
		histname << "gg_chan" << i+1 << "_rate_jesUp";
		susy_gg_chan_rate_jesUp[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
		histname.str("");
		histname << "gg_chan" << i+1 << "_rate_jesDown";
		susy_gg_chan_rate_jesDown[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
		histname.str("");
		histname << "gg_chan" << i+1 << "_rate_jerUp";
		susy_gg_chan_rate_jerUp[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
		histname.str("");
		histname << "gg_chan" << i+1 << "_rate_jerDown";
		susy_gg_chan_rate_jerDown[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
		histname.str("");
		histname << "gg_chan" << i+1 << "_rate_xsUp";
		susy_gg_chan_rate_xsUp[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
		histname.str("");
		histname << "gg_chan" << i+1 << "_rate_esfUp";
		susy_gg_chan_rate_esfUp[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
															
		histname.str("");
		histname << "gg_chan" << i+1 << "_syserr_jes";
		susy_gg_chan_syserr_jes[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
		histname.str("");
		histname << "gg_chan" << i+1 << "_syserr_jer";
		susy_gg_chan_syserr_jer[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
		histname.str("");
		histname << "gg_chan" << i+1 << "_syserr_esf";
		susy_gg_chan_syserr_esf[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
		histname.str("");
		histname << "gg_chan" << i+1 << "_syserr_scale";
		susy_gg_chan_syserr_scale[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
		histname.str("");
		histname << "gg_chan" << i+1 << "_syserr_eleshape";
		susy_gg_chan_syserr_eleshape[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
		histname.str("");
		histname << "gg_chan" << i+1 << "_syserr_jetshape";
		susy_gg_chan_syserr_jetshape[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
		histname.str("");
		histname << "gg_chan" << i+1 << "_syserr_qcdshape";
		susy_gg_chan_syserr_qcdshape[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
		histname.str("");
		histname << "gg_chan" << i+1 << "_syserr_xs";
		susy_gg_chan_syserr_xs[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
		histname.str("");
		histname << "gg_chan" << i+1 << "_syserr_lumi";
		susy_gg_chan_syserr_lumi[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
		histname.str("");
		histname << "gg_chan" << i+1 << "_syserr_isr";
		susy_gg_chan_syserr_isr[i] = new TH2D(histname.str().c_str(),histname.str().c_str(),NBIN_2DX, XBIN_low, XBIN_high, NBIN_2DY, YBIN_low, YBIN_high);
	}

		
	p_susy_MET_signal_mg = new TH1D("p_susy_MET_signal_1700_1000_mg","", nSigMETBins, sigMETBins);
	p_susy_HT_signal_mg  = new TH1D("p_susy_HT_signal_1700_1000_mg","",  nSigHTBins, sigHTBins);
	p_susy_PhoEt_signal_mg = new TH1D("p_susy_PhoEt_signal_1700_1000_mg","",nSigEtBins, sigEtBins);
	p_susy_MET_signal_eg = new TH1D("p_susy_MET_signal_1700_1000_eg","", nSigMETBins, sigMETBins);
	p_susy_HT_signal_eg  = new TH1D("p_susy_HT_signal_1700_1000_eg","",  nSigHTBins, sigHTBins);
	p_susy_PhoEt_signal_eg = new TH1D("p_susy_PhoEt_signal_1700_1000_eg","",nSigEtBins, sigEtBins);

}

void write_histo(){		
	outputfile_susy->Write();
	outputfile_susy->Close();
}


