#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TChain.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "../../include/bk_tdrstyle.C"
int RunYear = 2017;

float getPUESF(int nvertex){
	
	float pileupweights[100];
	if(RunYear == 2016){
	pileupweights[0] = 0;
        pileupweights[1] = 0;
        pileupweights[2] = 1.04404;
        pileupweights[3] = 0.671719;
        pileupweights[4] = 0.551154;
        pileupweights[5] = 0.571996;
        pileupweights[6] = 0.489252;
        pileupweights[7] = 0.600793;
        pileupweights[8] = 0.648089;
        pileupweights[9] = 0.70624;
        pileupweights[10] = 0.78362;
        pileupweights[11] = 0.775689;
        pileupweights[12] = 0.822737;
        pileupweights[13] = 0.900937;
        pileupweights[14] = 0.910491;
        pileupweights[15] = 0.918496;
        pileupweights[16] = 0.946654;
        pileupweights[17] = 0.964116;
        pileupweights[18] = 0.966068;
        pileupweights[19] = 1.00622;
        pileupweights[20] = 1.03616;
        pileupweights[21] = 1.03452;
        pileupweights[22] = 1.07318;
        pileupweights[23] = 1.05786;
        pileupweights[24] = 1.07856;
        pileupweights[25] = 1.10927;
        pileupweights[26] = 1.13568;
        pileupweights[27] = 1.15534;
        pileupweights[28] = 1.26312;
        pileupweights[29] = 1.22121;
        pileupweights[30] = 1.3091;
        pileupweights[31] = 1.38049;
        pileupweights[32] = 1.40534;
        pileupweights[33] = 1.62435;
        pileupweights[34] = 1.71987;
        pileupweights[35] = 1.69665;
        pileupweights[36] = 1.79234;
        pileupweights[37] = 2.21024;
        pileupweights[38] = 1.8312;
        pileupweights[39] = 2.66367;
        pileupweights[40] = 2.94597;
        pileupweights[41] = 2.70223;
        pileupweights[42] = 3.94469;
        pileupweights[43] = 3.60809;
        pileupweights[44] = 3.68486;
        pileupweights[45] = 3.55326;
        pileupweights[46] = 6.55086;
        pileupweights[47] = 3.37779;
        pileupweights[48] = 11.6687;
        pileupweights[49] = 0;
        pileupweights[50] = 3.5825;
        pileupweights[51] = 9.51922;
        pileupweights[52] = 1.99597;
        pileupweights[53] = 0;
        pileupweights[54] = 0;
        pileupweights[55] = 0;
        pileupweights[56] = 0;
        pileupweights[57] = 0;
        pileupweights[58] = 0;
        pileupweights[59] = 0;
        pileupweights[60] = 0;
        pileupweights[61] = 0;
        pileupweights[62] = 0.307072;
        pileupweights[63] = 0;
        pileupweights[64] = 0;
        pileupweights[65] = 0;
        pileupweights[66] = 0;
        pileupweights[67] = 0;
        pileupweights[68] = 0;
        pileupweights[69] = 0;
        pileupweights[70] = 0;
        pileupweights[71] = 0;
        pileupweights[72] = 0;
        pileupweights[73] = 0;
        pileupweights[74] = 0;
        pileupweights[75] = 0;
        pileupweights[76] = 0;
        pileupweights[77] = 0;
        pileupweights[78] = 0;
        pileupweights[79] = 0;
        pileupweights[80] = 0;
        pileupweights[81] = 0;
        pileupweights[82] = 0;
        pileupweights[83] = 0;
        pileupweights[84] = 0;
        pileupweights[85] = 0;
        pileupweights[86] = 0;
        pileupweights[87] = 0;
        pileupweights[88] = 0;
        pileupweights[89] = 0;
        pileupweights[90] = 0;
        pileupweights[91] = 0;
        pileupweights[92] = 0;
        pileupweights[93] = 0;
        pileupweights[94] = 0;
        pileupweights[95] = 0;
        pileupweights[96] = 0;
        pileupweights[97] = 0;
        pileupweights[98] = 0;
        pileupweights[99] = 0; }

	if(RunYear == 2017){
	pileupweights[0] = 0;
        pileupweights[1] = 0;
        pileupweights[2] = 0.368589;
        pileupweights[3] = 5.52883;
        pileupweights[4] = 0.652119;
        pileupweights[5] = 0.700319;
        pileupweights[6] = 0.425295;
        pileupweights[7] = 0.455316;
        pileupweights[8] = 0.466681;
        pileupweights[9] = 0.436581;
        pileupweights[10] = 0.463855;
        pileupweights[11] = 0.483533;
        pileupweights[12] = 0.458275;
        pileupweights[13] = 0.476209;
        pileupweights[14] = 0.479196;
        pileupweights[15] = 0.528833;
        pileupweights[16] = 0.498161;
        pileupweights[17] = 0.592429;
        pileupweights[18] = 0.610549;
        pileupweights[19] = 0.616739;
        pileupweights[20] = 0.628302;
        pileupweights[21] = 0.659463;
        pileupweights[22] = 0.711438;
        pileupweights[23] = 0.738332;
        pileupweights[24] = 0.712789;
        pileupweights[25] = 0.770604;
        pileupweights[26] = 0.808979;
        pileupweights[27] = 0.822348;
        pileupweights[28] = 0.902373;
        pileupweights[29] = 0.903288;
        pileupweights[30] = 0.987908;
        pileupweights[31] = 1.03126;
        pileupweights[32] = 1.06761;
        pileupweights[33] = 1.10689;
        pileupweights[34] = 1.20031;
        pileupweights[35] = 1.2554;
        pileupweights[36] = 1.28501;
        pileupweights[37] = 1.39552;
        pileupweights[38] = 1.43056;
        pileupweights[39] = 1.47498;
        pileupweights[40] = 1.58264;
        pileupweights[41] = 1.67398;
        pileupweights[42] = 1.76129;
        pileupweights[43] = 1.65561;
        pileupweights[44] = 1.95346;
        pileupweights[45] = 2.0641;
        pileupweights[46] = 2.35761;
        pileupweights[47] = 2.31189;
        pileupweights[48] = 2.51311;
        pileupweights[49] = 2.64964;
        pileupweights[50] = 2.74274;
        pileupweights[51] = 3.17832;
        pileupweights[52] = 3.48279;
        pileupweights[53] = 4.42307;
        pileupweights[54] = 4.15277;
        pileupweights[55] = 4.2729;
        pileupweights[56] = 4.33371;
        pileupweights[57] = 5.15272;
        pileupweights[58] = 7.05585;
        pileupweights[59] = 5.57639;
        pileupweights[60] = 5.38363;
        pileupweights[61] = 7.39021;
        pileupweights[62] = 4.79166;
        pileupweights[63] = 9.70617;
        pileupweights[64] = 12.5781;
        pileupweights[65] = 15.9968;
        pileupweights[66] = 10.8734;
        pileupweights[67] = 9.32003;
        pileupweights[68] = 28.9342;
        pileupweights[69] = 45.705;
        pileupweights[70] = 9.58331;
        pileupweights[71] = 14.0064;
        pileupweights[72] = 16.2179;
        pileupweights[73] = 7.49464;
        pileupweights[74] = 14.1907;
        pileupweights[75] = 12.1634;
        pileupweights[76] = 8.29325;
        pileupweights[77] = 0;
        pileupweights[78] = 0;
        pileupweights[79] = 11.4263;
        pileupweights[80] = 0;
        pileupweights[81] = 0;
        pileupweights[82] = 0;
        pileupweights[83] = 5.16024;
        pileupweights[84] = 0;
        pileupweights[85] = 0;
        pileupweights[86] = 0;
        pileupweights[87] = 0;
        pileupweights[88] = 2.21153;
        pileupweights[89] = 0;
        pileupweights[90] = 0;
        pileupweights[91] = 0;
        pileupweights[92] = 0;
        pileupweights[93] = 0;
        pileupweights[94] = 2.58012;
        pileupweights[95] = 0;
        pileupweights[96] = 0;
        pileupweights[97] = 0;
        pileupweights[98] = 0;
        pileupweights[99] = 0;}

	if(RunYear == 2018){
	pileupweights[0] = 0;
	pileupweights[1] = 0;
        pileupweights[2] = 1.28953;
        pileupweights[3] = 0.421578;
        pileupweights[4] = 0.82745;
        pileupweights[5] = 0.523873;
        pileupweights[6] = 0.800097;
        pileupweights[7] = 0.575908;
        pileupweights[8] = 0.681356;
        pileupweights[9] = 0.631048;
        pileupweights[10] = 0.735895;
        pileupweights[11] = 0.767841;
        pileupweights[12] = 0.733907;
        pileupweights[13] = 0.854597;
        pileupweights[14] = 0.71963;
        pileupweights[15] = 0.768989;
        pileupweights[16] = 0.7891;
        pileupweights[17] = 0.811136;
        pileupweights[18] = 0.795095;
        pileupweights[19] = 0.80085;
        pileupweights[20] = 0.831406;
        pileupweights[21] = 0.835328;
        pileupweights[22] = 0.875346;
        pileupweights[23] = 0.847136;
        pileupweights[24] = 0.88979;
        pileupweights[25] = 0.901946;
        pileupweights[26] = 0.902386;
        pileupweights[27] = 0.950084;
        pileupweights[28] = 0.9262;
        pileupweights[29] = 0.936081;
        pileupweights[30] = 0.968535;
        pileupweights[31] = 0.989897;
        pileupweights[32] = 1.01124;
        pileupweights[33] = 1.03438;
        pileupweights[34] = 0.998615;
        pileupweights[35] = 1.04456;
        pileupweights[36] = 1.17434;
        pileupweights[37] = 1.1385;
        pileupweights[38] = 1.20681;
        pileupweights[39] = 1.29955;
        pileupweights[40] = 1.36529;
        pileupweights[41] = 1.5202;
        pileupweights[42] = 1.49041;
        pileupweights[43] = 1.67056;
        pileupweights[44] = 1.82047;
        pileupweights[45] = 2.1055;
        pileupweights[46] = 2.1692;
        pileupweights[47] = 2.59291;
        pileupweights[48] = 2.65146;
        pileupweights[49] = 2.43148;
        pileupweights[50] = 3.33404;
        pileupweights[51] = 3.04206;
        pileupweights[52] = 3.85581;
        pileupweights[53] = 4.45841;
        pileupweights[54] = 4.90291;
        pileupweights[55] = 5.29153;
        pileupweights[56] = 5.35156;
        pileupweights[57] = 5.23522;
        pileupweights[58] = 6.20588;
        pileupweights[59] = 9.15568;
        pileupweights[60] = 6.47697;
        pileupweights[61] = 5.22013;
        pileupweights[62] = 8.3551;
        pileupweights[63] = 9.86493;
        pileupweights[64] = 19.343;
        pileupweights[65] = 11.6058;
        pileupweights[66] = 6.57662;
        pileupweights[67] = 9.83269;
        pileupweights[68] = 4.5456;
        pileupweights[69] = 0;
        pileupweights[70] = 0;
        pileupweights[71] = 6.77005;
        pileupweights[72] = 14.6684;
        pileupweights[73] = 12.4118;
        pileupweights[74] = 11.2834;
        pileupweights[75] = 0;
        pileupweights[76] = 0;
        pileupweights[77] = 0;
        pileupweights[78] = 6.12528;
        pileupweights[79] = 2.65966;
        pileupweights[80] = 0;
        pileupweights[81] = 0;
        pileupweights[82] = 0;
        pileupweights[83] = 0;
        pileupweights[84] = 0;
        pileupweights[85] = 0;
        pileupweights[86] = 0;
        pileupweights[87] = 0;
        pileupweights[88] = 0;
        pileupweights[89] = 0;
        pileupweights[90] = 0;
        pileupweights[91] = 0;
        pileupweights[92] = 0;
        pileupweights[93] = 0;
        pileupweights[94] = 0;
        pileupweights[95] = 0;
        pileupweights[96] = 0;
        pileupweights[97] = 0;
        pileupweights[98] = 0;
        pileupweights[99] = 0;

	}

	if(nvertex > 99)return 0;
	else return pileupweights[nvertex]; 

}
void plot_PUReweight(){
	setTDRStyle();
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFillColor(0);
        gROOT->SetBatch(kTRUE);	
	TH1D *p_PU_data = new TH1D("p_PU_data",";N_{vtx};",100,0,100);
	TH1D *p_PU_MC = new TH1D("p_PU_MC","",100,0,100);
	TH1D *p_PU_MC_raw = new TH1D("p_PU_MC_raw","",100,0,100);

	TChain *sigtree = new TChain("signalTree");
	sigtree->Add(Form("/eos/uscms/store/group/lpcsusyhad/Tribeni/eg_mg_trees/resTree_mgsignal_MuonEG_%d_NEW.root",RunYear));
	sigtree->Draw("nVertex >> p_PU_data");
 
	TChain *mctree;
  mctree = new TChain("mgTree","mgTree");
	mctree->Add(Form("/eos/uscms/store/user/tmishra/VGamma/resTree_VGamma_DYJetsToLL_%d.root",RunYear));
  int   nVertex(0);
	float PUweight(1);
	float puWeight;
  mctree->SetBranchAddress("nVertex",   &nVertex);
	mctree->SetBranchAddress("PUweight",  &PUweight);
	for(unsigned ievt(0); ievt < mctree->GetEntries(); ievt++){
		mctree->GetEntry(ievt);
		puWeight = getPUESF(nVertex);
		p_PU_MC->Fill(nVertex,puWeight);
		p_PU_MC_raw->Fill(nVertex);
	}

	
	p_PU_data->Sumw2();
	p_PU_data->Scale(1.0/p_PU_data->Integral(1,100));
	p_PU_MC->Scale(1.0/p_PU_MC->Integral(1,100));
	p_PU_MC_raw->Scale(1.0/p_PU_MC_raw->Integral(1,100));
	
	gStyle->SetOptStat(0);
	TCanvas *c_PU = new TCanvas("PU", "PU",600,600);
	c_PU->cd();
	TPad *PU_pad1 = new TPad("PU_pad1", "PU_pad1", 0, 0.3, 1, 1.0);
	PU_pad1->SetBottomMargin(0);
	PU_pad1->Draw();  
	PU_pad1->cd(); 
	p_PU_data->SetMaximum(1.2*p_PU_data->GetBinContent(p_PU_data->GetMaximumBin()));
	p_PU_data->SetLineColor(1);
	p_PU_data->SetMarkerStyle(20);
	p_PU_data->GetXaxis()->SetTitle("N_{vtx}");
	p_PU_MC_raw->Draw("hist");
	p_PU_data->Draw("P same");
	p_PU_MC->SetLineColor(kMagenta);
	p_PU_MC_raw->SetLineColor(kBlue);
	p_PU_MC->Draw("hist same");
	TLegend *leg_PU =  new TLegend(0.5,0.65,0.9,0.9);
	leg_PU->SetFillStyle(0);
	leg_PU->AddEntry(p_PU_data,"Data");
	leg_PU->AddEntry(p_PU_MC_raw,"Simulation (no correction)");
	leg_PU->AddEntry(p_PU_MC,"Simulation (corrected)");
	leg_PU->Draw("same");
	prelim->Draw();
	if(RunYear==2016)lumitex->Draw();
	if(RunYear==2017)lumitex17->Draw();
	if(RunYear==2018)lumitex18->Draw();

	c_PU->cd();
	TPad *PU_pad2 = new TPad("PU_pad2", "PU_pad2", 0, 0.05, 1, 0.3);
	PU_pad2->SetTopMargin(0);
	PU_pad2->SetBottomMargin(0.3);
	PU_pad2->Draw();
	PU_pad2->cd();
  TLine *flatratio_PU = new TLine(0,1,100,1);
	TH1F *ratio_PU=(TH1F*)p_PU_data->Clone("transfer factor");
	ratio_PU->SetMarkerStyle(20);
	ratio_PU->SetLineColor(kBlack);
	ratio_PU->GetXaxis()->SetRangeUser(0,100);
	ratio_PU->GetYaxis()->SetRangeUser(0,4);
	ratio_PU->SetMinimum(0);
	ratio_PU->SetMaximum(3);
	//ratio_PU->Divide(p_PU_MC_raw); //to get the pileupweights
	ratio_PU->Divide(p_PU_MC);
	
	cout<< RunYear <<endl;
	for(int bin=0;bin<100;bin++){
		cout<<"        pileupweights["<<bin<<"] = "<<ratio_PU->GetBinContent(bin)<<";\n";
		//cout<<ratio_PU->GetBinContent(bin)<<"; ";
	}
	cout<<endl;
	ratio_PU->SetTitle("");
	ratio_PU->GetYaxis()->SetTitle("Data/MC");
	ratio_PU->Draw();
	flatratio_PU->Draw("same");
	c_PU->SaveAs(Form("/eos/uscms/store/user/tmishra/Plots_myAN/PLOT_PUreweight_%d.pdf",RunYear));

}
