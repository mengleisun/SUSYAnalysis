#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "../../../include/tdrstyle.C"

void plotFake(){//main  
gStyle->SetOptStat(0);
int PtBins[]={25,30,35,40,50,60,70,90};
float EtaBins[] = {0.0,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4};
float VertexBins[]={0,4,8,10,11,12,13,14,15,16,17,18,20, 22};
unsigned nPtBins = sizeof(PtBins)/sizeof(int);
unsigned nEtaBins = sizeof(EtaBins)/sizeof(float);
unsigned nVertexBins = sizeof(VertexBins)/sizeof(float);

TGraphErrors *fr_bothcount_pt = new TGraphErrors(nPtBins-1);
TGraphErrors *fr_bothcount_eta = new TGraphErrors(nEtaBins-1);
TGraphErrors *fr_bothcount_vtx = new TGraphErrors(nVertexBins-1);

double bothcount_num[] = {1811.68 ,2743.94 ,4484.98 ,5865.94 ,1063.21 ,316.909 ,362.357};
double bothcount_numerror[] = {89.1805 ,105.796 ,136.779 ,117.965 ,60.1579 ,27.0541 ,26.3509};
double bothcount_den[] = {88599.3 ,150759 ,235616 ,456774 ,91937.1 ,28127.8 ,25741.8};
double bothcount_denerror[] = {366.088 ,482.562 ,485.402 ,927.565 ,369.975 ,200.697 ,187.72};
                                                                      
double bothcount_eta_num[] = {1998.48 ,1122.35 ,993.577 ,969.641 ,894.844 ,826.551 ,961.078 ,741.578 ,677.432 ,593.952 ,492.355 ,246.616 ,241.061 ,431.28};
double bothcount_eta_numerror[] = {48.3773 ,37.4573 ,86.3437 ,34.9478 ,118.579 ,32.6808 ,73.7561 ,60.5476 ,29.7681 ,27.9801 ,25.8867 ,30.7935 ,54.6724 ,24.9071};
double bothcount_eta_den[] = {66416 ,67355.3 ,66413.3 ,62050.9 ,59206.1 ,58840.9 ,56211.1 ,49627.9 ,47983.9 ,41692 ,34122.8 ,23157 ,19390.4 ,23690};
double bothcount_eta_denerror[] = {260.183 ,262.157 ,260.924 ,252.388 ,246.19 ,244.842 ,240.092 ,227.239 ,223.342 ,207.912 ,188.941 ,155.586 ,142.739 ,158.035};

 double bothcount_vtx_num[] = {86.1448 ,2521.46 ,3233.16 ,1999.81 ,1754.09 ,1771.93 ,1309.58 ,1051.85 ,795.709 ,540.555 ,406.979 ,448.213 ,250.753};
 double bothcount_vtx_numerror[]={13.6964 ,87.0778 ,111.225 ,84.0988 ,71.1755 ,74.2046 ,107.329 ,52.3201 ,45.5077 ,40.458 ,37.8934 ,46.8365 ,25.0772};
 double bothcount_vtx_den[] = {7930.85 ,186939 ,227883 ,124658 ,119988 ,103696 ,86346.8 ,66832.6 ,48997.1 ,34452.3 ,23006.8 ,24173.5 ,13021.1};
 double bothcount_vtx_denerror[]={91.2509 ,440.953 ,486.002 ,360.401 ,352.044 ,329.385 ,298.967 ,263.399 ,224.839 ,188.176 ,154.428 ,158.028 ,115.929};

for(int i(0); i<nPtBins-1; i++){
  fr_bothcount_pt->SetPoint(i,(PtBins[i]+PtBins[i+1])/2.0,bothcount_num[i]/bothcount_den[i]);
  fr_bothcount_pt->SetPointError(i, (PtBins[i+1]-PtBins[i])/2.0, bothcount_num[i]*bothcount_denerror[i]/(bothcount_den[i]*bothcount_den[i])+bothcount_numerror[i]/bothcount_den[i]);
}

for(int i(0); i<nEtaBins-1; i++){
  fr_bothcount_eta->SetPoint(i,(EtaBins[i]+EtaBins[i+1])/2.0,bothcount_eta_num[i]/bothcount_eta_den[i]);
  fr_bothcount_eta->SetPointError(i, (EtaBins[i+1]-EtaBins[i])/2.0, bothcount_eta_num[i]*bothcount_eta_denerror[i]/(bothcount_eta_den[i]*bothcount_eta_den[i])+bothcount_eta_numerror[i]/bothcount_eta_den[i]);
}

for(int i(0); i<nVertexBins-1; i++){
  fr_bothcount_vtx->SetPoint(i,(VertexBins[i]+VertexBins[i+1])/2.0,bothcount_vtx_num[i]/bothcount_vtx_den[i]);
  fr_bothcount_vtx->SetPointError(i, (VertexBins[i+1]-VertexBins[i])/2.0, bothcount_vtx_num[i]*bothcount_vtx_denerror[i]/(bothcount_vtx_den[i]*bothcount_vtx_den[i])+bothcount_vtx_numerror[i]/bothcount_vtx_den[i]);
}

setTDRStyle();
gStyle->SetOptFit(0);
fr_bothcount_pt->SetMarkerStyle(20);
fr_bothcount_pt->SetMarkerColor(kBlue);
fr_bothcount_pt->SetLineColor(kBlue);
fr_bothcount_pt->SetLineWidth(2);
fr_bothcount_pt->SetFillColor(0);

fr_bothcount_eta->SetMarkerStyle(20);
fr_bothcount_eta->SetMarkerColor(kBlue);
fr_bothcount_eta->SetLineColor(kBlue);
fr_bothcount_eta->SetLineWidth(2);
fr_bothcount_eta->SetFillColor(0);

fr_bothcount_vtx->SetMarkerStyle(20);
fr_bothcount_vtx->SetMarkerColor(kBlue);
fr_bothcount_vtx->SetLineColor(kBlue);
fr_bothcount_vtx->SetLineWidth(2);
fr_bothcount_vtx->SetFillColor(0);

TH1F *dummy_pt = new TH1F("",";Pt(GeV);fake rate",10,20,90);
TH1F *dummy_eta = new TH1F("",";|eta|;fake rate",14,0,1.4);
TH1F *dummy_vtx = new TH1F("",";nVtx;fake rate",11,0,22);

//TFile* file = TFile::Open("crosscheck.root","read");
//TProfile *p_fakerate = (TProfile*)file->Get("p_fakerate");

TCanvas *canvas_pt = new TCanvas("fake-rate vs. p_{T}"," fake rate",600,600);
canvas_pt->cd();

dummy_pt->SetMaximum(0.03);
dummy_pt->Draw();
TF1 *f1 = new TF1("f1", "pow([0]*x+[1],[2])",25,130);
f1->SetParameter(0, 1);
f1->SetParameter(1, 0);
f1->SetParameter(2, -2.0);

fr_bothcount_pt->Draw("P same");
//fr_bothcount_pt->Fit("f1","R");
TCanvas *canvas_eta = new TCanvas("fake-rate vs. eta"," fake rate",600,600);
canvas_eta->cd();
dummy_eta->SetMaximum(0.04);
dummy_eta->Draw();
fr_bothcount_eta->Draw("P same");
TCanvas *canvas_vtx = new TCanvas("fake-rate vs. N_{vtx}"," fake rate",600,600);
canvas_vtx->cd();
dummy_vtx->SetMaximum(0.03);
dummy_vtx->Draw();
fr_bothcount_vtx->Fit("pol1");
fr_bothcount_vtx->Draw("P same");
//fr_bothcount_MC->Draw("P same");
//p_fakerate->Draw("same");
//fr_random->Draw("same");
//fr_hightag->Draw("same");
//TLegend *leg =  new TLegend(0.4,0.7,0.9,0.9);
//leg->SetFillStyle(0);
//gStyle->SetLegendBorderSize(1);
//gStyle->SetLegendFillColor(0);
//gStyle->SetTextSize(4);
//leg->AddEntry(fr_bothcount,"data fake rate: DoubleEG");
//leg->AddEntry(fr_bothcount_MC,"MC fake rate: DYJetToLL");
//leg->Draw("same");

canvas_pt->SaveAs("EleFakePho_pt_data_76X.pdf");
canvas_pt->SaveAs("EleFakePho_pt_data_76X.png");
//canvas_eta->SaveAs("EleFakePho_eta_data_76X.pdf");
//canvas_eta->SaveAs("EleFakePho_eta_data_76X.png");
canvas_vtx->SaveAs("EleFakePho_vtx_data_76X.pdf");
canvas_vtx->SaveAs("EleFakePho_vtx_data_76X.png");

}


