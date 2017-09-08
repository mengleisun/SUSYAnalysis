#include<string>
#include<iostream>
#include<fstream>
#include<sstream>

#include "TCanvas.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TText.h"
#include "TLatex.h"

void plotXS(){
std::ifstream inputFile;
std::ostringstream filename;

gStyle->SetOptStat(0);
TCanvas *can1 = new TCanvas("cross section","cross section",800,600);
can1->SetFillStyle(4000);
can1->SetFrameFillStyle(4000); 
TPad *pad1 = new TPad("pad1","",0,0,1,1);
pad1->SetFillStyle(4000); 
pad1->SetFrameFillStyle(4000); 
pad1->Draw();

TH2D *p_XS = new TH2D("p_XS","13Tev cross section; M_{3}(GeV); M_{2}(GeV)",18,700,1600,28,200,1600);
TH1F *dummy = new TH1F("dummy",";m_{average} [GeV]",1,100,1500);
TGraph *p_nn8XS = new TGraph;
TGraph *p_gg8XS = new TGraph;
TGraph *p_nn13XS = new TGraph;
TGraph *p_gg13XS = new TGraph;

std::ofstream nn13XS;
nn13XS.open("13TeVEWcs.dat");
std::ofstream gg13XS;
gg13XS.open("13TeVStrongcs.dat");
std::ofstream nn8XS;
nn8XS.open("8TeVEWcs.dat");
std::ofstream gg8XS;
gg8XS.open("8TeVStrongcs.dat");


string runtype("");
float datagg[15];
float datancp[15];
float datancn[15];

for(unsigned m3(715); m3 <= 1615;){
  for(unsigned m2(105); m2 <= m3;){

    runtype="";
    for(unsigned ii(0); ii<15; ii++)datagg[ii]=0;
 
    filename.str("");
    filename << "./NEW/13TeVXS/gg/13TeVgg_M3_" << m3 << "_M2_" << m2 << ".dat";
    inputFile.open(filename.str());

    inputFile >> runtype;
    for(unsigned jj(0); jj < 15; jj++)inputFile >> datagg[jj];

    p_XS->Fill(m3, m2, datagg[14]);
    if(m2==105){
      p_gg13XS->SetPoint((m3-715)/50, (datagg[5]+datagg[6])/2, datagg[14]);
      gg13XS << runtype << " " << (datagg[5]+datagg[6])/2 << " " << datagg[14] << "\n"; 
    }
    inputFile.close();
   
    if(m3==1615){
      runtype="";
      for(unsigned ii(0); ii<15; ii++){datancp[ii]=0;datancn[ii]=0;}

      filename.str("");
      filename << "./NEW/13TeVXS/nn/13TeVncp_M3_" << m3 << "_M2_" << m2 << ".dat";
      inputFile.open(filename.str());
     
      inputFile >> runtype;
      for(unsigned jj(0); jj < 15; jj++)inputFile >> datancp[jj];
      inputFile.close(); 
      
      filename.str("");
      filename << "./NEW/13TeVXS/nn/13TeVncn_M3_" << m3 << "_M2_" << m2 << ".dat";
      inputFile.open(filename.str());

      inputFile >> runtype;
      for(unsigned jj(0); jj < 15; jj++)inputFile >> datancn[jj];
      inputFile.close(); 

      for(unsigned m(715); m<=m2+25; m+=50)p_XS->Fill(m, m2, datancp[14]+ datancn[14]);
      p_nn13XS->SetPoint((m2-105)/50, (datancp[5]+datancp[6])/2, datancp[14]+datancn[14]);
      nn13XS << runtype << " " << (datancp[5]+datancp[6])/2 << " " << datancp[14] << " " << datancn[14]<< "\n"; 
   }
   
   m2+=50;
  }
  m3+=50;
}


for(unsigned m3(715); m3 <= 1615;){
  for(unsigned m2(105); m2 <= m3;){

    runtype="";
    for(unsigned ii(0); ii<15; ii++)datagg[ii]=0;
 
    filename.str("");
    filename << "./NEW/8TeVXS/gg/8TeVgg_M3_" << m3 << "_M2_" << m2 << ".dat";
    inputFile.open(filename.str());

    inputFile >> runtype;
    for(unsigned jj(0); jj < 15; jj++)inputFile >> datagg[jj];

    if(m2==105){
      p_gg8XS->SetPoint((m3-715)/50, (datagg[5]+datagg[6])/2, datagg[14]);
      gg8XS << runtype << " " << (datagg[5]+datagg[6])/2 << " " << datagg[14] << "\n"; 
    }
    inputFile.close();
   
    if(m3==1615){
      runtype="";
      for(unsigned ii(0); ii<15; ii++){datancp[ii]=0;datancn[ii]=0;}

      filename.str("");
      filename << "./NEW/8TeVXS/nn/8TeVncp_M3_" << m3 << "_M2_" << m2 << ".dat";
      inputFile.open(filename.str());
     
      inputFile >> runtype;
      for(unsigned jj(0); jj < 15; jj++)inputFile >> datancp[jj];
      inputFile.close(); 
      
      filename.str("");
      filename << "./NEW/8TeVXS/nn/8TeVncn_M3_" << m3 << "_M2_" << m2 << ".dat";
      inputFile.open(filename.str());

      inputFile >> runtype;
      for(unsigned jj(0); jj < 15; jj++)inputFile >> datancn[jj];
      inputFile.close(); 

      p_nn8XS->SetPoint((m2-105)/50, (datancp[5]+datancp[6])/2, datancp[14]+datancn[14]);
      nn8XS << runtype << " " << (datancp[5]+datancp[6])/2 << " " << datancp[14] << " " << datancn[14]<< "\n"; 
   }
   
   m2+=50;
  }
  m3+=50;
}


pad1->cd();
//gPad->SetLogy();
//dummy->SetMaximum(50);
dummy->SetMaximum(20);
dummy->SetMinimum(0.001);
dummy->Draw();

p_nn8XS->SetLineWidth(3);
p_nn8XS->SetLineColor(kBlue);
p_nn8XS->SetLineStyle(3);

p_gg8XS->SetLineWidth(3);
p_gg8XS->SetLineColor(kRed);
p_gg8XS->SetLineStyle(3);

p_gg8XS->Draw("same");
p_nn8XS->Draw("same");

p_nn13XS->SetLineWidth(3);
p_nn13XS->SetLineColor(kBlue);
p_gg13XS->SetLineWidth(3);
p_gg13XS->SetLineColor(kRed);
p_nn13XS->Draw("same");
p_gg13XS->Draw("same");

TLatex *xlabel = new TLatex();
TLatex *xlabe2 = new TLatex(); 
TLatex *xlabe3 = new TLatex();
TLatex *xlabe4 = new TLatex();

xlabel->SetTextSize(0.05);
xlabel->SetTextColor(kBlue);
xlabel->DrawLatex(500, 0.004, "#tilde{#chi}_{1}^{0}#tilde{#chi}_{1}^{#pm} 8TeV");
xlabe2->SetTextSize(0.05);
xlabe2->SetTextColor(kBlue);
xlabe2->DrawLatex(500, 0.07,"#tilde{#chi}_{1}^{0}#tilde{#chi}_{1}^{#pm} 13TeV");
xlabe3->SetTextSize(0.05);
xlabe3->SetTextColor(kRed);
xlabe3->DrawLatex(1100, 0.007,"#tilde{g}#tilde{g} 8TeV");
xlabe4->SetTextSize(0.05);
xlabe4->SetTextColor(kRed);
xlabe4->DrawLatex(1100, 0.3,"#tilde{g}#tilde{g} 13TeV");

TLatex *title = new TLatex();
title->DrawLatex(900,10,"#sigma_{tot}[pb]: pp #rightarrow SUSY");

//can1->SaveAs("cross.pdf");
TCanvas *can2 = new TCanvas("cross section 2","cross section",600,600);
can2->cd();
gPad->SetLogz();
p_XS->Draw("colz");
//can2->SaveAs("GMSB_XS.pdf");  

}
  
    


