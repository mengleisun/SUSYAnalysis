void VariableComp( const TString &graphicsFormat = "png")
 {
    int RunYear = 2018;
    const char* channel = "TChiWG";

    TFile *f2 = new TFile(Form("T5Wg_%d_slim.root",RunYear));
    TFile *f1 = new TFile(Form("TChiWG_%d_slim.root",RunYear));
    const unsigned int kNDists = 24;
    TH1* hT5wg[kNDists];
    TH1* hTchiwg[kNDists];
   // TH1* h2018[kNDists];
    TH1* hRatio[kNDists];
    //for(int i = 0; i < kNDists; ++i) {    //scanning the root file
    for(int i = 0; i < 12; ++i) {    //scanning the root file
    	TString name = "",name1="";
    	if( i == 0){ name = "njets",name1="nJets";}
 	else if (i == 1){ name = "nelectrons",name1 = "nElectrons";}
 	else if (i == 2){ name = "nmuons",name1 = "nMuons";}
 	else if (i == 3){ name = "nphotons",name1 = "nPhotons";}
	else if( i == 4){ name = "jetpt",name1="Jet Pt";}
 	else if (i == 5){ name = "electronpt",name1 = "Electron Pt";}
 	else if (i == 6){ name = "muonpt",name1 = "Muon Pt";}
 	else if (i == 7){ name = "photonpt",name1 = "Photon Pt";}
	else if( i == 8){ name = "jetMass",name1="Jet Mass";}
 	else if (i == 9){ name = "electronMass",name1 = "Electron Mass";}
 	else if (i == 10){ name = "muonMass",name1 = "Muon Mass";}
 	else if (i == 11){ name = "photonMass",name1 = "Photon Mass";}
	
	
	
	
    	hTchiwg[i] = static_cast<TH1*>(f1->Get("h"+name));
    	hT5wg[i] = static_cast<TH1*>(f2->Get("h"+name));
    	hTchiwg[i]->Scale(1/(hTchiwg[i]->Integral()));
   	hT5wg[i]->Scale(1/hT5wg[i]->Integral());
	hTchiwg[i]->SetTitle("");
    	hTchiwg[i]->SetMarkerStyle(20);
    	hTchiwg[i]->SetLineColor(kRed);
    	hT5wg[i]->SetLineColor(kBlue);
    	hTchiwg[i]->SetLineWidth(2);   
    	hT5wg[i]->SetLineWidth(2); 

    	hTchiwg[i]->SetMarkerSize(0.2);
    	hTchiwg[i]->GetYaxis()->SetTitleOffset(0.65);
    	hTchiwg[i]->GetXaxis()->SetTitleOffset(0.65);

    	hTchiwg[i]->GetYaxis()->SetTitleFont(42);
   	hTchiwg[i]->GetYaxis()->SetLabelSize(0.05);
    	hTchiwg[i]->GetYaxis()->SetLabelFont(42);

	hT5wg[i]->SetMarkerStyle(222);
    	hT5wg[i]->GetXaxis()->SetTitle(name1);
    	hTchiwg[i]->GetXaxis()->SetTitle(name1);
    	hT5wg[i]->GetXaxis()->SetTitleSize(.08);
    	hTchiwg[i]->GetXaxis()->SetTitleSize(.08);
    	hT5wg[i]->SetTitleSize(0.08);
   	hTchiwg[i]->SetTitleSize(0.08);
   	hT5wg[i]->GetYaxis()->SetTitle("# of events    ");
  	hTchiwg[i]->GetYaxis()->SetTitle("# of events     ");
  	hT5wg[i]->GetYaxis()->SetTitleSize(.07);
	hTchiwg[i]->GetYaxis()->SetTitleSize(.07);
    	hT5wg[i]->SetMarkerSize(0.01);
    	gStyle->SetStatY(0.98);
    	gStyle->SetStatX(0.98);
    	gStyle->SetStatW(0.3);
    	gStyle->SetStatH(0.25);
    	
   	hTchiwg[i]->GetXaxis()->SetLabelSize(0.07);

	hRatio[i] = (TH1*) hTchiwg[i]->Clone("hRatio[i]");
    	hRatio[i]->Divide(hT5wg[i]);
 

	hRatio[i]->GetYaxis()->SetTitle(" #frac{Tchiwg}{T5wg}      ");
    	hRatio[i]->GetYaxis()->SetRangeUser(0,2);
    	hRatio[i]->SetStats(0);
    	hRatio[i]->SetMarkerStyle(20);
    	hRatio[i]->SetLineWidth(2);
    	hRatio[i]->SetLineColor(kBlack);
    	hRatio[i]->SetMarkerColor(kBlack);
    	hRatio[i]->GetYaxis()->SetTitleSize(0.1);

    	hRatio[i]->GetYaxis()->SetTitleOffset(0.4);
    	hRatio[i]->GetYaxis()->SetTitleFont(42);
    	hRatio[i]->GetYaxis()->SetLabelSize(0.09);
    	hRatio[i]->GetYaxis()->SetLabelFont(42);
    	hRatio[i]->GetXaxis()->SetLabelOffset(0.02);
    	hRatio[i]->GetXaxis()->SetLabelFont(42);
    	hRatio[i]->GetXaxis()->SetLabelSize(0.1);
    	hRatio[i]->GetXaxis()->SetTitleSize(0.15);
    	hRatio[i]->GetXaxis()->SetTitleFont(42);
    	hRatio[i]->GetXaxis()->SetTitleOffset(0.8);
	
}


    TCanvas *can[12];
    for(unsigned int i = 0; i<12; i++ ) {
                   
                           can[i]= new TCanvas(Form("c%d",i),"",1100,1200);
			   can[i]->cd();
			   const char* name = ""; 
		             	if( i == 0){ name="NJets";}
 			     	else if (i == 1){ name = "nElectrons";}
 				else if (i == 2){ name = "nMuons";}
 				else if (i == 3){ name = "nPhotons";}
				else if( i == 4){ name ="JetPt";}
 				else if (i == 5){ name = "ElectronPt";}
 				else if (i == 6){ name = "MuonPt";}
 				else if (i == 7){ name = "PhotonPt";}
				else if( i == 8){ name ="JetMass";}
 				else if (i == 9){ name = "ElectronMass";}
 				else if (i == 10){name = "MuonMass";}
 				else if (i == 11){name = "PhotonMass";}

			   hT5wg[i]->GetYaxis()->SetTitleFont(42);
			   float   hmax=0, hmin=0;
   
			   Float_t h1=hTchiwg[i]->GetMaximum();
			   Float_t h2=hT5wg[i]->GetMaximum();
        
			   if(h1>=h2)                             hmax=10*h1;     
			   else if(h1<h2)                         hmax=10*h2;
   
			   hT5wg[i]->SetAxisRange(0.003,hmax,"Y");
			   hTchiwg[i]->SetAxisRange(0.003,hmax,"Y");
			 //  can[j]->cd(i+1-4*j);
			 //
			   TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1.0, 1.0);
			   pad1->SetTopMargin(0.12);   
			   pad1->SetBottomMargin(0); // Upper and lower plot are joined
			   pad1->SetLeftMargin(0.10);
			   pad1->SetRightMargin(0.02);
   
			   pad1->Draw();
			   pad1->cd();     
			   pad1->SetLogy();
			   hTchiwg[i]->SetStats(0);
			   hT5wg[i]->SetStats(0);
			   hTchiwg[i]->Draw("hist");
			   hT5wg[i]->Draw("hist same");

			   TLegend *leg = new TLegend(.65,.68,.95,.88,"");
			   leg->SetFillColor(0);
			   leg->AddEntry(hTchiwg[i],"TChiWG ","l");
			   leg->AddEntry(hT5wg[i],"T5Wg ","l");
			   leg->SetFillStyle(0);
			   leg->SetBorderSize(0);
			   leg->SetTextSize(0.05);
			   leg->Draw("same");

	//		   can[j]->cd(i+1-4*j);
				can[i]->cd();
			   TPad *pad2 = new TPad("pad2", "pad2", 0, 0.01, 1.0, 0.3);
			   pad2->SetTopMargin(0);
			   pad2->SetBottomMargin(0.3);
	//		   pad2->SetLeftMargin(0.1);
			   pad2->SetRightMargin(0.02);
			   pad2->Draw();
			   pad2->cd();       // pad2 becomes the current pad


			   hRatio[i]->SetMarkerSize(0.7);
			   hRatio[i]->SetMarkerStyle(20);
			   hRatio[i]->GetYaxis()->SetRangeUser(-0.1,2.1);
			   hRatio[i]->SetMarkerColor(kBlack);
			   hRatio[i]->Draw("P SAME");

                           TLegend *legRatio = new TLegend(.3,.32,.6,.52,"");
                           legRatio->SetFillColor(0);
                           legRatio->AddEntry(hRatio[i],"#frac{TChiWG}{T5Wg} ","l");
                           legRatio->SetFillStyle(0);
                           legRatio->SetBorderSize(0);
			   legRatio->SetTextSize(0.09);
                         //  legRatio->Draw("same");


			   Double_t Gmax=hRatio[i]->GetXaxis()->GetXmax();
			   Double_t Gmin=hRatio[i]->GetXaxis()->GetXmin();
			   TLine *line = new TLine(Gmin,1,Gmax,1);
			   line->SetLineColor(kBlack);
			   line->Draw();
    			   can[i]->SaveAs(Form("TChiWG_vs_T5Wg_%s_%d.png",name,RunYear));
 		  }
   
}

