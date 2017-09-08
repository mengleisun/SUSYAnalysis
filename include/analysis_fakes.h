//double etaRates[]={0.0491495,0.0316066,0.0253055,0.024221,0.0223884,0.0222271,0.0210704,0.0203349,0.0207721,0.0198487,0.0193914,0.0194545,0.0181131,0.0244214,0.0241276,0.0219053,0.0204027,0.0195225,0.0198563,0.0191962,0.0173755,0.0181452,0.0299651,0.0222873,0.0188254,0.0185798,0.0200921,0.0235761,0.0296801};   //nominal value

//double etaRates[]={0.220502,0.140268,0.10809,0.103675,0.097537,0.0952774,0.0908136,0.0850506,0.0866651,0.0835787,0.0757995,0.0753234,0.068137,0.0889137,0.084842,0.0757298,0.0660832,0.0620414,0.0606887,0.0575167,0.0498608,0.049543,0.0793691,0.0565387,0.0459593,0.0437259,0.0464076,0.0531607,0.0645416};  // with 0.02 FSR

double etaRates[]={0.0492472,0.0317083,0.0249878,0.0237811,0.0225391,0.0217952,0.021218,0.0204555,0.0209017,0.0213058,0.0195457,0.0195127,0.0182495,0.0245207,0.0241966,0.0217942,0.0204352,0.0195351,0.0197452,0.0191388,0.0175107,0.0181074,0.0299934,0.0222525,0.0188657,0.0184809,0.0201496,0.0238605,0.0297451,0.363044,0.023322,0.0245266,0.0229644,0.0256796,0.0238304,0.0238648,0.023646,0.0245623,0.020066,0.0159119,0.0142946,0.0147044,0.0145453,0.0136805,0.0139513,0.0154788,0.0385528};

Double_t fakerate_func(Double_t *x, Double_t *par)
{
	double weight_pt(0.02);
	double weight_nvtx(0.02);
	double weight_eta(0.02);

	double normfactor = par[0];
  double slope = par[1];
  double constant = par[2];
  double index = par[3];
  double coeff = 1.0; 
	double vtx_constant = par[5];
	double vtx_slope = par[6];

  double pt = TMath::Max(x[0],0.000001);
	double nvtx=x[1];
	double eta= x[2];

	//if(pt >200)pt= 200.0;

	double arg = 0;
	arg = slope*pt + constant; 
	double fitval = pow(arg, index)*coeff; 
	weight_pt = fitval;

	weight_nvtx = vtx_constant + vtx_slope*nvtx;


	for(int ieta(0); ieta < 29; ieta++)
		if(eta > ieta*0.05 && eta < (ieta+1)*0.05)weight_eta = etaRates[ieta];
	if(eta <= 0 || eta > 1.4)weight_eta = 0.025;

	double totalfakerate= normfactor*weight_pt*weight_nvtx*weight_eta;
	//return totalfakerate/(1-totalfakerate);
	return totalfakerate;

//1275 1565
}

Double_t jetfake_func(Double_t *x, Double_t *par)
{
	double pt_low = x[0] - 1.0/2.0;
	double pt_high = x[0]+ 1.0/2.0;

	double c1 = par[0];
	double c2 = par[1];
	double lamda1 = par[2];
	double lamda2 = par[3];

	double jetfakes_lowedge = c1*exp(lamda1*pt_low)/lamda1 + c2*exp(lamda2*pt_low)/lamda2;
	double jetfakes_highedge =  c1*exp(lamda1*pt_high)/lamda1 + c2*exp(lamda2*pt_high)/lamda2;
	//return (jetfakes_highedge + jetfakes_lowedge)/2.0*REBINSIZE;
	return (jetfakes_highedge - jetfakes_lowedge);
}

