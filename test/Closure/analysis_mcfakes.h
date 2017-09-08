double etaRates[]={0.0484152,0.0330344,0.0274351,0.0279825,0.0268659,0.0255559,0.0255363,0.0243439,0.0238033,0.0219459,0.0210441,0.0190666,0.0198119,0.019615,0.0174586,0.0184199,0.0174611,0.0153868,0.0161515,0.0126975,0.0125806,0.0127848,0.0127544,0.0121073,0.0132903,0.014223,0.0142164,0.0145537,0.0175709};

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
	if(eta <= 0 || eta > 1.4)weight_eta = 0.018;

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
