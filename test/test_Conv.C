#include "../include/analysis_commoncode.h"
#include "../include/analysis_cuts.h"

#define debugMCLevel 4
enum branchType{
	gg = 1,
	gZ = 2,
	gH = 3,
	gW = 4,
	ZZ = 5,
	HH = 6,
	noType = 7	
};

struct decayChain{
	std::vector<mcData>::iterator iter; 
  std::vector< std::vector<mcData>::iterator > daughter;
};

void test_Conv(){//main  

	gSystem->Load("/uscms/home/mengleis/work/SUSY2016/SUSYAnalysis/lib/libAnaClasses.so");

	TChain* es = new TChain("ggNtuplizer/EventTree");
	es->Add("root://cmseos.fnal.gov///store/user/msun/Signal/SMS-TChiNG_BF50N50G_TuneCUETP8M1.root");

	RunType datatype(MC); 
	
	rawData raw(es, datatype);
	std::vector<mcData>  MCData;
	std::vector<recoPhoton> Photon;
	std::vector<recoEle> Ele;
	std::vector<recoMuon>   Muon;
	std::vector<recoJet>   JetCollection;
	float met(0);
	float metPhi(0);
	float met_T1JERUp(0);
	float met_T1JERDo(0);
	float met_T1JESUp(0);
	float met_T1JESDo(0);	
	float	metPhi_T1JESUp(0);
	float	metPhi_T1JESDo(0);
	float	metPhi_T1UESUp(0);
	float	metPhi_T1UESDo(0);
	int nVtx(0);
	int jetNumber(0);

	const unsigned nEvts = es->GetEntries()/10; 
	std::cout << "total event : " << nEvts << std::endl;

	for(unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
  
		if (ievt%10000==0) std::cout << " -- Processing event " << ievt << std::endl;
		raw.GetData(es, ievt);
		MCData.clear();
		Photon.clear();
		Muon.clear();
		Ele.clear();
		JetCollection.clear();
		if(datatype == MC)for(int iMC(0); iMC < raw.nMC; iMC++){MCData.push_back(mcData(raw, iMC));}
		for(int iPho(0); iPho < raw.nPho; iPho++){Photon.push_back(recoPhoton(raw, iPho));}
		for(int iMu(0); iMu < raw.nMu; iMu++){Muon.push_back(recoMuon(raw, iMu));}
		for(int iEle(0); iEle < raw.nEle; iEle++){Ele.push_back(recoEle(raw, iEle));}
		for(int iJet(0); iJet < raw.nJet; iJet++){JetCollection.push_back(recoJet(raw, iJet));}

		std::vector< std::vector<mcData>::iterator > genNegList;
		std::vector< std::vector<mcData>::iterator > genPosList;
		std::vector< std::vector<mcData>::iterator > genTauList;
		std::vector< std::vector<mcData>::iterator > genPhotonList;
		genNegList.clear();
		genPosList.clear();
		genPhotonList.clear();
		genTauList.clear();
		std::vector< std::vector<mcData>::iterator > itMCList;
		itMCList.clear();	
		std::vector< pair<int, int> > finalState;
		finalState.clear();
		for(std::vector<mcData>::iterator itMC = MCData.begin(); itMC!= MCData.end(); itMC++){ 
				bool inList(false);
				for(unsigned imc(0); imc < itMCList.size(); imc++){
					if( itMC->getPID() == itMCList[imc]->getPID() && itMC->getMomPID() == itMCList[imc]->getMomPID() && DeltaR(itMC->getEta(), itMC->getPhi(), itMCList[imc]->getEta(), itMCList[imc]->getPhi()) < 0.3)inList=true;
				}
				if(!inList){
			itMCList.push_back(itMC);

			if(itMC->getPID() == 22)genPhotonList.push_back(itMC);
			if(itMC->getPID()	== 11)genNegList.push_back(itMC);
			if(itMC->getPID() == -11)genPosList.push_back(itMC);
			if(itMC->getPID() == 15 || itMC->getPID() == -15 || itMC->getPID() == 13 || itMC->getPID() == -13)genTauList.push_back(itMC); 
			if(fabs(itMC->getPID()) > 20 && fabs(itMC->getPID()) < 26 && fabs(itMC->getMomPID()) > 1000000){
				finalState.push_back( make_pair(fabs(itMC->getPID()), fabs(itMC->getMomPID())) );
			}
			}
		}//loop on MC particles
		int decaybranch(-1);
		bool hasBranchPho(false);
		bool hasBranchZ(false);
		bool hasBranchH(false);
		bool hasBranchW(false);
		for(int i(0); i < finalState.size(); i++){
			if(finalState[i].first == 22)hasBranchPho = true;
			if(finalState[i].first == 23)hasBranchZ = true;
			if(finalState[i].first == 24)hasBranchW = true;
			if(finalState[i].first == 25)hasBranchH = true;
		}
		if(hasBranchPho && hasBranchZ)decaybranch = 2; 
		else if(hasBranchPho && hasBranchH)decaybranch = 3;
		else if(hasBranchPho && hasBranchW)decaybranch = 4; 
		else if(hasBranchPho && !hasBranchZ && !hasBranchH && !hasBranchW)decaybranch = 1; 
   

		if(decaybranch == 1){ 
				int NrecoPho(0);
					for(std::vector<recoPhoton>::iterator itpho = Photon.begin() ; itpho != Photon.end(); ++itpho){
						if(!itpho->passSignalSelection())continue;
						bool PixelVeto = itpho->PixelSeed()==0? true: false;
						if(PixelVeto ){
							NrecoPho += 1;
						}
					}
					std::cout << "gg pho " << NrecoPho << std::endl;

		for(std::vector<recoEle>::iterator itEle = Ele.begin(); itEle != Ele.end(); itEle++){
			if((itEle->isEB() && itEle->getR9() < R9EBCut) || (itEle->isEE() && itEle->getR9() < R9EECut))continue;
	if(itEle->isEB() && itEle->getSigma() > 0.00998)continue;
	if(itEle->isEB() && fabs(itEle->getdEtaIn()) > 0.00311)continue;
	if(itEle->isEB() && fabs(itEle->getdPhiIn()) > 0.103)continue;
	if(itEle->isEB() && itEle->getHoverE() > 0.253)continue;
	if(itEle->isEB() && fabs(itEle->getEoverPInv()) > 0.134)continue;
	if(itEle->isEB() && itEle->getMissHits() > 1)continue;
	if(itEle->isEB() && itEle->getMiniIso() > 0.1)continue;
	if(itEle->isEE() && itEle->getSigma() > 0.0298)continue;
	if(itEle->isEE() && fabs(itEle->getdEtaIn()) > 0.00609)continue;
	if(itEle->isEE() && fabs(itEle->getdPhiIn()) > 0.045)continue;
	if(itEle->isEE() && itEle->getHoverE() > 0.0878)continue;
	if(itEle->isEE() && fabs(itEle->getEoverPInv()) > 0.13)continue;
	if(itEle->isEE() && itEle->getMissHits() > 1)continue;
	if(itEle->isEE() && itEle->getMiniIso() > 0.1)continue;

			if(itEle->getPt() > 25){
				bool matchEle(false);
				for(unsigned imc(0); imc < genNegList.size(); imc++){
					if(DeltaR(itEle->getEta(), itEle->getPhi(), genNegList[imc]->getEta(), genNegList[imc]->getPhi()) < 0.2){std::cout << "conv " << itEle->getConvVeto() << " match e-  Mom:" << genNegList[imc]->getMomPID() << std::endl; matchEle = true;}
				}
				for(unsigned imc(0); imc < genPosList.size(); imc++){
					if(DeltaR(itEle->getEta(), itEle->getPhi(), genPosList[imc]->getEta(), genPosList[imc]->getPhi()) < 0.2){std::cout << "conv " << itEle->getConvVeto() << " match e+  Mom:" << genPosList[imc]->getMomPID() << std::endl; matchEle = true;}
				}
				for(unsigned imc(0); imc < genTauList.size(); imc++){
					if(DeltaR(itEle->getEta(), itEle->getPhi(), genTauList[imc]->getEta(), genTauList[imc]->getPhi()) < 0.2){std::cout << "conv " << itEle->getConvVeto() << " match tau  Mom:" << genTauList[imc]->getMomPID() << std::endl; matchEle = true;}
				}
				if(!matchEle){
				for(unsigned imc(0); imc < genPhotonList.size(); imc++){
					if(DeltaR(itEle->getEta(), itEle->getPhi(), genPhotonList[imc]->getEta(), genPhotonList[imc]->getPhi()) < 0.05)std::cout << "conv " << itEle->getConvVeto() << "match pho  Mom:" << genPhotonList[imc]->getMomPID() << std::endl;
				}	
				}
			}
		}

		std::cout << std::endl;
		}

	}//loop on entries

}
