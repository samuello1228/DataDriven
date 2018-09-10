// C++
#include <iostream>
using namespace std;
// ROOT
#include "TChain.h"
#include "TString.h"
#include "TFile.h"
#include "TH2D.h"
// My packages
#include "FakeRate/susyEvts.h"

Float_t totweight;
Float_t lumiScaling;
int MCId;

int NlepBL;
int nSigLep;
int isTruthLep1;
int isTruthLep2;
int nJets20;
int nBJets20;
Float_t Pt_l;
Float_t Pt_subl;
Float_t DeltaEtaLep;
Float_t met;
Float_t mt;
Float_t meff;
Float_t mljj_comb;
Float_t MT2;

int main()
{
	const TString CHAIN_NAME = "WZ_nom";
	
	//read old tree
	TString path = "/eos/user/d/dkoeck/WHSS/HistFitterTrees/Trees/IncludingLepTruth/Backgrounds_inclTruth_tWZ_tH.36100.root";
	TChain* tree1 = new TChain(CHAIN_NAME);
	tree1->Add(path.Data());
	cout << "There are " << tree1->GetEntries() << " events" << endl;

	tree1->SetBranchAddress("totweight", &totweight);
	tree1->SetBranchAddress("lumiScaling", &lumiScaling);
	tree1->SetBranchAddress("MCId", &MCId);

	tree1->SetBranchAddress("NlepBL", &NlepBL);
	tree1->SetBranchAddress("nSigLep", &nSigLep);
	tree1->SetBranchAddress("isTruthLep1", &isTruthLep1);
	tree1->SetBranchAddress("isTruthLep2", &isTruthLep2);
	tree1->SetBranchAddress("nJets20", &nJets20);
	tree1->SetBranchAddress("nBJets20", &nBJets20);
	tree1->SetBranchAddress("Pt_l", &Pt_l);
	tree1->SetBranchAddress("Pt_subl", &Pt_subl);
	tree1->SetBranchAddress("DeltaEtaLep", &DeltaEtaLep);
	tree1->SetBranchAddress("met", &met);
	tree1->SetBranchAddress("mt", &mt);
	tree1->SetBranchAddress("meff", &meff);
	tree1->SetBranchAddress("mljj_comb", &mljj_comb);
	tree1->SetBranchAddress("MT2", &MT2);
	//tree1->SetBranchAddress("", &);

	std::vector<TString> cutflowList;
	cutflowList.push_back("=2BaseLep && =2SigLep && SS");
	cutflowList.push_back("=2BaseLep && =2SigLep && SS,w");
	cutflowList.push_back("isTruthLep1");
	cutflowList.push_back("isTruthLep1,w");
	cutflowList.push_back("isTruthLep2");
	cutflowList.push_back("isTruthLep2,w");
	cutflowList.push_back("nJet");
	cutflowList.push_back("nJet,w");
	cutflowList.push_back("nBJets20");
	cutflowList.push_back("nBJets20,w");
	cutflowList.push_back("pt1");
	cutflowList.push_back("pt1,w");
	cutflowList.push_back("pt2");
	cutflowList.push_back("pt2,w");
	cutflowList.push_back("DeltaEtaLep");
	cutflowList.push_back("DeltaEtaLep,w");
	cutflowList.push_back("met");
	cutflowList.push_back("met,w");
	cutflowList.push_back("mt");
	cutflowList.push_back("mt,w");
	cutflowList.push_back("meff");
	cutflowList.push_back("meff,w");
	cutflowList.push_back("mljj");
	cutflowList.push_back("mljj,w");
	cutflowList.push_back("MT2");
	cutflowList.push_back("MT2,w");
	//cutflowList.push_back("");
	//cutflowList.push_back(",w");

	TH1D* hCutflow = new TH1D("cutflow", "cutflow", cutflowList.size(), 0, cutflowList.size());
	for(unsigned int i=0;i<cutflowList.size();i++)    
	{
		hCutflow->GetXaxis()->SetBinLabel(i+1,cutflowList[i].Data());
	}

	//for (long i = 0; i < 100000; i++)
	for (long i = 0; i < tree1->GetEntries(); i++)
	{
		tree1->GetEntry(i);

		if(NlepBL !=2) continue;
		if(nSigLep !=2) continue;
		double weight = totweight*lumiScaling;
		hCutflow->Fill("=2BaseLep && =2SigLep && SS",1);
		hCutflow->Fill("=2BaseLep && =2SigLep && SS,w",weight);

		if(isTruthLep1 != 1) continue;
		hCutflow->Fill("isTruthLep1",1);
		hCutflow->Fill("isTruthLep1,w",weight);

		if(isTruthLep2 != 1) continue;
		hCutflow->Fill("isTruthLep2",1);
		hCutflow->Fill("isTruthLep2,w",weight);

		//int SR = 1;
		int SR = 2;
		if(SR == 1)
		{
			if(nJets20 !=1) continue;
			hCutflow->Fill("nJet",1);
			hCutflow->Fill("nJet,w",weight);
		}
		else if(SR == 2)
		{
			if(nJets20 !=2 && nJets20 !=3) continue;
			hCutflow->Fill("nJet",1);
			hCutflow->Fill("nJet,w",weight);
		}

		if(nBJets20 !=0) continue;
		hCutflow->Fill("nBJets20",1);
		hCutflow->Fill("nBJets20,w",weight);

		if(Pt_l < 25000) continue;
		hCutflow->Fill("pt1",1);
		hCutflow->Fill("pt1,w",weight);

		if(Pt_subl < 25000) continue;
		hCutflow->Fill("pt2",1);
		hCutflow->Fill("pt2,w",weight);

		if(SR == 1)
		{
			if(fabs(DeltaEtaLep)>1.5) continue;
			hCutflow->Fill("DeltaEtaLep",1);
			hCutflow->Fill("DeltaEtaLep,w",weight);

			if(met<100000) continue;
			hCutflow->Fill("met",1);
			hCutflow->Fill("met,w",weight);

			if(mt<140000) continue;
			hCutflow->Fill("mt",1);
			hCutflow->Fill("mt,w",weight);

			if(meff<260000) continue;
			hCutflow->Fill("meff",1);
			hCutflow->Fill("meff,w",weight);

			if(mljj_comb>=180000) continue;
			hCutflow->Fill("mljj",1);
			hCutflow->Fill("mljj,w",weight);

			if(MT2<80000) continue;
			hCutflow->Fill("MT2",1);
			hCutflow->Fill("MT2,w",weight);
		}
		else if(SR == 2)
		{
			hCutflow->Fill("DeltaEtaLep",1);
			hCutflow->Fill("DeltaEtaLep,w",weight);

			if(met<100000) continue;
			hCutflow->Fill("met",1);
			hCutflow->Fill("met,w",weight);

			if(mt<120000) continue;
			hCutflow->Fill("mt",1);
			hCutflow->Fill("mt,w",weight);

			if(meff<240000) continue;
			hCutflow->Fill("meff",1);
			hCutflow->Fill("meff,w",weight);

			if(mljj_comb>=130000) continue;
			hCutflow->Fill("mljj",1);
			hCutflow->Fill("mljj,w",weight);

			if(MT2<70000) continue;
			hCutflow->Fill("MT2",1);
			hCutflow->Fill("MT2,w",weight);
		}

		//cout<<"MCID: "<<MCId<<endl;
		//cout<<"isTruthLep1: "<<isTruthLep1<<endl;
		//cout<<"isTruthLep2: "<<isTruthLep2<<endl;
	}

	delete tree1;

	for(unsigned int j=1;j<=cutflowList.size();j++)
	{
		cout<<hCutflow->GetXaxis()->GetBinLabel(j)<<": ";
		if(j%2 == 1) cout<<int(hCutflow->GetBinContent(j))<<endl;
		if(j%2 == 0) cout<<std::fixed<<std::setprecision(6)<<hCutflow->GetBinContent(j)<<endl;
	}
	delete hCutflow;
}

