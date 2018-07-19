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

bool ptEtaRequirement(double pt, double eta, int ID);
double getEff(double pt,double eta,int ID,TH2D* heEff,TH2D* huEff);

int main()
{
	const TString CHAIN_NAME = "evt2l";
	
	//read old tree
	TString path = "/srv/SUSY/ntuple/AnalysisBase-02-04-39-4171b36f/user.clo.v12.1.data_myOutput.root/";
	path += "*.root*";
	TChain* tree1 = new TChain(CHAIN_NAME);
	tree1->Add(path.Data());
	cout << "There are " << tree1->GetEntries() << " events" << endl;
	susyEvts* evt1 = new susyEvts(tree1);

	//create new tree
	TFile* output = new TFile("output/gen_ntuple/fake.root","RECREATE");
	susyEvts* evt2 = new susyEvts();
	evt2->makeTree(CHAIN_NAME);
	evt2->tree2->SetDirectory(output);

	//get 2D hist for real efff
	TFile* file_real = new TFile("output/tag_probe/fake_rate.root","READ");
	TH2D *hRealeEff = (TH2D*) file_real->Get("El_hEff");
	TH2D *hRealuEff = (TH2D*) file_real->Get("Mu_hEff");

	//get 2D hist for fake efff
	TFile* file_fake = new TFile("../scripts/getFakeEffs/fake_eff.root","READ");
	TH2D *hFakeeEff = (TH2D*) file_fake->Get("El_hEff");
	TH2D *hFakeuEff = (TH2D*) file_fake->Get("Mu_hEff");

	//for (long i = 0; i < 100000; i++)
	for (long i = 0; i < tree1->GetEntries(); i++)
	{
		evt1->GetEntry(i);

		// two leptons
		if(evt1->leps.size() != 2) continue;
		
		// SS
		int ID1 = int(evt1->leps[0].ID/1000);
		int ID2 = int(evt1->leps[1].ID/1000);
		int product = ID1 * ID2;
		if(product < 0) continue;

		//pt and eta requirement
		ID1 = abs(ID1);
		ID2 = abs(ID2);
		double pt1  = evt1->leps[0].pt;
		double pt2  = evt1->leps[1].pt;
		double eta1 = fabs(evt1->leps[0].eta);
		double eta2 = fabs(evt1->leps[1].eta);
		if(!ptEtaRequirement(pt1,eta1,ID1)) continue;
		if(!ptEtaRequirement(pt2,eta2,ID2)) continue;

		//get eff
		double realEff1 = getEff(pt1,eta1,ID1,hRealeEff,hRealuEff);
		double realEff2 = getEff(pt2,eta2,ID2,hRealeEff,hRealuEff);
		double fakeEff1 = getEff(pt1,eta1,ID1,hFakeeEff,hFakeuEff);
		double fakeEff2 = getEff(pt2,eta2,ID2,hFakeeEff,hFakeuEff);

		//calculate fake weight
		bool isTight1 = evt1->leps[0].lFlag & IS_SIGNAL;
		bool isTight2 = evt1->leps[1].lFlag & IS_SIGNAL;

		//For the formula see:
		//http://live.sympy.org/?evaluate=e1%2Ce2%2Cf1%2Cf2%20%3D%20symbols%28%22e1%20e2%20f1%20f2%22%29%0A%23--%0Ae1%0A%23--%0Am%20%3D%20Matrix%28[[e1%2Ce2]%2C[f1%2Cf2]]%29%0A%23--%0Am%0A%23--%0Am.inv%28Abs%28%29%0A%23--%0Am%0A%23--%0Am.inv%0A%23--%0Am.inv%28%29%0A%23--%0Am.inv%28%29*m%0A%23--%0Asimplify%28m.inv%28%29*m%29%0A%23--%0Asimplify%28m.inv%28%29%29%0A%23--%0Am%20%3D%20Matrix%28[[e1*e2%2Ce1*f2%2Cf1*e2%2Cf1*f2]%2C[e1*%281-e2%29%2Ce1*%281-f2%29%2Cf1*%281-e2%29%2Cf1*%281-f2%29]%2C[%281-e1%29*e2%2C%281-e1%29*f2%2C%281-f1%29*e2%2C%281-f1%29*f2]%2C[%281-e1%29*%281-e2%29%2C%281-e1%29*%281-f2%29%2C%281-f1%29*%281-e2%29%2C%281-f1%29*%281-f2%29]]%29%0A%23--%0Am%0A%23--%0Am2%20%3D%20m.inv%28%29%0A%23--%0Asimplify%28m2%29%0A%23--%0Am2s%20%3D%20simplify%28m2%29%0A%23--%0A
		double d = realEff1*realEff2 + fakeEff1*fakeEff2 - realEff1*fakeEff2 - fakeEff1*realEff2;
		d = 1/d;
		double Nrf, Nfr, Nff;

		if      ( isTight1 &&  isTight2)
		{
			//Ntt
			Nrf = d * ( -realEff2*fakeEff1 + realEff2 + fakeEff1 -1. );
			Nfr = d * ( -realEff1*fakeEff2 + realEff1 + fakeEff2 -1. );
			Nff = d * (  realEff1*realEff2 - realEff1 - realEff2 +1. );
		}
		else if ( isTight1 && !isTight2)
		{
			//Ntl
			Nrf = d * ( -realEff2*(fakeEff1-1.) );
			Nfr = d * ( -fakeEff2*(realEff1-1.) );
			Nff = d * (  realEff2*(realEff1-1.) );
		}
		else if (!isTight1 &&  isTight2)
		{
			//Nlt
			Nrf = d * ( -fakeEff1*(realEff2-1.) );
			Nfr = d * ( -realEff1*(fakeEff2-1.) );
			Nff = d * (  realEff1*(realEff2-1.) );
		}
		else if (!isTight1 && !isTight2)
		{
			//Nll
			Nrf = d * ( -realEff2*fakeEff1 );
			Nfr = d * ( -realEff1*fakeEff2 );
			Nff = d * (  realEff1*realEff2 );
		}
		else
		{
			cout<<"Impossible case"<<endl;;
			Nrf = 0.;
			Nfr = 0.;
			Nff = 0.;
			continue;
		}
		
		//copy
		evt2->evt = evt1->evt;
		evt2->leps = evt1->leps;
		evt2->l12 = evt1->l12;
		evt2->jets = evt1->jets;
		evt2->truths = evt1->truths;
		evt2->sig = evt1->sig;

		//set fake weight
		evt2->evt.fLwt = Nrf * realEff1 * fakeEff2 + Nfr * fakeEff1 * realEff2 + Nff * fakeEff1 * fakeEff2 ;
		
		//fill event
		evt2->fill();
	}

	output->cd();	
	evt2->tree2->Write(CHAIN_NAME);

	delete evt1;
	delete tree1;
	delete evt2;
	delete output;
	delete file_real;
	delete file_fake;
}

bool ptEtaRequirement(double pt, double eta, int ID)
{
	if(pt < 25) return false;
	if(ID == 11)
	{
		if(eta>=2.47) return false;
	}
	else if(ID == 13)
	{
		if(eta>=2.4) return false;
	}
	return true;
}

double getEff(double pt,double eta,int ID,TH2D* heEff,TH2D* huEff)
{
	TH2D* hEff;
	if(ID == 11) hEff = heEff;
	else if(ID == 13) hEff = huEff;
	else {hEff = 0; return 0;}

	int binx = hEff->GetXaxis()->FindBin(pt);
	if(binx == 0)
	{
		cout<<"pt underflow"<<endl;
		return 0;
	}

	int biny = hEff->GetYaxis()->FindBin(eta);
	if(biny == 0)
	{
		cout<<"eta underflow"<<endl;
		return 0;
	}
	else if(biny == hEff->GetYaxis()->GetXbins()->GetSize())
	{
		cout<<"eta overflow"<<endl;
		return 0;
	}

	return hEff->GetBinContent(binx,biny);
}
