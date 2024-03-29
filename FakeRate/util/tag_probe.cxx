// C++
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <cstdlib>
// ROOT
#include "TSystem.h"
#include "TChain.h"
#include "TString.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TH1D.h"
#include "TRegexp.h"
#include "TH2D.h"
#include "TDirectory.h"
// ATLAS
#include "SUSYTools/SUSYCrossSection.h"
#include "MCTruthClassifier/MCTruthClassifierDefs.h"
// My packages
#include "FakeRate/susyEvts.h"

#define PRINT(x) {std::cout << #x << " = " << x << std::endl;}
using namespace std;

// constants
const TString CHAIN_NAME = "evt2l";
const double M_Z         = 91.1876;
const double M_Z_WIDTH   = 10;
const double LUMI        = 36074.56;
const int N_EL_SOURCE    = 6;
const int N_MU_SOURCE    = 5;
const int N_PROC         = 8;

enum STUDY
{
	SS = 0,
	THREEL,
};

enum LEP_TYPE 
{
	ELEC,
	MUON,
};

const double ETA_EL[] = {0, 0.8, 1.37, 1.52, 2.47};
const double ETA_MU[] = {0, 0.6, 1.2, 1.8, 2.4};
const double PT[]  = {25, 35, 45, 55, 65, 75, 85, 95, 7000};
const unsigned int NETA_EL = sizeof(ETA_EL) / sizeof(ETA_EL[0]) - 1;
const unsigned int NETA_MU = sizeof(ETA_MU) / sizeof(ETA_MU[0]) - 1;
const unsigned int NPT     = sizeof(PT)  / sizeof(PT[0]) - 1;

struct Histo 
{
	Histo(TString name, LEP_TYPE e) 
	{
		_e = e;
		if (e == LEP_TYPE::ELEC)
		{
			hTight = new TH2D(name+"_hTight", ";p_{T} [GeV];|#eta|", NPT, PT, NETA_EL, ETA_EL);
			hLoose = new TH2D(name+"_hLoose", ";p_{T} [GeV];|#eta|", NPT, PT, NETA_EL, ETA_EL);
			hEff   = new TH2D(name+"_hEff",   ";p_{T} [GeV];|#eta|", NPT, PT, NETA_EL, ETA_EL);
		}
		else 
		{
			hTight = new TH2D(name+"_hTight", ";p_{T} [GeV];|#eta|", NPT, PT, NETA_MU, ETA_MU);
			hLoose = new TH2D(name+"_hLoose", ";p_{T} [GeV];|#eta|", NPT, PT, NETA_MU, ETA_MU);
			hEff   = new TH2D(name+"_hEff",   ";p_{T} [GeV];|#eta|", NPT, PT, NETA_MU, ETA_MU);
		}
		hTight->Sumw2();
		hLoose->Sumw2();
		hEff->Sumw2();
	}
	void calEff()
	{
		hEff->Add(hTight); 
		hEff->Divide(hLoose); 
		
		unsigned int NETA;
		if (_e == LEP_TYPE::ELEC) NETA = NETA_EL;
		else NETA = NETA_MU;
		
		cout<<"result for real eff:"<<endl;
		for(unsigned int i=1;i<=NPT+1;i++)
		{
			for(unsigned int j=1;j<=NETA;j++)
			{
				cout<<hEff->GetBinContent(i,j)<<", "<<hEff->GetBinError(i,j)<<endl;
			}
			cout<<endl;
		}
		
	}
	TH2D *hLoose;
	TH2D *hTight;
	TH2D *hEff;
	LEP_TYPE _e;
};

Histo *real_mu;
Histo *real_el;

// global settings
bool gDEBUG      = true;
STUDY gStudy     = SS;
TString gPWD     = "$ROOTCOREBIN/..";
bool gIsDBReady  = false;

// global variables
SUSY::CrossSectionDB* xsecDB;
TFile* outFile;
stringstream monitor;

///////// methods ///////////
TChain* loadData(TString fileList, TString prePath = "./");
bool ptEtaRequirement(double pt, double eta, int ID);
bool passRatesCR(susyEvts* tree, int& flag);
void fillHist(susyEvts* tree, int index);
void sigRate(susyEvts* tree);
void initialize();
void finalize();
void calDivideErr(const double a, const double da, const double b, const double db, double &s, double &ds);
char* itoa(int num, char* str,int radix);
/////////////////////////////

bool initializeDB()
{
	xsecDB = new SUSY::CrossSectionDB(gSystem->ExpandPathName("$ROOTCOREBIN/data/SUSYTools/mc15_13TeV/"));
	if (!xsecDB)
	{
		cout << "CrossSectionDB could not be loaded. " << endl;
		return false;
	}
	else 
	{
		gIsDBReady = true;
		cout << "CrossSectionDB is loaded successfullly." << endl;
		return true;
	}
}

void initialize()
{
	cout << "Initializing ... " << endl;
	monitor.clear();

	initializeDB();

	TString output = gPWD + "/FakeRate/run/output/tag_probe/fake_rate.root";

	outFile = new TFile(output, "RECREATE");
	
	real_el = new Histo("El", LEP_TYPE::ELEC);
	real_mu = new Histo("Mu", LEP_TYPE::MUON);

}

void finalize()
{
	TString output = gPWD + "/FakeRate/run/output/tag_probe/summary.table";
	ofstream out(gSystem->ExpandPathName(output.Data()));
	out << monitor.str();

	// check the number of el, mu components
	out.close();	

	real_el->calEff();
	real_mu->calEff();

	outFile->Write();
	outFile->Close();

}

int main()
{
	initialize();

	/*
	//TString pre_path = "/srv/SUSY/ntuple/AnalysisBase-02-04-31/";
	TString pre_path = "/eos/user/c/clo/ntuple/AnalysisBase-02-04-31-fake/";
	TString files = "../share/data.txt";
	TChain *tc = loadData(files, pre_path);
	*/

	TString pre_path = "/eos/user/c/clo/ntuple/AnalysisBase-02-04-31-6ecc6eb7/user.clo.v13.5.";
	TString path = pre_path + "data_myOutput.root/user.clo.data*.*.myOutput.root*";
	TChain *tc = new TChain( CHAIN_NAME );
	tc->Add(path);

	susyEvts *evts = new susyEvts(tc);
	cout << "There are " << evts->tree1->GetEntries() << " events" << endl;
	sigRate(evts);
	delete evts;

	finalize();

	cout << "Finished !" << endl;

	return 0;
}


TChain* loadData(TString fileList, TString prePath)
{
	TChain* tc = new TChain( CHAIN_NAME );

	ifstream in(fileList.Data());
	vector<TString> allFiles;
	allFiles.clear();

	for (string line; getline( in, line ); )
	{
		if (gDEBUG) cout << line << endl;
		if (line[0] == '#') continue;
		allFiles.push_back(prePath + line);
	}

	cout << "Loading NTuples from " << fileList << endl;
	for (auto fname : allFiles)
	{
		if (gSystem->AccessPathName(fname, kFileExists)) //if file exists, return false 
		{
			cout << ">> File: '" << fname << "' DO NOT exist!" << endl;
			continue;
		}
		if (tc->Add(fname))
		{
			cout << ">> File: '" << fname << "' :: TTree '" << tc->GetName() << "' has been loaded" << endl;
		}
		else
		{
			cout << ">> File: '" << fname << "' :: TTree '" << tc->GetName() << "' cannot be loaded" << endl;
		}
	}
	return tc;
}

void sigRate(susyEvts* tree)
{
	unsigned int nEntries = tree->tree1->GetEntries();
	if (gDEBUG) cout << "Begin calculating signal rates" << endl;
	//for (unsigned int i = 0; i < 10000; i++)
	for (unsigned int i = 0; i < nEntries; i++)
	{
		tree->GetEntry(i); 
		
		int flag = 0;
		if (passRatesCR(tree, flag))
		{
			if(flag & 1) fillHist(tree,1);
			if(flag & 2) fillHist(tree,0);
		}
	}
}

bool passRatesCR(susyEvts* tree, int &flag)
{
	// two leptons
	if(tree->leps.size() != 2) return false;
	
	// SFOS
	if(int(tree->leps[0].ID /1000) + int(tree->leps[1].ID /1000) != 0) return false;
	
	// pt and eta requirement 
	if(!ptEtaRequirement(tree->leps[0].pt, tree->leps[0].eta, tree->leps[0].ID)) return false;
	if(!ptEtaRequirement(tree->leps[1].pt, tree->leps[1].eta, tree->leps[1].ID)) return false;
	
	// in Z-mass window
	if(fabs(tree->l12.m - M_Z) > M_Z_WIDTH) return false;
	
	/*
	// no more than three jets
	if(tree->jets.size() > 3) return false;
	
	// no b-jets
	for (unsigned int i = 0; i < tree->jets.size(); i++) 
	{
		if(tree->jets[i].jFlag & JT_BJET) return false;
	}
	*/
	
	// tag with signal requirement
	if(tree->leps[0].pt > 25 && (tree->leps[0].lFlag & IS_SIGNAL))
	{
		flag |= 1;
	}
	
	if(tree->leps[1].pt > 25 && (tree->leps[1].lFlag & IS_SIGNAL))
	{
		flag |= 2;
	}
	
	if(flag == 0) return false;
	else return true;
}

void fillHist(susyEvts* tree, int index)
{
	//For probe lepton
	int id     = int(abs(tree->leps[index].ID) / 1000);
	double pt  = tree->leps[index].pt;
	double eta = fabs(tree->leps[index].eta);
	bool isTight = tree->leps[index].lFlag & IS_SIGNAL;

	if (id == 11) //electron 
	{
		real_el->hLoose->Fill(pt, eta);
		if (isTight)
		{
			real_el->hTight->Fill(pt, eta);
		}
	}
	else if (id == 13) //muon
	{
		real_mu->hLoose->Fill(pt, eta);
		if (isTight)
		{
			real_mu->hTight->Fill(pt, eta);
		}
	}
}

void calDivideErr(const double a, const double da, const double b, const double db, double &s, double &ds)
{
	//a -- numerator
	//b -- denominator
	double EPS = 1e-6;
	if (fabs(a) > EPS && fabs(b) > EPS)
	{
		double k = da/a;
		double j = db/b;
		s  = a/b;
		ds = s * sqrt(k*k + j*j);
	}
	else 
	{
		s  = 0;
		ds = 0;
	}
}

char* itoa(int num, char* str,int radix)
{
	char index[] = "0123456789ABCDEF";
	unsigned unum;
	int i = 0, j, k;
	if (radix == 10 && num < 0)
	{
		unum     = (unsigned)-num;
		str[i++] = '-';
	}
	else 
	{
		unum     = (unsigned)num;
	}
	do
	{
		str[i++] = index[unum%(unsigned)radix];
		unum    /=radix;
	}
	while(unum);
	str[i] = '\0';
	if (str[0] == '-')
		k = 1;
	else 
		k = 0;
	char temp;
	for (j = k; j <= (i-1)/2; j++)
	{
		temp         = str[j];
		str[j]       = str[i-1+k-j];
		str[i-1+k-j] = temp;
	}
	return str;
}

bool ptEtaRequirement(double pt, double eta, int ID)
{
	//if (pt < PT[0] || pt > PT[NPT])
	if (pt < PT[0])
	{
		return false;
	}
	if (int(abs(ID)/1000) == 11)
	{
		if (fabs(eta) < ETA_EL[0] || fabs(eta) > ETA_EL[NETA_EL])
		{
			return false;
		}
		else 
		{
			return true;
		}
	}
	else if (int(abs(ID)/1000) == 13)
	{
		if (fabs(eta) < ETA_MU[0] || fabs(eta) > ETA_MU[NETA_MU])
		{
			return false;
		}
		else 
		{
			return true;
		}
	}
	else
	{
		return false;
	}
}
