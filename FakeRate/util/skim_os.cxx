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
const int N_MAX_FILES    = 20;

// global settings
bool gDEBUG      = true;
TString gPWD     = "$ROOTCOREBIN/..";
bool gIsDBReady  = false;

// global variables
SUSY::CrossSectionDB* xsecDB;
TFile* outFile;
TTree* outTree;
stringstream monitor;


int    g_pid1, g_pid2, g_charge1, g_charge2; 
double g_pt1, g_pt2, g_eta1, g_eta2;
double g_mt1, g_mt2, g_mt;
double g_mlj, g_detall, g_met_rel, g_meff;
double g_met, g_mttwo, g_mtmax;
int    g_njet;
int    g_source; //0-ee; 1-emu; 2-mumu
int    g_type;   //0-tt; 1-lt; 2-tl; 3-ll


///////// methods ///////////
TChain* loadData(TString fileList, TString prePath, bool isMC);
//bool passMuonCR(susyEvts* tree);
//bool passElectronCR(susyEvts* tree);
bool passPreSelection(susyEvts* tree);
bool sigRate(susyEvts* tree, bool isMC);
void initialize();
void finalize();
void calDivideErr(const double a, const double da, const double b, const double db, double &s, double &ds);
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

}

void finalize()
{
	TString output = gPWD + "/FakeRate/run/output/skim/summary.table";
	ofstream out(gSystem->ExpandPathName(output.Data()));
	out << monitor.str();


	out.close();	
}

int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		cerr << "Wrong number of arguments. " << argc - 1 << " provided, but 1 is required." << endl;
		cerr << "Command format : " << argv[0] << " <data_type>" << endl;
		cerr << "data_type: data or mc " << endl;
		return -1;
	}
	initialize();

	//TString pre_path  = "/srv/SUSY/ntuple/AnalysisBase-02-04-31/";
	TString pre_path = "";
	TString data_type = TString(argv[1]); 
	TString output    = gPWD + "/FakeRate/run/output/skim/" + data_type + "_os_skim.root";

	outFile = new TFile(output, "RECREATE");
	outTree = new TTree("signal", "signal");

	outTree->Branch("pt1",      &g_pt1,      "pt1/D");
	outTree->Branch("pt2",      &g_pt2,      "pt2/D");
	outTree->Branch("mt1",      &g_mt1,      "mt1/D");
	outTree->Branch("mt2",      &g_mt2,      "mt2/D");
	outTree->Branch("mt",       &g_mt,       "mt/D");
	outTree->Branch("eta1",     &g_eta1,     "eta1/D");
	outTree->Branch("eta2",     &g_eta2,     "eta2/D");
	outTree->Branch("detall",   &g_detall,   "detall/D");
	outTree->Branch("mttwo",    &g_mttwo,    "mttwo/D");
	outTree->Branch("mtmax",    &g_mtmax,    "mtmax/D");
	outTree->Branch("meff",     &g_meff,     "meff/D");
	outTree->Branch("met",      &g_met,      "met/D");
	outTree->Branch("mlj",      &g_mlj,      "mlj/D");
	//outTree->Branch("met_rel",  &g_met_rel,  "met_rel/D");

	outTree->Branch("source",   &g_source,   "source/I");
	outTree->Branch("charge1",  &g_charge1,  "charge1/I");
	outTree->Branch("charge2",  &g_charge2,  "charge2/I");
	outTree->Branch("pid1",     &g_pid1,     "pid1/I");
	outTree->Branch("pid2",     &g_pid2,     "pid2/I");
	outTree->Branch("type",     &g_type,     "type/I");
	outTree->Branch("njet",     &g_njet,     "nget/I");

	bool isMC;
	vector<TString> files;
	files.clear();
	if (data_type == "data")
	{
		//files.push_back("../share/inFileList-data.txt");
		files.push_back("../share/data.txt");
		
		isMC = false;
	}
	else if (data_type == "mc")
	{
		files.push_back( "../share/inFileList-ttbar.txt" );
		files.push_back( "../share/inFileList-DY.txt" );
		files.push_back( "../share/inFileList-ttV.txt" );
		files.push_back( "../share/inFileList-ZPowheg.txt" );
		files.push_back( "../share/inFileList-Wjets.txt" );
		files.push_back( "../share/inFileList-VV.txt" );
		files.push_back( "../share/inFileList-Vgamma.txt" );
		files.push_back( "../share/inFileList-SingleTop.txt" );
		isMC = true;
	}
	else 
	{
		cerr << "<data_type> should be data or mc provided." << endl; 
		return -1;	
	}

	susyEvts *evts[N_MAX_FILES];

	for (unsigned int i = 0; i < files.size(); i++)
	{
		TChain *tc = loadData(files[i], pre_path, isMC);
		evts[i] = new susyEvts(tc);

		cout << "In " << files[i] << " control samples, there are " ;
		cout << evts[i]->tree1->GetEntries() << " events" << endl;
		monitor << "After " << files[i] << " samples: " << endl;
		sigRate(evts[i], isMC);
		delete evts[i];
	}


	finalize();

	outFile->Write();
	outFile->Close();
	cout << "Finished !" << endl;

	return 0;
}


TChain* loadData(TString fileList, TString prePath, bool isMC)
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
	for (auto f : allFiles)
	{
		TRegexp expr  = ".root.?[0-9]?$";
		TString fname;
		if (isMC) 
		{
			fname  = TString(f(0, f.Index(expr)));
			fname += "_WEIGHTED.root";
		}
		else 
		{
			fname = f;
		}
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

bool sigRate(susyEvts* tree, bool isMC)
{
	long nEntries = tree->tree1->GetEntries();
	if (gDEBUG) cout << "Begin calculating signal rates" << endl;
	//for (long i = 0; i < 20000; i++)
	for (long i = 0; i < nEntries; i++)
	{
		tree->GetEntry(i); 

		double w = 1.0;
		if (isMC)
		{
			w *= ((TChain*)(tree->tree1))->GetTree()->GetWeight();
			w *= tree->evt.weight;
			w *= tree->evt.pwt;
			w *= tree->evt.ElSF;
			w *= tree->evt.MuSF;
		}

		if (passPreSelection(tree))
		{ 
			g_pt1 = tree->leps[0].pt;
			g_pt2 = tree->leps[1].pt;
			g_eta1 = tree->leps[0].eta;
			g_eta2 = tree->leps[1].eta;	
			g_njet = tree->jets.size();
			g_detall = fabs(g_eta1 - g_eta2);
			g_mttwo = tree->sig.mT2;
			g_met = tree->sig.Met;
			//g_met_rel = tree->sig.MetRel;
			g_mlj = tree->sig.mlj;
			g_pid1 = (int)tree->leps[0].ID/1000;
			g_pid2 = (int)tree->leps[1].ID/1000;
			g_meff = tree->sig.HT + tree->sig.Met; 

			double mt1 = tree->leps[0].mT;
			double mt2 = tree->leps[1].mT;
			g_mt1 = mt1;
			g_mt2 = mt2;
			if (g_pt1 > g_pt2)
			{
				g_mt = mt1;
			}
			else 
			{
				g_mt = mt2;
			}
			if (mt1 > mt2)
			{
				g_mtmax = mt1;
			}
			else 
			{
				g_mtmax = mt2;
			}
			// requirement
			if (g_pt1 < 25) continue;
			if (g_pt2 < 25) continue;
			//if (g_detall > 1.5) continue;
			//if (g_meff < 200) continue;

			int product = fabs(g_pid1 * g_pid2);
			if (product == 121)
			{
				g_source = 0;
			}
			else if (product == 143)
			{
				g_source = 1;
			}
			else if (product == 169)
			{
				g_source = 2;
			}
			else 
			{
				g_source = -1;
			}

			bool pass1 = (tree->leps[0].lFlag & IS_SIGNAL);
			bool pass2 = (tree->leps[1].lFlag & IS_SIGNAL);
			if (pass1 && pass2)
			{
				g_type = 0;
			}
			else if (pass1 && !pass2)
			{
				g_type = 1;
			}
			else if (!pass1 && pass2)
			{
				g_type = 2;
			}
			else 
			{
				g_type = 3;
			}
			outTree->Fill();
		}
	}

	return true;
}

bool passPreSelection(susyEvts* tree)
{
	bool pass = true;
	// trigger
	//pass *= tree->sig.trigCode > 0;
	//pass *= ( tree->sig.trigCode & tree->sig.trigMask ) != 0;
	// two leptons
	pass *= tree->leps.size() == 2;
	// same sign
	int cp = int(tree->leps[0].ID/1000) * int(tree->leps[1].ID/1000);
	pass *= (cp	== -121);

	// no b-jets
	bool fJet = false;
	for (unsigned int i = 0; i < tree->jets.size(); i++) 
	{
		fJet = fJet || (tree->jets[i].jFlag & JT_BJET);
	}
	pass *= !fJet;

	// 1 < njets < 4
	pass *= (tree->jets.size() > 0 && tree->jets.size() < 4);

	return pass;
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

