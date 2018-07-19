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

const double ETA_EL[] = {0, 1.37, 1.52, 2.47};
const double ETA_MU[] = {0, 1.37, 1.52, 2.4};
const double PT_EL[]  = {25, 35, 45, 120, 200};
const double PT_MU[]  = {25, 30, 45, 120, 200};
const unsigned int NETA_EL = sizeof(ETA_EL) / sizeof(ETA_EL[0]) - 1;
const unsigned int NETA_MU = sizeof(ETA_MU) / sizeof(ETA_MU[0]) - 1;
const unsigned int NPT_EL  = sizeof(PT_EL)  / sizeof(PT_EL[0])  - 1;
const unsigned int NPT_MU  = sizeof(PT_MU)  / sizeof(PT_MU[0])  - 1;

struct Histo 
{
	Histo(TString name, LEP_TYPE e) 
	{
		_e = e;
		if (e == LEP_TYPE::ELEC)
		{
			hTight = new TH2D(name+"_hTight", ";p_{T} [GeV];|#eta|", NPT_EL, PT_EL, NETA_EL, ETA_EL);
			hLoose = new TH2D(name+"_hLoose", ";p_{T} [GeV];|#eta|", NPT_EL, PT_EL, NETA_EL, ETA_EL);
		}
		else 
		{
			hTight = new TH2D(name+"_hTight", ";p_{T} [GeV];|#eta|", NPT_MU, PT_MU, NETA_MU, ETA_MU);
			hLoose = new TH2D(name+"_hLoose", ";p_{T} [GeV];|#eta|", NPT_MU, PT_MU, NETA_MU, ETA_MU);
		}
		hTight->Sumw2();
		hLoose->Sumw2();
	}
	TH2D *hLoose;
	TH2D *hTight;
	LEP_TYPE _e;
};

Histo *real_mu;
Histo *real_el;
Histo *prompt_mu;
Histo *prompt_el;

// global settings
bool gDEBUG      = true;
STUDY gStudy     = SS;
TString gPWD     = "$ROOTCOREBIN/..";
bool gIsDBReady  = false;
bool DoCutflow = false;

// global variables
SUSY::CrossSectionDB* xsecDB;
TFile* outFile;
stringstream monitor;

///////// methods ///////////
std::vector<TChain*> loadData(TString fileList, TString prePath, bool isMC);
bool ptEtaRequirement(double pt, double eta, LEP_TYPE e);
bool passMuonCR(susyEvts* tree);
bool passElectronCR(susyEvts* tree);
bool sigRate(susyEvts* tree, bool isMC, double treeWeight, TH1D* hCutflow);
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

}

void finalize()
{
	TString output = gPWD + "/FakeRate/run/output/fr_mxm/summary.table";
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
	TString pre_path  = "/srv/SUSY/ntuple/AnalysisBase-02-04-39-4171b36f/";
	TString data_type = TString(argv[1]); 
	TString output    = gPWD + "/FakeRate/run/output/fr_mxm/" + data_type + "_fr_mxm.root";

	outFile = new TFile(output, "RECREATE");
	real_el = new Histo("El", LEP_TYPE::ELEC);
	real_mu = new Histo("Mu", LEP_TYPE::MUON);
	prompt_el = new Histo("El_prompt", LEP_TYPE::ELEC);
	prompt_mu = new Histo("Mu_prompt", LEP_TYPE::MUON);

	bool isMC;
	vector<TString> files;
	files.clear();
	//TString path = "../share/";
	TString path = "../share_new/";
	if (data_type == "data")
	{
		files.push_back( path + "data.txt");
		isMC = false;
	}
	else if (data_type == "mc")
	{
		//files.push_back( path + "inFileList-ttbar.txt" );
		//files.push_back( path + "inFileList-DY.txt" );
		//files.push_back( path + "inFileList-ZPowheg.txt" );
		//files.push_back( path + "inFileList-Wjets.txt" );
		//files.push_back( path + "inFileList-Vgamma.txt" );
		//files.push_back( path + "inFileList-SingleTop.txt" );
		files.push_back( path + "inFileList-higgs.txt" );
		files.push_back( path + "inFileList-VV.txt" );
		files.push_back( path + "inFileList-VVV.txt" );
		files.push_back( path + "inFileList-ttV.txt" );
		files.push_back( path + "inFileList-multitop.txt" );
		isMC = true;
	}
	else 
	{
		cerr << "<data_type> should be data or mc provided." << endl; 
		return -1;	
	}
	
	//hist for cutflow
	std::vector<TString> cutflowList;
	cutflowList.push_back("ntuple");
	cutflowList.push_back("=2BaseLep");
	
	TH1D* hCutflow = new TH1D("cutflow", "cutflow", cutflowList.size(), 0, cutflowList.size());
	for(unsigned int i=0;i<cutflowList.size();i++)    
	{
		hCutflow->GetXaxis()->SetBinLabel(i+1,cutflowList[i].Data());
	}
	
	for (unsigned int i = 0; i < files.size(); i++)
	{
		std::vector<TChain*> tc = loadData(files[i], pre_path, isMC);
		for (unsigned int j = 0; j < tc.size(); j++)
		{
			susyEvts *evts = new susyEvts(tc[j]);
			cout << "There are " << evts->tree1->GetEntries() << " events" << endl;

			double treeWeight = 1;
			if (isMC)
			{ 
				double sumW = 0;
				TObjArray* fileList = tc[j]->GetListOfFiles();
				for(int k=0;k<fileList->GetEntries();k++)
				{
					TFile *file = TFile::Open(fileList->At(k)->GetTitle());
					TH1D *hCutFlow = (TH1D*) file->Get("hCutFlow");
					sumW += hCutFlow->GetBinContent(2);
					delete file;
				}
                
				TString f = fileList->At(0)->GetTitle();
				int sampleID = TString(f(f.Index("TeV")+4,6)).Atoi();
				double xSecxEff = xsecDB->xsectTimesEff(sampleID);
				treeWeight = xSecxEff / sumW * LUMI;
				cout<<"treeWeight: "<<treeWeight<<endl;
			}

			sigRate(evts, isMC, treeWeight, hCutflow);
			delete evts;
		}
	}
	
	//print cutflow
	if(DoCutflow)
	{
		for(unsigned int j=1;j<=cutflowList.size();j++)
		{
			cout<<hCutflow->GetXaxis()->GetBinLabel(j)<<": ";
			cout<<int(hCutflow->GetBinContent(j))<<endl;
		}
	}
	
	finalize();
	outFile->Write();
	outFile->Close();
	cout << "Finished !" << endl;
	
	return 0;
}


std::vector<TChain*> loadData(TString fileList, TString prePath, bool isMC)
{
	std::vector<TChain*> tc;

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

	int sampleID1 = 0;
	TChain* tc2 = new TChain( CHAIN_NAME );
	for (auto f : allFiles)
	{
		if (isMC) 
		{
			int sampleID2 = TString(f(f.Index("TeV")+4,6)).Atoi();
			if(sampleID2 != sampleID1 && sampleID1 != 0)
			{
				tc.push_back(tc2);
				tc2 = new TChain( CHAIN_NAME );
			}
			sampleID1 = sampleID2;
		}

		if (gSystem->AccessPathName(f, kFileExists)) //if file exists, return false 
		{
			cout << ">> File: '" << f << "' DO NOT exist!" << endl;
			continue;
		}
		if (tc2->Add(f))
		{
			cout << ">> File: '" << f << "' :: TTree '" << tc2->GetName() << "' has been loaded" << endl;
		}
		else
		{
			cout << ">> File: '" << f << "' :: TTree '" << tc2->GetName() << "' cannot be loaded" << endl;
		}
	}
	tc.push_back(tc2);
	return tc;
}

bool sigRate(susyEvts* tree, bool isMC, double treeWeight, TH1D* hCutflow)
{
	using namespace MCTruthPartClassifier;
	long nEntries = tree->tree1->GetEntries();
	if (gDEBUG) cout << "Begin calculating signal rates" << endl;
	for (long i = 0; i < nEntries; i++)
	{
		tree->GetEntry(i); 
		
		if(DoCutflow) hCutflow->Fill("ntuple",1);
		
		// trigger
		//For old ntuple
		//if(tree->sig.trigCode<=0) {cout<<"Trigger Error!!!!!"<<endl; continue;}
		
		// two baseline leptons
		if(tree->leps.size() != 2) continue;
		if(DoCutflow) hCutflow->Fill("=2BaseLep",1);
		
		int product = int(tree->leps[0].ID/1000) * int(tree->leps[1].ID/1000);
		// SS emu or SS mumu
		if(product != 143 && product != 169) continue;
		
		// more than one jets
		if(tree->jets.size() == 0) continue;
		
		// at least one b-jets
		bool fJet = false;
		for (unsigned int i = 0; i < tree->jets.size(); i++) 
		{
			fJet = fJet || (tree->jets[i].jFlag & JT_BJET);
		}
		if(!fJet) continue;
        	
		//if(tree->sig.Met + tree->sig.HT <= 200) continue;
		
		double w = 1.0;
		if (isMC)
		{
			w *= treeWeight;
			w *= tree->evt.weight;
			w *= tree->evt.pwt;
			w *= tree->evt.ElSF;
			w *= tree->evt.MuSF;
			w *= tree->evt.BtagSF;
			w *= tree->evt.JvtSF;
			w *= tree->evt.trigSF;
		}
		//electron sample	
		if (product == 143)
		{ 
			// muon tag 
			int tag_idx, probe_idx;
			int id = int(tree->leps[0].ID/1000);
			if (abs(id) == 13)
			{
				tag_idx   = 0;
				probe_idx = 1;
			}
			else 
			{
				tag_idx   = 1;
				probe_idx = 0;
			}
			if(tree->leps[tag_idx].pt <= 40) continue;
			if(!(tree->leps[tag_idx].lFlag & IS_SIGNAL)) continue;
			//tag assigned
			{
				double pt  = tree->leps[probe_idx].pt;
				double eta = fabs(tree->leps[probe_idx].eta); 
				if(!ptEtaRequirement(pt, eta, LEP_TYPE::ELEC)) continue;
				{
					bool isTight = tree->leps[probe_idx].lFlag & IS_SIGNAL;
					real_el->hLoose->Fill(pt, eta, w);
					if (isTight) real_el->hTight->Fill(pt, eta, w);
					if (isMC)
					{
						//do something here to keep prompt-lepton info
						ParticleType   type = static_cast<ParticleType>(tree->leps[probe_idx].truthType);
						if (type == IsoElectron)
						{
							prompt_el->hLoose->Fill(pt, eta, w);
							if (isTight) prompt_el->hTight->Fill(pt, eta, w);
						}
					}
				}
			}
		}
		//muon sample
		else if (product == 169)
		{ 
			// muon tag 
			for (unsigned int j = 0; j < 2; j++)
			{
				int tag_idx, probe_idx;
				if (j == 0)
				{
					tag_idx   = 0;
					probe_idx = 1;
				}
				else 
				{
					tag_idx   = 1;
					probe_idx = 0;
				}
				if(tree->leps[tag_idx].pt <= 40) continue;
				if(!(tree->leps[tag_idx].lFlag & IS_SIGNAL)) continue;
				//tag assigned
				{
					double pt  = tree->leps[probe_idx].pt;
					double eta = fabs(tree->leps[probe_idx].eta); 
					if(!ptEtaRequirement(pt, eta, LEP_TYPE::MUON)) continue;
					{
						bool isTight = tree->leps[probe_idx].lFlag & IS_SIGNAL;
						real_mu->hLoose->Fill(pt, eta, w);
						if (isTight) real_mu->hTight->Fill(pt, eta, w);
						if (isMC)
						{
							ParticleType   type = static_cast<ParticleType>(tree->leps[probe_idx].truthType);
							if (type == IsoMuon)
							{
								prompt_mu->hLoose->Fill(pt, eta, w);
								if (isTight) prompt_mu->hTight->Fill(pt, eta, w);
							}
						}
					}
				}
			}
		}
		//end of loop for entries 
	}
	
	monitor << real_el->hLoose->GetName() << "\t";
	monitor << real_el->hLoose->GetEntries() << "\t";
	monitor << real_el->hLoose->GetEffectiveEntries() << endl;
	monitor << prompt_el->hLoose->GetName() << "\t";
	monitor << prompt_el->hLoose->GetEntries() << "\t";
	monitor << prompt_el->hLoose->GetEffectiveEntries() << endl;
	monitor << real_el->hTight->GetName() << "\t";
	monitor << real_el->hTight->GetEntries() << "\t";
	monitor << real_el->hTight->GetEffectiveEntries() << endl;
	monitor << prompt_el->hTight->GetName() << "\t";
	monitor << prompt_el->hTight->GetEntries() << "\t";
	monitor << prompt_el->hTight->GetEffectiveEntries() << endl;
	monitor << real_mu->hLoose->GetName() << "\t";
	monitor << real_mu->hLoose->GetEntries() << "\t";
	monitor << real_mu->hLoose->GetEffectiveEntries() << endl;
	monitor << prompt_mu->hLoose->GetName() << "\t";
	monitor << prompt_mu->hLoose->GetEntries() << "\t";
	monitor << prompt_mu->hLoose->GetEffectiveEntries() << endl;
	monitor << real_mu->hTight->GetName() << "\t";
	monitor << real_mu->hTight->GetEntries() << "\t";
	monitor << real_mu->hTight->GetEffectiveEntries() << endl;
	monitor << prompt_mu->hTight->GetName() << "\t";
	monitor << prompt_mu->hTight->GetEntries() << "\t";
	monitor << prompt_mu->hTight->GetEffectiveEntries() << endl;

	return true;
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

bool ptEtaRequirement(double pt, double eta, LEP_TYPE e)
{
	if (e == LEP_TYPE::ELEC)
	{
		if (pt < PT_EL[0] || pt > PT_EL[NPT_EL]) return false;
		if (eta < ETA_EL[0] || eta > ETA_EL[NETA_EL]) return false;
		return true;
	}
	else 
	{
		if (pt < PT_MU[0] || pt > PT_MU[NPT_MU]) return false;
		if (eta < ETA_MU[0] || eta > ETA_MU[NETA_MU]) return false;
		return true;
	}
}
