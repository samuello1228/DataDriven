// Author    : mat@ihep.ac.cn
// Date      : May 13, 2018 
// Descrption: It is used to re-weigh the selcected Z->ee control samples by using opposite-sign ee events to predict same-sign ee events
// Inputs    : ROOT files pre-selected by multiLepSearch package (see https://gitlab.cern.ch/hku/SUSY2L/tree/master/code/multiLepSearch)
// Outputs   : .root file which stores the corresponding histograms with calculated weights 
//           : locate at the directory of ../run/checks_weight
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include "TChain.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TMath.h"
#include "ChargeFlip/susyEvts.h"

#define PRINT(x) {std::cout << #x << " = " << x << std::endl;}

using namespace std;

typedef vector<vector<double> > matrix;
typedef vector<double> row;

row VETAS = {0, 0.5, 1, 1.37, 1.52, 1.8, 2.0, 2.47};
row VPTS  = {25, 60, 90, 130, 150, 1000};
int NETA = VETAS.size();
int NPT  = VPTS.size();
int SIZE = (NETA - 1)*(NPT - 1);

const string CHAIN_NAME = "evt2l";
const double M_Z       = 91.1876;
const double M_Z_WIDTH = 10;

double *corr = new double[SIZE];
double *corr_down = new double[SIZE];
double *corr_up = new double[SIZE];

//------ methods for simplicity -----
bool   is_out_zmass(double mll, double lZcand_M, double rZcand_M, double bl, double br);
bool   is_out_eta_pt(double e1_eta, double e1_pt, double e2_eta, double e2_pt, 
		double min_eta, double max_eta, double min_pt, double max_pt);
double getMll(susyEvts* mEvts, int i = 0, int j = 1);
double row_max(row i);
double row_min(row i);
int    eta_bin(double eta);
int    pt_bin(double pt);
int    bin_id(int etaB, int ptB, int npt);
TChain* LoadData(string file);
void   processEvents(
	string   data_type, 
	susyEvts *mEvts,
	row&     etas,
	row&     pts,
	double   lZcand_M,
	double   rZcand_M,
	double   bl       = 0.0,
	double   br       = 0.0,
	bool     isBkgSub = true);
void   LoadCorr(string file);
//-----------------------------------

int main(int argc, char* argv[])
{
	if (argc != 3)
	{
		cerr << "Wrong number of arguments. " << argc - 1 << " provided, but 2 are required." << endl;
		cerr << "Command format : " << argv[0] << " <data_type> <input_file>" << endl;
		return -1;
	}
	string data_type = string(argv[1]);
	string file      = string(argv[2]);

	//string fit_data = "../scripts/Fit/fit_data_signal_80_100_20_20.txt";
	//string fit_mc   = "../scripts/Fit/fit_mc_signal_80_100_20_20.txt";
	string fit_data = "hFlipProb_data_SysAndStat";
	string fit_mc   = "hFlipProb_MC_SysAndStat";
	string input_corr;

	bool isBkgSub = true;
	pair<double, double> sideband = {20, 20};

	if (data_type == "data")
	{
		input_corr = fit_data;
	}
	else
	{
		input_corr = fit_mc;
	}

	LoadCorr(input_corr);

	TChain *tc = LoadData(file);

	susyEvts* mEvts = new susyEvts(tc);
	long long nEntries = (mEvts->tree1)->GetEntries();
	cout << "Total entries: " << nEntries << endl;

	processEvents(data_type, mEvts, VETAS, VPTS, 80, 100, sideband.first, sideband.second, isBkgSub);

	return 0;
}


TChain *LoadData(string file)
{
	TChain *tc = new TChain(CHAIN_NAME.c_str());
	ifstream in(file);
	if (in.is_open())
	{
		string line;
		while (getline(in, line))
		{
			if (gSystem->AccessPathName(line.c_str(), kFileExists)) //if file exists, return false 
			{
				cout << ">> File: '" << line << "' DO NOT exist!" << endl;
				continue;
			}
			if (tc->Add(line.c_str()))
			{
				cout << ">> File: '" << line << "' :: TTree '" << tc->GetName() << "' has been loaded" << endl;
			}
			else
			{
				cout << ">> File: '" << line << "' :: TTree '" << tc->GetName() << "' cannot be loaded" << endl;
			}
		}
	}
	in.close();	
	return tc;
}

double getMll(susyEvts *mEvts, int i, int j)
{
	TLorentzVector p1, p2;
	p1.SetPtEtaPhiM(mEvts->leps[i].pt, mEvts->leps[i].eta, mEvts->leps[i].phi, 0.000511);
	p2.SetPtEtaPhiM(mEvts->leps[j].pt, mEvts->leps[j].eta, mEvts->leps[j].phi, 0.000511);

	return (p1+p2).M();
}

void processEvents(
	string   data_type,
	susyEvts *mEvts,
	row&     etas,
	row&     pts,
	double   lZcand_M,
	double   rZcand_M,
	double   bl,
	double   br,
	bool     isBkgSub)
{
	char checkName[200];
	sprintf(checkName, "../run/checks_weight/%s_%d_%d_%d_%d.root", data_type.c_str(), (int)lZcand_M, (int)rZcand_M, (int)bl, (int)br);
	TFile *f = new TFile(checkName, "RECREATE");

	struct Variable 
	{
		Variable(TString name_new, int nbin_new, int xlow_new, int xup_new)
		{
			name = name_new;
			nbin = nbin_new;
			xlow = xlow_new;
			xup  = xup_new;
		}
		
		TH1D *h_nom;
		TH1D *h_down;
		TH1D *h_up;
		TH1D *h_nom_ss;
		TString name;
		int nbin;
		int xlow;
		int xup;
		double var;
	};

	vector<Variable> hist;
	hist.push_back( Variable("mll",50,20,200) );
	hist.push_back( Variable("eta_1",50,-2.47,2.47) );
	hist.push_back( Variable("eta_2",50,-2.47,2.47) );
	hist.push_back( Variable("pt_1",50,0,500) );
	hist.push_back( Variable("pt_2",50,0,500) );
	hist.push_back( Variable("phi_1",50,-TMath::Pi(),TMath::Pi()) );
	hist.push_back( Variable("phi_2",50,-TMath::Pi(),TMath::Pi()) );

	for (unsigned int i = 0; i < hist.size(); i++)
	{
		TString hName = "h_" + hist[i].name;
		hist[i].h_nom = new TH1D(hName.Data(), hName.Data(), hist[i].nbin, hist[i].xlow, hist[i].xup);

		TString hName2 = hName + "_down";
		hist[i].h_down = new TH1D(hName2.Data(), hName2.Data(), hist[i].nbin, hist[i].xlow, hist[i].xup);

		hName2 = hName + "_up";
		hist[i].h_up = new TH1D(hName2.Data(), hName2.Data(), hist[i].nbin, hist[i].xlow, hist[i].xup);

		hName2 = hName + "_ss";
		hist[i].h_nom_ss = new TH1D(hName2.Data(), hName2.Data(), hist[i].nbin, hist[i].xlow, hist[i].xup);
	}

	TH2D *h_m_eta_1    = new TH2D("h_m_eta_1",   "h_m_eta", 50, 20, 200, 50, -2.47, 2.47);
	TH2D *h_m_eta_1_ss = new TH2D("h_m_eta_1_ss","h_m_eta", 50, 20, 200, 50, -2.47, 2.47);
	TH2D *h_m_eta_2    = new TH2D("h_m_eta_2",   "h_m_eta", 50, 20, 200, 50, -2.47, 2.47);
	TH2D *h_m_eta_2_ss = new TH2D("h_m_eta_2_ss","h_m_eta", 50, 20, 200, 50, -2.47, 2.47);
	TH2D *h_m_pt_1     = new TH2D("h_m_pt_1",    "h_m_pt",  50, 20, 200, 50, 0, 500);
	TH2D *h_m_pt_1_ss  = new TH2D("h_m_pt_1_ss", "h_m_pt",  50, 20, 200, 50, 0, 500);
	TH2D *h_m_pt_2     = new TH2D("h_m_pt_2",    "h_m_pt",  50, 20, 200, 50, 0, 500);
	TH2D *h_m_pt_2_ss  = new TH2D("h_m_pt_2_ss", "h_m_pt",  50, 20, 200, 50, 0, 500);

	long long nEntries = (mEvts->tree1)->GetEntries();
	long long num[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	double MINETA = row_min(etas);
	double MAXETA = row_max(etas);
	double MINPT  = row_min(pts);
	double MAXPT  = row_max(pts);

	double n_ss[3] = {0, 0, 0};
	double n_os[3] = {0, 0, 0};


	for (long long i = 0; i < nEntries; i++)
	{
		mEvts->GetEntry(i);
		const double weight = mEvts->evt.weight * mEvts->evt.pwt * mEvts->evt.ElSF * mEvts->evt.MuSF * mEvts->evt.BtagSF * mEvts->evt.JvtSF;
		//const double weight = mEvts->evt.weight * mEvts->evt.pwt * mEvts->evt.ElSF * mEvts->evt.MuSF * mEvts->evt.BtagSF * mEvts->evt.JvtSF * mEvts->evt.trigSF;

		double mll = mEvts->l12.m;
		double e1_eta = mEvts->leps[0].eta;
		double e1_pt  = mEvts->leps[0].pt;
		double e1_phi = mEvts->leps[0].phi;
		double e2_eta = mEvts->leps[1].eta;
		double e2_pt  = mEvts->leps[1].pt;
		double e2_phi = mEvts->leps[1].phi;
		int e1_charge = (mEvts->leps[0].ID > 0) - (mEvts->leps[0].ID < 0);
		int e2_charge = (mEvts->leps[1].ID > 0) - (mEvts->leps[1].ID < 0);

		int bid1, bid2;

		num[0]++;
		if (mEvts->sig.trigCode <= 0) continue;
		if (mEvts->leps.size() != 2) continue;
		if (int(fabs(mEvts->leps[0].ID/1000))!=11 || int(fabs(mEvts->leps[1].ID/1000))!=11) continue;
		if (!((mEvts->leps[0].lFlag & 2)/2)) continue;
		if (!((mEvts->leps[1].lFlag & 2)/2)) continue;
		if (is_out_eta_pt(fabs(e1_eta), e1_pt, fabs(e2_eta), e2_pt, MINETA, MAXETA, MINPT, MAXPT)) continue;
		num[1]++;

		bid1 = bin_id(eta_bin(fabs(e1_eta)), pt_bin(e1_pt), NPT);
		bid2 = bin_id(eta_bin(fabs(e2_eta)), pt_bin(e2_pt), NPT);

		hist[0].var = mll;
		hist[1].var = e1_eta;
		hist[2].var = e2_eta;
		hist[3].var = e1_pt;
		hist[4].var = e2_pt;
		hist[5].var = e1_phi;
		hist[6].var = e2_phi;

		if (e1_charge == e2_charge) //ss event
		{
			if (mll > 70 && mll < 80)
			{
				n_ss[0]++;
			}
			else if (mll > 100 && mll < 110)
			{
				n_ss[1]++;
			}
			else if (mll > 80 && mll < 100)
			{
				n_ss[2]++;
			}

			for (unsigned int i = 0; i < hist.size(); i++)
			{
				hist[i].h_nom_ss->Fill(hist[i].var, weight);
			}

			h_m_pt_1_ss->Fill(mll, e1_pt, weight);
			h_m_pt_2_ss->Fill(mll, e2_pt, weight);
			h_m_eta_1_ss->Fill(mll, e1_eta, weight);
			h_m_eta_2_ss->Fill(mll, e2_eta, weight);
		}
		else
		{
			double w_1  = corr[bid1];
			double w_2  = corr[bid2];
			double qFwt = w_1 * (1-w_2) + w_2 * (1-w_1);
			qFwt = qFwt / (1-qFwt);
			double weight_nom = weight * qFwt;

			if (mll > 70 && mll < 80)
			{
				n_os[0]+= weight_nom;
			}
			else if (mll > 100 && mll < 110)
			{
				n_os[1]+= weight_nom;
			}
			else if (mll > 80 && mll < 100)
			{
				n_os[2]+= weight_nom;
			}

			for (unsigned int i = 0; i < hist.size(); i++)
			{
				hist[i].h_nom->Fill(hist[i].var, weight_nom);
			}

			h_m_pt_1->Fill(mll, e1_pt, weight_nom);
			h_m_pt_2->Fill(mll, e2_pt, weight_nom);
			h_m_eta_1->Fill(mll, e1_eta, weight_nom);
			h_m_eta_2->Fill(mll, e2_eta, weight_nom);
		}
	}

	for (unsigned int i = 0; i < hist.size(); i++)
	{
		hist[i].h_nom->Write();
		hist[i].h_nom_ss->Write();
	}

	h_m_pt_1->Write();
	h_m_pt_1_ss->Write();
	h_m_pt_2->Write();
	h_m_pt_2_ss->Write();
	h_m_eta_1->Write();
	h_m_eta_1_ss->Write();
	h_m_eta_2->Write();
	h_m_eta_2_ss->Write();
	
	PRINT(num[0]);
	PRINT(num[1]);

	cout << "sideband and signal counting: " << endl;
	cout << "SS events: " << endl;
	cout << n_ss[0] << "\t" << n_ss[1] << "\t" << n_ss[2] << endl;
	cout << "OS events: " << endl;
	cout << n_os[0] << "\t" << n_os[1] << "\t" << n_os[2] << endl;

	f->Close();
}

bool is_out_zmass(double mll, double lZcand_M, double rZcand_M, double bl, double br)
{
	return (mll < (lZcand_M - bl)) || (mll > (rZcand_M + br));
}

bool is_out_eta_pt(double e1_eta, double e1_pt, double e2_eta, double e2_pt, double min_eta, double max_eta, double min_pt, double max_pt)
{
	return
		(e1_eta < min_eta) || (e1_eta > max_eta) ||
		(e2_eta < min_eta) || (e2_eta > max_eta) ||
		(e1_pt < min_pt) || (e1_pt > max_pt) ||
		(e2_pt < min_pt) || (e2_pt > max_pt);
}

double row_min(row i)
{
	return *min_element(i.begin(), i.end());
}

double row_max(row i) 
{
	return *max_element(i.begin(), i.end());
}

int eta_bin(double eta)
{
	for (row::size_type i = 1 ; i < VETAS.size(); i++)
	{
		if (eta <= VETAS[i])
			return i-1;
	}
	return -1;
}

int pt_bin(double pt)
{
	for (row::size_type i = 1; i < VPTS.size(); i++)
	{
		if (pt <= VPTS[i])
			return (i-1);
	}
	return -1;
}

int bin_id(int etaB, int ptB, int npt)
{
	return (etaB * (npt-1)) + ptB;
}

void LoadCorr(string file)
{
	/*
	ifstream in(file.c_str());
	if (in.is_open())
	{
		string line;
		unsigned int i = 0;
		while(getline(in, line))
		{
			stringstream tmp(line);
			tmp >> corr[i];
			i++;
		}
	}
	else
	{
		cerr << "Cannot open file: " << file << endl;
		exit(1);
	}
	in.close();
	*/

	TFile *f = new TFile("../scripts/GenNTuple/cfRatesSignal.root");
	TH2D *h = (TH2D*) f->Get(file.c_str());

	for (int i = 1; i < NETA; i++)
	{
		for (int j = 1; j < NPT; j++)
		{
			//cout<<h->GetBinContent(j,i)<<", "<<h->GetBinError(j,i)<<endl;
			int bid = (i-1)*(NPT-1) + j-1;
			corr[bid] = h->GetBinContent(j,i);

			if(h->GetBinError(j,i) < h->GetBinContent(j,i)) corr_down[bid] = h->GetBinContent(j,i) - h->GetBinError(j,i);
			else corr_down[bid] = 0;
			
			corr_up[bid]   = h->GetBinContent(j,i) + h->GetBinError(j,i);
			cout<<corr_down[bid]<<", "<<corr[bid]<<", "<<corr_up[bid]<<endl;
		}
	}
	delete f;
}

