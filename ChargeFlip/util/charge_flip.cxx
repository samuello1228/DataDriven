#include "TChain.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TH1D.h"
#include "TH2D.h"
#include "ChargeFlip/susyEvts.h"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdio>

#define PRINT(x) {std::cout << #x << " = " << x << std::endl;}

using namespace std;

typedef vector<vector<double> > matrix;
typedef vector<double> row;

row VETAS = {0, 0.5, 1, 1.37, 1.52, 1.8, 2.0, 2.47};
row VPTS  = {25, 60, 90, 130, 150, 1000};

const string CHAIN_NAME = "evt2l";
const double M_Z       = 91.1876;
const double M_Z_WIDTH = 10;


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
string file_format(string data_type, string elec_type, string name, double lZcand_M, double rZcand_M, double bl, double br, string extension);
void   bin_zmass(double mll, int bid1, int bid2, double lZcand_M, double rZcand_M, 
		matrix &central,
		matrix &left,
		matrix &right);
void   print_setting(double lZcand_M, double rZcand_M, double bl = 0.0, double br = 0.0, 
		bool isBkgSub = true);
void   write_to_file(string outFile, string data_type, string elec_type, string n_name, string nss_name, string des,	
		matrix n,
		matrix nss,
		double lZcand_M,
		double rZcand_M,
		double bl = 0.0,
		double br = 0.0); 
void   write_matrix(string ofn, const matrix& M);
void   print_write_info(string header, string info, const matrix& central, string name);
matrix s_mul(double s, const matrix& A);
matrix m_add(const matrix& A, const matrix& B);
matrix m_sub(const matrix& A, const matrix& B);
matrix negative_to_zero(const matrix& M);
matrix sub_bkg(double bl, double br, double width, 
		const matrix& central, 
		const matrix& left,
		const matrix& right);
TChain* LoadData(string file);
void    processEvents(string outFile,
		string   data_type, 
		string   elec_type,
		susyEvts *mEvts,
		row&     etas,
		row&     pts,
		double   lZcand_M,
		double   rZcand_M,
		double   bl       = 0.0,
		double   br       = 0.0,
		bool     isBkgSub = true);
int    get_charge(susyEvts *mEvts, int i);
double howNear(const TLorentzVector &w1, const TLorentzVector &w2);
int    match(susyEvts *mEvts, int i);
int    preJudge(susyEvts *mEvts, int i, int j, int &t_char_1, int &t_char_2);
//-----------------------------------


int main(int argc, char* argv[])
{
	if (argc != 5)
	{
		cerr << "Wrong number of arguments. " << argc - 1 << " provided, but 4 are required." << endl;
		cerr << "Command format : " << argv[0] << " <data_type> <elec_type> <input_file> <output_file> " << endl;
		cerr << "data_type: data or mc " << endl;
		cerr << "elec_type: signal or baseline " << endl;
		return -1;
	}
	string data_type = string(argv[1]);
	string elec_type = string(argv[2]);
	string file      = string(argv[3]);
	string output    = string(argv[4]);

	TChain *tc = LoadData(file);

	susyEvts* mEvts = new susyEvts(tc);
	long long nEntries = (mEvts->tree1)->GetEntries();
	cout << "Total entries: " << nEntries << endl;

	if (data_type == "data")
	{
		processEvents(output, data_type, elec_type, mEvts, VETAS, VPTS, 80, 100, 20, 20, true);
		processEvents(output, data_type, elec_type, mEvts, VETAS, VPTS, 80, 100, 25, 25, true);
		processEvents(output, data_type, elec_type, mEvts, VETAS, VPTS, 80, 100, 15, 15, true);
		processEvents(output, data_type, elec_type, mEvts, VETAS, VPTS, 75, 105, 20, 20, true);
		processEvents(output, data_type, elec_type, mEvts, VETAS, VPTS, 80, 100, 0, 0, false);
	}
	else
	{
		processEvents(output, data_type, elec_type, mEvts, VETAS, VPTS, 80, 100, 0, 0, false);
	}

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

void processEvents(string outFile,
		string   data_type,
		string   elec_type,
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
	bool isData;
	if (data_type == "data")
	{
		isData = true;
	}
	else
	{
		isData = false;
	}
	sprintf(checkName, "../run/checks/%s_%s_%d_%d_%d_%d.root", data_type.c_str(), elec_type.c_str(), (int)lZcand_M, (int)rZcand_M, (int)bl, (int)br);
	TFile *f = new TFile(checkName, "RECREATE");
	vector<TH1D*> hEta;
	vector<TH1D*> hPt;
	vector<TH1D*> hMass;
	hEta.clear();
	hPt.clear();
	hMass.clear();
	TH1I* hCutflow = new TH1I("hCutflow", "Number of events passed", 12, 0, 12);
	hCutflow->GetXaxis()->SetBinLabel(1, "Total");
	hCutflow->GetXaxis()->SetBinLabel(2, "Trig+GRL+2e");
	hCutflow->GetXaxis()->SetBinLabel(3, "Pass LooseBaseline");
	hCutflow->GetXaxis()->SetBinLabel(4, "Pass EID");
	hCutflow->GetXaxis()->SetBinLabel(5, "Pass Pt+Eta");
	hCutflow->GetXaxis()->SetBinLabel(6, "Signal+SB");
	hCutflow->GetXaxis()->SetBinLabel(7, "Signal");
	hCutflow->GetXaxis()->SetBinLabel(8, "Left SB");
	hCutflow->GetXaxis()->SetBinLabel(9, "Right SB");
	hCutflow->GetXaxis()->SetBinLabel(10,"SS ZSB");

	hEta.push_back(new TH1D("hEtaAll", "Eta distribution of all electrons", 200, -2.47, 2.47));
	hEta.push_back(new TH1D("hEtaPreselected", "Eta distribution of preselected electrons", 200, -2.47, 2.47));
	hEta.push_back(new TH1D("hEtaLoose", "Eta distribution of electrons in LooseBaseline pairs", 200, -2.47, 2.47));
	hEta.push_back(new TH1D("hEtaSignal", "Eta distribution of electrons passing EID", 200, -2.47, 2.47));
	hEta.push_back(new TH1D("hEtaSignalPtEta", "Eta distribution of electrons in Signal pairs passing Pt+Eta", 200, -2.47, 2.47));
	hEta.push_back(new TH1D("hEtaSignalZSB", "Eta distribution of electrons in Signal pairs within Zmass+SB window", 200, -2.47, 2.47));
	hEta.push_back(new TH1D("hEtaSignalZ", "Eta distribution of electrons in Signal pairs within Zmass window", 200, -2.47, 2.47));
	hEta.push_back(new TH1D("hEtaSignalLSB", "Eta distribution of electrons in Signal pairs within left SB", 200, -2.47, 2.47));
	hEta.push_back(new TH1D("hEtaSignalRSB", "Eta distribution of electrons in Signal pairs within right SB", 200, -2.47, 2.47));
	hEta.push_back(new TH1D("hEtaSSZSB", "Eta distribution of electrons in SS pairs within Zmass+SB", 200, -2.47, 2.47));

	hPt.push_back(new TH1D("hPtAll", "Pt distribution of all electrons", 200, 0, 500));
	hPt.push_back(new TH1D("hPtPreselected", "Pt distribution of preselected electrons", 200, 0, 500));
	hPt.push_back(new TH1D("hPtLoose", "Pt distribution of electrons in LooseBaseline pairs", 200, 0, 500));
	hPt.push_back(new TH1D("hPtSignal", "Pt distribution of electrons passing EID", 200, 0, 500));
	hPt.push_back(new TH1D("hPtSignalPtEta", "Pt distribution of electrons in Signal pairs passing Pt+Eta", 200, 0, 500));
	hPt.push_back(new TH1D("hPtSignalZSB", "Pt distribution of electrons in Signal pairs within Zmass+SB window", 200, 20, 200));
	hPt.push_back(new TH1D("hPtSignalZ", "Pt distribution of electrons in Signal pairs within Zmass window", 200, 20, 200));
	hPt.push_back(new TH1D("hPtSignalLSB", "Pt distribution of electrons in Signal pairs within left SB", 200, 20, 200));
	hPt.push_back(new TH1D("hPtSignalRSB", "Pt distribution of electrons in Signal pairs within right SB", 200, 20, 200));
	hPt.push_back(new TH1D("hPtSSZSB", "Pt distribution of electrons in SS pairs within Zmass+SB window", 200, 20, 200));

	hMass.push_back(new TH1D("hMassAll", "Mass distribution of all electrons", 200, 20, 200));
	hMass.push_back(new TH1D("hMassPreselected", "Mass distribution of preselected electrons", 200, 20, 200));
	hMass.push_back(new TH1D("hMassLoose", "Mass distribution of electrons in LooseBaseline pairs", 200, 20, 200));
	hMass.push_back(new TH1D("hMassSignal", "Mass distribution of electrons in EID", 200, 20, 200));
	hMass.push_back(new TH1D("hMassSignalQID", "Mass distribution of electrons in Signal pairs passing Pt+Eta", 200, 20, 200));
	hMass.push_back(new TH1D("hMassSignalZSB", "Mass distribution of electrons in Signal pairs within Zmass+SB window", 200, 20, 200));
	hMass.push_back(new TH1D("hMassSignalZ", "Mass distribution of electrons in Signal pairs within Zmass window", 200, 20, 200));
	hMass.push_back(new TH1D("hMassSignalLSB", "Mass distribution of electrons in Signal pairs within left SB", 200, 20, 200));
	hMass.push_back(new TH1D("hMassSignalRSB", "Mass distribution of electrons in Signal pairs within right SB", 200, 20, 200));
	hMass.push_back(new TH1D("hMassSSZSB", "Mass distribution of electrons in SS pairs within Zmass+SB window", 200, 20, 200));

	long long nEntries = (mEvts->tree1)->GetEntries();
	long long num[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	double MINETA = row_min(etas);
	double MAXETA = row_max(etas);
	double MINPT  = row_min(pts);
	double MAXPT  = row_max(pts);
	int NETA = etas.size();
	int NPT  = pts.size();
	int SIZE = (NETA - 1)*(NPT - 1);

	matrix nC(SIZE, row(SIZE, 0.0));
	matrix nL(SIZE, row(SIZE, 0.0));
	matrix nR(SIZE, row(SIZE, 0.0));

	matrix nssC(SIZE, row(SIZE, 0.0));
	matrix nssL(SIZE, row(SIZE, 0.0));
	matrix nssR(SIZE, row(SIZE, 0.0));

	TH2D *hTruth = 0;
	TH2D *hN = 0;
	TH2D *hNflipped = 0;
	stringstream monitor;
	monitor.clear();
	int n_monitor = 0;
	if (!isData)
	{
		hTruth    = new TH2D("hMCTruthRate", "MC Truth misID rate", etas.size()-1, &etas[0], pts.size()-1, &pts[0]);
		hN        = (TH2D*) hTruth->Clone("hMC_N");
		hNflipped = (TH2D*) hTruth->Clone("hMC_Nflipped");
	}

#define CUTFLOW(i)                                                                          \
	hCutflow->Fill(i);                                                                      \
	hEta[i]->Fill(mEvts->leps[0].eta, weight); hEta[i]->Fill(mEvts->leps[1].eta, weight);   \
	hPt[i]->Fill(mEvts->leps[0].pt, weight); hPt[i]->Fill(mEvts->leps[1].pt, weight);       \
	hMass[i]->Fill(mll, weight);                                                            \
	num[i]++;

	for (long long i = 0; i < nEntries; i++)
	{
		mEvts->GetEntry(i);
		double weight = 1.0; 
		double mll = getMll(mEvts);
		double e1_eta = fabs(mEvts->leps[0].eta);
		double e1_pt  = mEvts->leps[0].pt;
		double e2_eta = fabs(mEvts->leps[1].eta);
		double e2_pt  = mEvts->leps[1].pt;
		int bid1, bid2;

		CUTFLOW(0); //total

		if (mEvts->sig.trigCode <= 0) continue;
		if (mEvts->leps.size() != 2) continue;
		if (int(fabs(mEvts->leps[0].ID/1000))!=11 || int(fabs(mEvts->leps[1].ID/1000))!=11) continue;
		CUTFLOW(1); //trig+2e

		CUTFLOW(2);

		if (elec_type == "signal")
		{
			if (!((mEvts->leps[0].lFlag & 2)/2)) continue;
			if (!((mEvts->leps[1].lFlag & 2)/2)) continue;
			CUTFLOW(3); //signal e
		}

		if (is_out_eta_pt(e1_eta, e1_pt, e2_eta, e2_pt, MINETA, MAXETA, MINPT, MAXPT)) continue;
		CUTFLOW(4); //eta, pt requirement

		if (is_out_zmass(mll, lZcand_M, rZcand_M, bl, br)) continue;
		CUTFLOW(5); //Zmass+SB

		if (mll > lZcand_M && mll < rZcand_M) 
		{
			CUTFLOW(6);
		}
		if (mll > lZcand_M-bl && mll < lZcand_M) 
		{
			CUTFLOW(7);
		}
		if (mll > rZcand_M && mll < rZcand_M+br) 
		{
			CUTFLOW(8);
		}

		bid1 = bin_id(eta_bin(e1_eta), pt_bin(e1_pt), NPT);
		bid2 = bin_id(eta_bin(e2_eta), pt_bin(e2_pt), NPT);

		bin_zmass(mll, bid1, bid2, lZcand_M, rZcand_M, nC, nL, nR);

		int e1_charge = (mEvts->leps[0].ID > 0) - (mEvts->leps[0].ID < 0);
		int e2_charge = (mEvts->leps[1].ID > 0) - (mEvts->leps[1].ID < 0);

		if (e1_charge == e2_charge)
		{
			CUTFLOW(9);
			bin_zmass(mll, bid1, bid2, lZcand_M, rZcand_M, nssC, nssL, nssR);
		}

		if (!isData)
		{
			int e1_charge_truth, e2_charge_truth;
			int status = preJudge(mEvts, 0, 1, e1_charge_truth, e2_charge_truth);
			if (status == -1)
			{
				e1_charge_truth = get_charge(mEvts, 0);
				e2_charge_truth = get_charge(mEvts, 1);
				num[10]++;
			}
			// judge on MC truth

			hN->Fill(e1_eta, e1_pt, weight);
			hN->Fill(e2_eta, e2_pt, weight);
			if (e1_charge == e2_charge) 
			{
				// solve no mother Z info 
				if (e1_charge_truth == 77 && abs(e2_charge_truth) == 1)
				{
					e1_charge_truth = -e2_charge_truth;
				}
				if (e2_charge_truth == 77 && abs(e1_charge_truth) == 1)
				{
					e2_charge_truth = -e1_charge_truth;
				}
				// cannot judge both missing match info 
				if (abs(e1_charge_truth) != 1 && abs(e2_charge_truth) != 1)
				{
					//num[11]++;
					monitor << i << "\t";
					monitor << e1_charge << "\t" << e1_charge_truth << "\t";
					monitor << e2_charge << "\t" << e2_charge_truth << "\t" << "X" << endl;
					n_monitor++;
					continue;
				}

				if (e1_charge != e1_charge_truth)
				{
					hNflipped->Fill(e1_eta, e1_pt, weight);
					monitor << i << "\t";
					monitor << e1_charge << "\t" << e1_charge_truth << "\t";
					monitor << e2_charge << "\t" << e2_charge_truth << endl;
				}
				else if (e2_charge != e2_charge_truth)
				{
					hNflipped->Fill(e2_eta, e2_pt, weight);
					monitor << i << "\t";
					monitor << e1_charge << "\t" << e1_charge_truth << "\t";
					monitor << e2_charge << "\t" << e2_charge_truth << endl;
				}
				else
				{
					monitor << i << "\t";
					monitor << e1_charge << "\t" << e1_charge_truth << "\t";
					monitor << e2_charge << "\t" << e2_charge_truth << "\t" << "X" << endl;
				}
			}
		}
	}
#undef CUTFLOW

	cout << endl << endl << "=======================================> RESULTS: " << endl << endl;
	print_setting(lZcand_M, rZcand_M, bl, br, isBkgSub);

	if (!isBkgSub)
	{
		write_to_file(outFile, data_type, elec_type, "n_nbg", "nss_nbg", "NO BKG SUBTRACTION", nC, nssC, lZcand_M, rZcand_M, bl, br);
	}
	else
	{
		matrix n_s   = sub_bkg(bl, br, rZcand_M - lZcand_M, nC, nL, nR);
		matrix nss_s = sub_bkg(bl, br, rZcand_M - lZcand_M, nssC, nssL, nssR);
		write_to_file(outFile, data_type, elec_type, "n_bg", "nss_bg", "BKG SUBTRACTION", n_s, nss_s, lZcand_M, rZcand_M, bl, br);

	}

	if (!isData)
	{
		hN->Sumw2();
		hNflipped->Sumw2();
		hTruth->Sumw2();
		hTruth->Add(hNflipped);
		hTruth->Divide(hN);
		hTruth->Write();
		hNflipped->Write();
		hN->Write();

		cout << "monitor number: " << n_monitor << endl;
		ofstream fm("../run/csv/mc_monitor.txt");
		fm << monitor.str();
		fm.close();
	}

	for (unsigned int i = 0; i < 12; i++)
	{
		PRINT(num[i]);
	}

	for(auto h : hPt) h->Write();
	for(auto h : hEta) h->Write();
	for(auto h : hMass) h->Write();
	hCutflow->Write();

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

void bin_zmass(double mll, int bid1, int bid2, double lZcand_M, double rZcand_M, 
		matrix &central,
		matrix &left,
		matrix &right)
{
	double weight = 1;
	if ((mll > lZcand_M) && (mll < rZcand_M))
		central[bid1][bid2] += weight;
	else if (mll < lZcand_M)
		left[bid1][bid2] += weight;
	else if (mll > rZcand_M)
		right[bid1][bid2] += weight;
	else return;
}

void print_setting(double lZcand_M, double rZcand_M, double bl, double br, bool isBkgSub)
{
	cout << "Setting -- LMLL: " << lZcand_M << " -- RMLL: " << rZcand_M;
	if (isBkgSub)
		cout << " -- BL: " << bl << " -- BR: " << br << endl;
	else
		cout << endl;
}

void write_to_file(string outFile,
		string data_type,
		string elec_type,
		string n_name,
		string nss_name,
		string des,
		matrix n,
		matrix nss,
		double lZcand_M,
		double rZcand_M,
		double bl,
		double br) 
{
	ofstream binFile;
	binFile.open(outFile, ios_base::app);

	string n_n = file_format(data_type, elec_type, n_name, lZcand_M, rZcand_M, bl, br, "csv");
	n_n = "../run/csv/" + n_n;
	binFile << n_n << endl;
	print_write_info(des, "N", n, n_n);

	string nss_n = file_format(data_type, elec_type, nss_name, lZcand_M, rZcand_M, bl, br, "csv");
	nss_n = "../run/csv/" + nss_n;
	binFile << nss_n << endl;
	print_write_info(des, "NSS", nss, nss_n);
}

string file_format(string data_type, string elec_type, string name, double lZcand_M, double rZcand_M, double bl, double br, string extension)
{
	string slZcand_M = to_string((int) lZcand_M);
	string srZcand_M = to_string((int) rZcand_M);
	string sbl       = to_string((int) bl);
	string sbr       = to_string((int) br);
	string s         = "_";
	string output    = data_type + s + elec_type + s + name + s + slZcand_M + s + srZcand_M + s + sbl + s + sbr + "." + extension;

	return output;
}

void write_matrix(string ofn, const matrix& M)
{
	ofstream text;
	text.open(ofn);
	text << std::fixed << std::setprecision(1);
	for (size_t i = 0; i < M.size(); i++)
	{
		text << M[i][0];
		for (size_t j = 1; j < M[i].size(); j++)
		{
			text << " " << M[i][j];
		}
		text << endl;
	}
}

void print_write_info(string header, string info, const matrix& central, string name)
{
	cout << endl << endl;
	cout << "------------------------------------" << endl;
	cout << "------------------------------------" << endl;
	cout << header << endl;
	cout << "----------" << endl;
	write_matrix(name, central);
	cout << info << " has been written to file " << name << endl;
	cout << "----------" << endl;
}

matrix s_mul(double s, const matrix& A)
{
	matrix M(A.size(), row(A[0].size(), 0.0));
	for(size_t i = 0; i < M.size(); i++)
	{
		for (size_t j = 0; j < M[i].size(); j++)
			M[i][j] = A[i][j] * s;
	}
	return M;

}

matrix m_add(const matrix& A, const matrix& B)
{
	matrix M(A.size(), row(A[0].size(), 0.0));
	for (size_t i = 0; i < M.size(); i++)
	{
		for (size_t j = 0; j < M[i].size(); j++)
			M[i][j] = A[i][j] + B[i][j];
	}
	return M;
}

matrix m_sub(const matrix& A, const matrix& B)
{
	return m_add(A, s_mul(-1.0, B));
}

matrix negative_to_zero(const matrix& M)
{
	matrix N(M.size(), row(M[0].size(), 0.0));
	for (size_t i = 0; i < M.size(); i++)
	{
		for (size_t j = 0; j < M[i].size(); j++)
		{
			if (M[i][j] > 0)
				N[i][j] = M[i][j];
		}
	}
	return N;
}

matrix sub_bkg(double bl, double br, double width,
		const matrix& central,
		const matrix& left,
		const matrix& right)
{
	matrix M = m_sub(central, s_mul(width/(bl+br), m_add(left, right)));
	return negative_to_zero(M);
}

int get_charge(susyEvts *mEvts, int i)
{
	// Find matching truth particle
	int tp = mEvts->leps[i].truthI;
	if (tp < 0) return 99; 
	int id = mEvts->truths[tp].pdgId;
	if (id != 11 && id != -11)
	{
		return 88;
	}
	else 
	{
		int mother = mEvts->truths[tp].motherI;
		int m_id   = mEvts->truths[mother].pdgId;
		bool f_z     = false;
		bool f_gamma = false;
		int index_gamma;
		if (mother == -1)
		{
			return 77;
		}
		while (mother != -1)
		{
			if (m_id == 23) //tag z
			{
				f_z = true;
			}
			if (m_id == 22)
			{
				f_gamma     = true;
				index_gamma = mother;
			}
			mother = mEvts->truths[mother].motherI;
			m_id   = mEvts->truths[mother].pdgId;
		}
		if (!f_z) //no z tag 
		{
			return 66;
		}
		if (f_gamma) //gamma conversion 
		{
			mother   = mEvts->truths[index_gamma].motherI;
			m_id     = mEvts->truths[mother].pdgId;
			if (abs(m_id) == 11) //has e mother 
			{
				id = m_id;
			}
			else 
			{
				return 55;
			}
		}
		if (id == -11)
		{
			return 1;
		}
		else 
		{
			return -1;
		}

	}
}

int match(susyEvts *mEvts, int i)
{
	TLorentzVector p1, p2;
	p1.SetPtEtaPhiM(mEvts->leps[i].pt, mEvts->leps[i].eta, mEvts->leps[i].phi, 0.000511);
	double min = 999.;
	int flag = -1;

	for (unsigned int i = 0; i < mEvts->truths.size(); i++)
	{
		p2.SetPtEtaPhiM(mEvts->truths[i].pt, mEvts->truths[i].eta, mEvts->truths[i].phi, 0.000511);
		double tmp = howNear(p1, p2);
		if (tmp < min)
		{
			min  = tmp;
			flag = i;
		}
	}
	return flag;
}

double howNear(const TLorentzVector &w1, const TLorentzVector &w2)
{
	double wdw = std::fabs(w1.Vect().Dot(w2.Vect()) + 0.25*((w1.T()+w2.T())*(w1.T()*w2.T())));
	double delta = (w1.Vect() - w2.Vect()).Mag2() + (w1.T()-w2.T())*(w1.T()-w2.T());
	if ( (wdw > 0) && (delta < wdw)  ) 
	{
		return std::sqrt(delta/wdw);
	} 
	else if ( (wdw == 0) && (delta == 0) ) 
	{
		return 0;
	} 
	else 
	{
		return 1;
	}
}

int preJudge(susyEvts *mEvts, int i, int j, int &t_char_1, int &t_char_2)
{
	//0, 1-----match to which electron 
	//-1  -----no match 
	int tp_1 = mEvts->leps[i].truthI;
	int tp_2 = mEvts->leps[j].truthI;
	int mother_1, m_id_1;
	int mother_2, m_id_2;
	int id_1, id_2;
	if (tp_1 < 0 && tp_2 < 0) 
	{	
		return -1;
	}
	else if(tp_1 >= 0 && tp_2 < 0)
	{
		id_1 = mEvts->truths[tp_1].pdgId;
		if (abs(id_1) == 11) //electron
		{
			mother_1 = mEvts->truths[tp_1].motherI;
			m_id_1   = mEvts->truths[mother_1].pdgId;
			while (abs(m_id_1) == 11)
			{
				id_1     = m_id_1;
				mother_1 = mEvts->truths[mother_1].motherI;
				m_id_1   = mEvts->truths[mother_1].pdgId;	
			}
			if (m_id_1 == 23) //z0
			{
				if (id_1 == 11)
				{
					t_char_1 = -1;
					t_char_2 = 1;
				}
				else 
				{
					t_char_1 = 1;
					t_char_2 = -1;
				}
				return 0;
			}
			else 
			{
				return -1;
			}
		}
		else 
		{
			return -1;
		}
	}
	else if(tp_1 < 0 && tp_2 >= 0)
	{
		id_2 = mEvts->truths[tp_2].pdgId;
		if (abs(id_2) == 11) //electron
		{
			mother_2 = mEvts->truths[tp_2].motherI;
			m_id_2   = mEvts->truths[mother_2].pdgId;
			while (abs(m_id_2) == 11)
			{
				id_2     = m_id_2;
				mother_2 = mEvts->truths[mother_2].motherI;
				m_id_2   = mEvts->truths[mother_2].pdgId;	
			}
			if (m_id_2 == 23) //z0
			{
				if (id_2 == 11)
				{
					t_char_2 = -1;
					t_char_1 = 1;
				}
				else 
				{
					t_char_2 = 1;
					t_char_1 = -1;
				}
				return 1;
			}
			else 
			{
				return -1;
			}
		}
		else 
		{
			return -1;
		}
	}
	else 
	{
		id_1 = mEvts->truths[tp_1].pdgId;
		id_2 = mEvts->truths[tp_2].pdgId;
		mother_1 = mEvts->truths[tp_1].motherI;
		m_id_1   = mEvts->truths[mother_1].pdgId;
		mother_2 = mEvts->truths[tp_2].motherI;
		m_id_2   = mEvts->truths[mother_2].pdgId;
		while (abs(m_id_1) == 11)
		{
			id_1     = m_id_1;
			mother_1 = mEvts->truths[mother_1].motherI;
			m_id_1   = mEvts->truths[mother_1].pdgId;	
		}
		while (abs(m_id_2) == 11)
		{
			id_2     = m_id_2;
			mother_2 = mEvts->truths[mother_2].motherI;
			m_id_2   = mEvts->truths[mother_2].pdgId;	
		}
		if (abs(id_1) == 11 && m_id_1 == 23)
		{
			if (id_1 == 11)
			{
				t_char_1 = -1;
				t_char_2 = 1;
			}
			else 
			{
				t_char_1 = 1;
				t_char_2 = -1;
			}
			return 0;
		}
		else if (abs(id_2) == 11 && m_id_2 == 23)
		{
			if (id_2 == 11)
			{
				t_char_2 = -1;
				t_char_1 = 1;
			}
			else 
			{
				t_char_2 = 1;
				t_char_1 = -1;
			}
			return 1;
		}
		else 
		{
			return -1;
		}
	}	
}
