// Author    : mat@ihep.ac.cn
// Date      : May 15, 2018 
// Descrption: It is used to obtain charge-flip rates by fitting based on selected Z->ee control samples
// Inputs    : .csv files generated by program ../../util/charge_flip.cxx 
//           : locate at the directory of ../../run/csv
// Outputs   : .txt files which records charge-flips rate in each (eta, pt) bin 

#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TMinuit.h"

using namespace std;

typedef vector<vector<double> > matrix;
typedef vector<double> row;

row VETAS = {0, 0.5, 1, 1.37, 1.52, 1.8, 2.0, 2.47};
row VPTS  = {25, 60, 90, 130, 150, 1000};
unsigned int NETA  = VETAS.size();
unsigned int NPTS  = VPTS.size();
unsigned int NVAR  = (NETA - 1) * (NPTS - 1);

matrix matrix_n(NVAR, row(NVAR, 0.0));
matrix matrix_nss(NVAR, row(NVAR, 0.0));

int loadmatrix(string path_n, string path_nss);
void showmatrix(const matrix& m);
double likelihood(const double *x);
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
int doFit(bool isData, string elec_type);

int fit()
{
	// fit sequence
	// data based on signal electron
	doFit(true, "signal");
	// data based on baseline electron
	doFit(true, "baseline");
	// mc based on signal electron
	doFit(false, "signal");
	// mc based on baseline electron
	doFit(false, "baseline");
	
	return 0;
}

int doFit(bool isData, string elec_type)
{
	vector<string> path;
	path.clear();
	string path_n, path_nss, output;
	if (isData)
	{
		path.push_back( "80_100_20_20" );
		path.push_back( "80_100_15_15" );
		path.push_back( "80_100_25_25" );
		path.push_back( "75_105_20_20" );
		path.push_back( "80_100_0_0" );
	}
	else
	{
		path.push_back( "80_100_0_0" );
	}

	for (unsigned int k = 0; k < path.size(); k++)
	{
		if (isData)
		{
			if (k != 4)
			{
				path_n   = string("../../run/csv/data_") + elec_type + string("_n_bg_") + path[k] + string(".csv");
				path_nss = string("../../run/csv/data_") + elec_type + string("_nss_bg_") + path[k] + string(".csv");
			}
			else 
			{
				path_n   = string("../../run/csv/data_") + elec_type + string("_n_nbg_") + path[k] + string(".csv");
				path_nss = string("../../run/csv/data_") + elec_type + string("_nss_nbg_") + path[k] + string(".csv");
			}
			output   = string("fit_data_") + elec_type + string("_") + path[k] + string(".txt");
		}
		else 
		{
			path_n   = string("../../run/csv/mc_") + elec_type + string("_n_nbg_") + path[k] + string(".csv");
			path_nss = string("../../run/csv/mc_") + elec_type + string("_nss_nbg_") + path[k] + string(".csv");
			output   = string("fit_mc_") + elec_type + string("_") + path[k] + string(".txt");
		}

		if (loadmatrix(path_n, path_nss) == -1)
			return -1;

		static Double_t p0=0;
		static Double_t p1=1;
		static Double_t p3=3;
		//static Double_t p500=500;

		TMinuit *gMinuit = new TMinuit(NVAR);
		gMinuit->SetFCN(fcn);

		Int_t iflag = 0;
		Double_t *startval = new double[NVAR];
		std::fill_n(startval, NVAR, 5e-5);

		gMinuit->mnexcm("SET ERR", &p1, 1, iflag);
		gMinuit->SetErrorDef(0.5);

		for (unsigned int i = 0; i < NVAR; i++)
		{
			if (elec_type == "signal")
			{
				gMinuit->mnparm(i, Form("x%d", i), startval[i], 1e-5, 1e-5, 0.1, iflag);
			}
			else 
			{
				gMinuit->mnparm(i, Form("x%d", i), startval[i], 1e-5, 1e-5, 0.2, iflag);
			}
		}
		gMinuit->mnexcm("CALL FCN", &p1 , 1, iflag);
		gMinuit->mnexcm("SET PRINT",&p0 , 1, iflag);
		gMinuit->mnexcm("MIGRAD",   &p0 , 0, iflag);
		gMinuit->mnexcm("MINOS",    &p0 , 0, iflag);
		gMinuit->mnexcm("MIGRAD",   &p0 , 0, iflag);
		gMinuit->mnexcm("MINOS",    &p0 , 0, iflag);
		gMinuit->mnexcm("CALL FCN", &p3 , 1, iflag);


		cout << "Fit Status: " << gMinuit->GetStatus() << endl;
		// output of fit parameters
		double *par     = new double[NVAR];
		double *par_err = new double[NVAR];
		ofstream out(output.c_str());
		for (unsigned int i = 0; i < NVAR; i++)
		{
			gMinuit->GetParameter(i, par[i], par_err[i]);
			out << par[i] << "\t";
			out << par_err[i] << endl;
		}
		out.close();
	}
	return 0;
}

double likelihood(const double *x)
{
	double result = 0.;
	for (unsigned int i = 0; i < NVAR; i++)
	{
		for (unsigned int j = 0; j< NVAR; j++)
		{
			double n       = matrix_n[i][j];
			double nss     = matrix_nss[i][j];;
			double nss_exp = n * ((1-x[i])*x[j] + (1-x[j])*x[i]);
			if (fabs(n) < 1e-5) continue;
			double item;
			item = nss * log(nss_exp) - nss_exp;
			result -= item;
		}
		//negative value
	}
	return result; //for accurate error estimation
}

int loadmatrix(string path_n, string path_nss)
{
	ifstream file_n(path_n.c_str());
	if (file_n.is_open())
	{
		string line;
		unsigned int i = 0;
		while(getline(file_n, line))
		{
			stringstream tmp(line);
			for (unsigned j = 0; j < NVAR; j++)
			{
				tmp >> matrix_n[i][j];
			}
			i++;
		}
	}
	else
	{
		cerr << "Cannot open file: " << path_n << endl;
		return -1;
	}
	file_n.close();
	ifstream file_nss(path_nss.c_str());
	if (file_nss.is_open())
	{
		string line;
		unsigned int i = 0;
		while(getline(file_nss, line))
		{
			stringstream tmp(line);
			for (unsigned j = 0; j < NVAR; j++)
			{
				tmp >> matrix_nss[i][j];
			}
			i++;
		}
	}
	else
	{
		cerr << "Cannot open file: " << path_nss << endl;
	}
	file_nss.close();

	return 0;
}

void showmatrix(const matrix& m)
{
	for (unsigned int i = 0; i < m.size(); i++)
	{
		cout << "Line : " << i << "\t";
		for (unsigned int j = 0; j < m[0].size(); j++)
		{
			cout << m[i][j] << "\t";
		}
		cout << endl;
	}
}

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
	f = likelihood(par);
}