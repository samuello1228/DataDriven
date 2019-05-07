// Author    : mat@ihep.ac.cn
// Date      : May 13, 2018 
// Descrption: Acquisition of MC Truth info selected by ChargeFlip package
// Inputs    : ../../run/checks/mc_signal_80_100_0_0.root
//           : ../../run/checks/mc_baseline_80_100_0_0.root
// Output    : ./mc_signal_truth.txt
//           : ./mc_baseline_truth.txt

#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <vector>
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"

using namespace std;

typedef vector<vector<double> > matrix;
typedef vector<double> row;

row VETAS = {0, 0.5, 1, 1.37, 1.52, 1.8, 2.0, 2.47};
row VPTS  = {25, 60, 90, 130, 150, 1000};

int NPTS  = VPTS.size() - 1;
int NETAS = VETAS.size() - 1;

int get_truth()
{
	vector<string> elec_type;
	elec_type.clear();
	elec_type.push_back( "signal" );
	//elec_type.push_back( "baseline" );
	for (unsigned int j = 0; j < elec_type.size(); j++)
	{
		string file   = string("../../run/checks/mc_") + elec_type[j] + string("_80_100_20_20.root");
		string output = string("./mc_") + elec_type[j] + string("_truth.txt");

		TFile *f = new TFile(file.c_str(), "OPEN");
		if (!f->IsOpen())
		{
			cout << "Cannot open file : " << file << endl;
			continue;
		}
		TH2D *h = (TH2D*)f->Get("hMCTruthRate");

		ofstream fout(output);

		for (unsigned int i = 0; i < NETAS; i++)
		{
			for (unsigned int j = 0; j < NPTS; j++)
			{
				fout << h->GetBinContent(i+1, j+1) << "\t";
				fout << h->GetBinError(i+1, j+1) << endl;
			}
		}

		fout.close();
		cout << "File ' " << output << "' has been generated successfully. " << endl; 

		f->Close();
	}
	return 0;
}
