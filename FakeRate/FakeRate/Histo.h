#ifndef __HISTO__
#define __HISTO__
#include "TH2D.h"
#include "TString.h"
#include "TDirectory.h"
#include "math.h"

const double PT[]  = {25, 40, 80, 150};
const double ETA[] = {0, 1.37, 2.5};
const int NPT  = sizeof(PT) / sizeof(PT[0]) - 1;
const int NETA = sizeof(ETA) / sizeof(ETA[0]) - 1;

struct lepInfo
{
	double   pt;
	double   eta;
	double   w;
	bool     isTight;
	bool     isPassCom;
};

class Histo
{
	public:
		Histo(TString name);
		~Histo();
		void    Fill(lepInfo &lep);
		//void    SetDir(TDirectory *dir);
		void    CalEff();
		TH2D*   GetHistEff()   { return hEff; }
		TH2D*   GetHistTight() { return hTight; }
		TH2D*   GetHistLoose() { return hLoose; }
		TH2D*   GetHistCom()   { return hCom; }

	private:
		TH2D*   hTight;
		TH2D*   hLoose;
		TH2D*   hCom;
		TH2D*   hEff;
		bool    isEffValid;

};

#endif
