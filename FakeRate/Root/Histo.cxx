#include "FakeRate/Histo.h"

Histo::Histo(TString name)
{
	isEffValid = false;

	hTight = new TH2D(name+"_hTight", ";p_{T} [GeV];|#eta|", NPT, PT, NETA, ETA);
		
	hLoose = new TH2D(name+"_hLoose", ";p_{T} [GeV];|#eta|", NPT, PT, NETA, ETA);
	hEff   = new TH2D(name+"_hEff",   ";p_{T} [GeV];|#eta|", NPT, PT, NETA, ETA);
	hCom   = new TH2D(name+"_hCom",   ";p_{T} [GeV];|#eta|", NPT, PT, NETA, ETA);
	hTight->Sumw2();
	hLoose->Sumw2();
	hEff->Sumw2();
	hCom->Sumw2();
}

Histo::~Histo()
{
	delete hTight;
	delete hLoose;
	delete hEff;
}

void Histo::Fill(lepInfo &lep)
{
	hLoose->Fill(lep.pt, fabs(lep.eta), lep.w);
	if (lep.isTight)
	{	
		hTight->Fill(lep.pt, fabs(lep.eta), lep.w);
	}
	if (lep.isPassCom)
	{	
		hCom->Fill(lep.pt, fabs(lep.eta), lep.w);
	}
}

void Histo::CalEff() 
{
	if (!isEffValid) 
	{
		hEff->Add(hTight); 
		hEff->Divide(hLoose); 
		isEffValid = true; 
	}
}

// void Histo::SetDir(TDirectory *dir)
// {
	// hLoose->SetDirectory(dir);
	// hTight->SetDirectory(dir);
	// hEff->SetDirectory(dir);
// }
