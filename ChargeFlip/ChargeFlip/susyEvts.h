#ifndef susyEvts_H
#define susyEvts_H
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TEntryList.h"
#include "obj_def.h"
#include <string>

class susyEvts
{
	public:
		susyEvts() { tree2 = 0; }
		susyEvts(TTree* tr); 
		~susyEvts();

		TString version;
		TTree* makeTree(TString treename);
		TTree* makeWeightOnlyTree(TString treename, susyEvts* parent);
		TTree* makePtCorrTree(TString treename, susyEvts* parent);
		TTree* makeKinematicsSysTree(TString treename, susyEvts* parent);
		Int_t GetEntry(Long64_t entry = 0, Int_t getall = 0){ return tree1->GetEntry(entry, getall); }
		int fill(){ return tree2->Fill(); }
		void writeTree(TString treename){ tree2->Write(treename);}
		Int_t Next();
		void resetIter(){ m_elist = tree1->GetEntryList(); m_el = 0; }

		TTree* tree1;
		TTree* tree2;

		EVT evt;
		SIGNATURE sig;
		LEPTONS leps;
		LEPTONS* vleps = &leps;
		R_PAR l12;
		JETS jets;
		JETS* vjets = &jets;
		TRUTHS truths;
		TRUTHS* vtruths = &truths;

	private:
		void getTree(TTree* t);
		TEntryList* m_elist;
		Long64_t m_el;
};
#endif //susyEvts_H
