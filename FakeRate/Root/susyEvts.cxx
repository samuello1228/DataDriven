#include "FakeRate/susyEvts.h"

susyEvts::susyEvts(TTree *tr)
{
	getTree(tr);
	tree2 = 0;
	resetIter();
}

susyEvts::~susyEvts() 
{ 
	tree1 = 0; 
	if (tree2)
	{
		delete tree2; 
		tree2=0;
	}
}

TTree* susyEvts::makeTree(TString treename)
{
	tree2 = new TTree(treename, "a angles tree");
	tree2->Branch("evt", &evt, EVT_s.c_str());
	tree2->Branch("leps", &vleps);
	tree2->Branch("l12", &l12, R_PAR_s.c_str());
	tree2->Branch("jets", &vjets);
	tree2->Branch("truths", &vtruths);
	tree2->Branch("sig", &sig, SIGNATURE_s.c_str());
	return tree2;
}

TTree* susyEvts::makeWeightOnlyTree(TString treename, susyEvts* parent)
{
	tree2 = new TTree(treename, "a angles tree");
	tree2->Branch("evt", &evt, EVT_s.c_str());

	if (parent->tree2)
	{ 
		tree2->SetDirectory(parent->tree2->GetDirectory());
		tree2->AddFriend( parent->tree2 );
	}
	else if (parent->tree1)
	{ 
		tree2->SetDirectory(parent->tree1->GetDirectory());
		tree2->AddFriend( parent->tree1 );
	} 

	return tree2;
}

TTree* susyEvts::makePtCorrTree(TString treename, susyEvts* parent)
{
	tree2 = new TTree(treename, "a angles tree");
	tree2->Branch("evt", &evt, EVT_s.c_str());
	tree2->Branch("leps", &vleps);
	tree2->Branch("l12", &l12, R_PAR_s.c_str());
	tree2->Branch("sig", &sig, SIGNATURE_s.c_str());

	if (parent->tree2)
	{ 
		tree2->SetDirectory(parent->tree2->GetDirectory());
		tree2->AddFriend( parent->tree2 );
	}
	else if (parent->tree1)
	{ 
		tree2->SetDirectory(parent->tree1->GetDirectory());
		tree2->AddFriend( parent->tree1 );
	} 

	return tree2;
}

TTree* susyEvts::makeKinematicsSysTree(TString treename, susyEvts* parent)
{
	// Make a Tree with every field
	// except the truth branch contents are identical to the parent NominalTree
	tree2 = new TTree(treename, "a angles tree");
	tree2->Branch("evt", &evt, EVT_s.c_str());
	tree2->Branch("leps", &vleps);
	tree2->Branch("l12", &l12, R_PAR_s.c_str());
	tree2->Branch("jets", &vjets);
	tree2->Branch("sig", &sig, SIGNATURE_s.c_str());

	if (parent->tree2)
	{ 
		tree2->SetDirectory(parent->tree2->GetDirectory());
		tree2->AddFriend( parent->tree2 );
	}
	else if (parent->tree1)
	{ 
		tree2->SetDirectory(parent->tree1->GetDirectory());
		tree2->AddFriend( parent->tree1 );
	} 

	return tree2;
}

void susyEvts::getTree(TTree* tr)
{
	tree1 = tr;
	tree1->SetBranchAddress("evt", (ULong64_t*)&evt); 
	tree1->SetBranchAddress("leps", &vleps);
	tree1->SetBranchAddress("l12", (float*)&l12);
	tree1->SetBranchAddress("jets", &vjets);
	tree1->SetBranchAddress("truths", &vtruths);
	tree1->SetBranchAddress("sig", (ULong64_t*)&sig);
}

Int_t susyEvts::Next()
{
	if (m_el >= tree1->GetEntriesFast()) return -1;
	int readBytes = tree1->GetEntry(m_el);
	m_el++;
	return readBytes;
}
