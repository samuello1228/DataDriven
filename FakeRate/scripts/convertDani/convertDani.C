// C++
#include <iostream>
using namespace std;
// ROOT
#include "TChain.h"
#include "TString.h"
#include "TFile.h"

//CommonSelection
Int_t NlepBL;
Int_t nSigLep;
Int_t nBJets20;

//SRSelection
Float_t Pt_l;
Float_t Pt_subl;

//other selection
Int_t nJets20;
Float_t DeltaEtaLep;
Float_t met;
Float_t mt;
Float_t meff;
Float_t mljj_comb;
Float_t MT2;

Float_t FakeWeight;

void convertDani()
{
    //read from Peter fake ntuple
    const TString TREE_NAME = "Fakes_nom";
    TChain* tree1 = new TChain(TREE_NAME.Data());
    TString path1 = "/eos/user/d/dkoeck/WHSS/HistFitterTrees/Trees/Fakes_HF.root";
    tree1->Add(path1.Data());
    cout << "There are " << tree1->GetEntries() << " events" << endl;

    tree1->SetBranchAddress("NlepBL", &NlepBL);
    tree1->SetBranchAddress("nSigLep", &nSigLep);
    tree1->SetBranchAddress("nBJets20", &nBJets20);
    tree1->SetBranchAddress("Pt_l", &Pt_l);
    tree1->SetBranchAddress("Pt_subl", &Pt_subl);
    tree1->SetBranchAddress("nJets20", &nJets20);
    tree1->SetBranchAddress("DeltaEtaLep", &DeltaEtaLep);
    tree1->SetBranchAddress("met", &met);
    tree1->SetBranchAddress("mt", &mt);
    tree1->SetBranchAddress("meff", &meff);
    tree1->SetBranchAddress("mljj_comb", &mljj_comb);
    tree1->SetBranchAddress("MT2", &MT2);
    tree1->SetBranchAddress("FakeWeight", &FakeWeight);

    /*
    tree1->GetEntry(0);
    cout<<"nSigLep: "<<nSigLep<<endl;
    cout<<"NlepBL: "<<NlepBL<<endl;
    cout<<"nJets20: "<<nJets20<<endl;
    cout<<"nBJets20: "<<nBJets20<<endl;
    cout<<"FakeWeight: "<<FakeWeight<<endl;
    cout<<"met: "<<met<<endl;
    cout<<"meff: "<<meff<<endl;
    cout<<"mljj_comb: "<<mljj_comb<<endl;
    cout<<"DeltaEtaLep: "<<DeltaEtaLep<<endl;
    cout<<"mt: "<<mt<<endl;
    cout<<"MT2: "<<MT2<<endl;
    cout<<"Pt_l: "<<Pt_l<<endl;
    cout<<"Pt_subl: "<<Pt_subl<<endl;
    */

    //create new tree
    TString path2 = "Fakes_HF_converted.root";
    TFile* file2 = new TFile(path2.Data(),"RECREATE");
    TTree* tree2 = new TTree(TREE_NAME.Data(),TREE_NAME.Data());

    tree2->Branch("NlepBL", &NlepBL, "NlepBL/I");
    tree2->Branch("nSigLep", &nSigLep, "nSigLep/I");
    tree2->Branch("nBJets20", &nBJets20, "nBJets20/I");
    tree2->Branch("Pt_l", &Pt_l, "Pt_l/F");
    tree2->Branch("Pt_subl", &Pt_subl, "Pt_subl/F");
    tree2->Branch("nJets20", &nJets20, "nJets20/I");
    tree2->Branch("DeltaEtaLep", &DeltaEtaLep, "DeltaEtaLep/F");
    tree2->Branch("met", &met, "met/F");
    tree2->Branch("mt", &mt, "mt/F");
    tree2->Branch("meff", &meff, "meff/F");
    tree2->Branch("mljj_comb", &mljj_comb, "mljj_comb/F");
    tree2->Branch("MT2", &MT2, "MT2/F");
    tree2->Branch("FakeWeight", &FakeWeight, "FakeWeight/F");

    //copy
    for (long i = 0; i < tree1->GetEntries(); i++)
    {
        tree1->GetEntry(i);
        tree2->Fill();
    }

    file2->cd();
    tree2->Write();

    //delete
    delete tree1;
    delete tree2;
    delete file2;
}
