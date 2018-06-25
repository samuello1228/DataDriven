#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

#include <TFile.h>
#include <TChain.h>
#include <TBranch.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TString.h>

const double ETA_EL[] = {0, 0.8, 1.37, 1.52, 2.00, 2.47};
const double ETA_MU[] = {0, 0.6, 1.2, 1.8, 2.5};
const double PT[]  = {20, 30, 40, 50, 60, 70, 80, 90, 100, 200};
const unsigned int NETA_EL = sizeof(ETA_EL) / sizeof(ETA_EL[0]) - 1;
const unsigned int NETA_MU = sizeof(ETA_MU) / sizeof(ETA_MU[0]) - 1;
const unsigned int NPT     = sizeof(PT)  / sizeof(PT[0]) - 1;

void Draw(TString treeName,unsigned int NETA)
{
    TH1D* h2[NETA+1];
    for(unsigned int j=1;j<=NETA;j++)
    {
        TString Name = "eff";
        Name += TString::Itoa(j,10);
        h2[j] = new TH1D(Name.Data(),Name.Data(),NPT, PT);
    }
    
    TFile* file = new TFile("fake_rate.root","READ");
    TH2D *h1 = (TH2D*) file->Get(treeName.Data());
    for(unsigned int i=1;i<=NPT;i++)
    {
        for(unsigned int j=1;j<=NETA;j++)
        {
            cout<<h1->GetBinContent(i,j)<<", "<<h1->GetBinError(i,j)<<endl;
            h2[j]->SetBinContent(i,h1->GetBinContent(i,j));
            h2[j]->SetBinError(i,h1->GetBinError(i,j));
        }
        cout<<endl;
    }
    delete file;
    
    h2[1]->SetLineColor(kBlack);
    h2[2]->SetLineColor(kRed);
    h2[3]->SetLineColor(kGreen);
    h2[4]->SetLineColor(kBlue);
    
    h2[1]->GetYaxis()->SetRangeUser(0.6,1.2);
    
    TCanvas* c2 = new TCanvas();
    c2->cd();
    gStyle->SetOptStat(0);
    
    h2[1]->Draw();
    h2[2]->Draw("same");
    h2[3]->Draw("same");
    h2[4]->Draw("same");
    
    TString NameTemp = treeName;
    NameTemp += ".eps";
    c2->Print(NameTemp.Data(),"eps");
    
    for(unsigned int j=1;j<=NETA;j++)
    {
        delete h2[j];
    }
    delete c2;
}

void plot()
{
    Draw("El_hEff",NETA_EL);
    Draw("Mu_hEff",NETA_MU);
}


