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
#include "AtlasLabels.C"
#include "AtlasStyle.C"

const double ETA_EL[] = {0, 0.8, 1.37, 1.52, 2.47};
const double ETA_MU[] = {0, 0.6, 1.2, 1.8, 2.4};
const double PT[]  = {25, 35, 45, 55, 65, 75, 85, 95};
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
    
    TFile* file = new TFile("../../run/output/tag_probe/fake_rate.root","READ");
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
    
    SetAtlasStyle();
    
    h2[1]->SetLineColor(kBlack);
    h2[2]->SetLineColor(kRed);
    h2[3]->SetLineColor(kGreen);
    h2[4]->SetLineColor(kBlue);
    
    h2[1]->SetMarkerStyle(20);
    h2[2]->SetMarkerStyle(21);
    h2[3]->SetMarkerStyle(22);
    h2[4]->SetMarkerStyle(23);
    
    h2[1]->SetMarkerColor(kBlack);
    h2[2]->SetMarkerColor(kRed);
    h2[3]->SetMarkerColor(kGreen);
    h2[4]->SetMarkerColor(kBlue);
    
    
    h2[1]->GetXaxis()->SetTitle("p_{T} [GeV]");
    h2[1]->GetYaxis()->SetRangeUser(0.6,1.1);
    
    if(treeName=="El_hEff") h2[1]->GetYaxis()->SetTitle("Electron Real Efficiency");
    else if(treeName=="Mu_hEff") h2[1]->GetYaxis()->SetTitle("Muon Real Efficiency");
    
    TCanvas* c2 = new TCanvas();
    c2->cd();
    gStyle->SetOptStat(0);
    
    h2[1]->Draw();
    h2[2]->Draw("same");
    h2[3]->Draw("same");
    h2[4]->Draw("same");
    
    {
        ATLASLabel(0.2,0.88,"Work in progress");
        TLatex lt2;
        TString NameTemp = "#sqrt{#it{s}} = 13 TeV, 36.1 fb^{-1}";
        lt2.DrawLatexNDC(0.2,0.83, NameTemp.Data());
    }
    
    Double_t xl1, yl1, xl2, yl2;
    xl2=0.85;
    yl2=0.4;
    xl1=xl2-0.2;
    yl1=yl2-0.2;
    
    TLegend* leg = new TLegend(xl1,yl1,xl2,yl2);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetBorderSize(0);
    
    for(unsigned int j=1;j<=4;j++)
    {
        TString Name;
        if(treeName=="El_hEff")
        {
            Name = TString::Format("%.2f",ETA_EL[j-1]);
            Name += "<|#eta|<";
            Name += TString::Format("%.2f",ETA_EL[j]);
        }
        else if(treeName=="Mu_hEff")
        {
            Name = TString::Format("%.2f",ETA_MU[j-1]);
            Name += "<|#eta|<";
            Name += TString::Format("%.2f",ETA_MU[j]);
        }
        
        leg->AddEntry(h2[j],Name.Data(),"pl");
    }
    
    leg->Draw();
    
    
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


