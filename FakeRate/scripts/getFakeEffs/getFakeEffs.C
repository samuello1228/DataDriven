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

const double ETA_EL[] = {0, 1.37, 2.47};
const double ETA_MU[] = {0, 2.47};
const double PT[]  = {20, 25, 30, 40, 200};
const unsigned int NETA_EL = sizeof(ETA_EL) / sizeof(ETA_EL[0]) - 1;
const unsigned int NETA_MU = sizeof(ETA_MU) / sizeof(ETA_MU[0]) - 1;
const unsigned int NPT     = sizeof(PT)  / sizeof(PT[0]) - 1;

void GetEffs(TString lepton_type,unsigned int NETA,TH2D* hEff)
{
    TFile* data = new TFile("../../run/output/fr_mxm/data_fr_mxm.root","READ");
    TH2D *hTight_data = (TH2D*) data->Get(lepton_type+"_hTight");
    TH2D *hLoose_data = (TH2D*) data->Get(lepton_type+"_hLoose");
    
    TFile* mc = new TFile("../../run/output/fr_mxm/mc_fr_mxm.root","READ");
    TH2D *hTight_prompt_mc = (TH2D*) mc->Get(lepton_type+"_prompt_hTight");
    TH2D *hLoose_prompt_mc = (TH2D*) mc->Get(lepton_type+"_prompt_hLoose");
    
    {
        //hEff = numerator = hTight_data - hTight_prompt_mc
        hEff->Add(hTight_data);
        TH2D *hTemp = (TH2D*) hTight_prompt_mc->Clone("hTemp");
        hTemp->Scale(-1);
        hEff->Add(hTemp);
        
        //hTemp = denominator = hLoose_data - hLoose_prompt_mc
        hTemp->Scale(0);
        hTemp->Add(hLoose_prompt_mc);
        hTemp->Scale(-1);
        hTemp->Add(hLoose_data);
        
        //hEff = hEff / denominator
        hEff->Divide(hTemp);
        
        delete hTemp;
    }
    
    delete data;
    delete mc;
    
    TH1D* h2[NETA+1];
    for(unsigned int j=1;j<=NETA;j++)
    {
        TString Name = "eff";
        Name += TString::Itoa(j,10);
        h2[j] = new TH1D(Name.Data(),Name.Data(),NPT, PT);
    }
    
    for(unsigned int i=1;i<=NPT;i++)
    {
        for(unsigned int j=1;j<=NETA;j++)
        {
            cout<<hEff->GetBinContent(i,j)<<", "<<hEff->GetBinError(i,j)<<endl;
            h2[j]->SetBinContent(i,hEff->GetBinContent(i,j));
            h2[j]->SetBinError(i,hEff->GetBinError(i,j));
        }
        cout<<endl;
    }
    
    SetAtlasStyle();
    
    h2[1]->SetLineColor(kBlue);
    if(lepton_type=="El") h2[2]->SetLineColor(kRed);
    
    h2[1]->SetMarkerStyle(20);
    if(lepton_type=="El") h2[2]->SetMarkerStyle(21);
    
    h2[1]->SetMarkerColor(kBlue);
    if(lepton_type=="El") h2[2]->SetMarkerColor(kRed);
    
    h2[1]->GetXaxis()->SetTitle("p_{T} [GeV]");
    h2[1]->GetYaxis()->SetRangeUser(0,0.5);
    
    if(lepton_type=="El") h2[1]->GetYaxis()->SetTitle("Electron Fake Efficiency");
    else if(lepton_type=="Mu") h2[1]->GetYaxis()->SetTitle("Muon Fake Efficiency");
    
    TCanvas* c2 = new TCanvas();
    c2->cd();
    gStyle->SetOptStat(0);
    
    h2[1]->Draw();
    if(lepton_type=="El") h2[2]->Draw("same");
    
    ATLASLabel(0.2,0.88,"Internal");
    
    Double_t xl1, yl1, xl2, yl2;
    xl2=0.85;
    yl2=0.85;
    xl1=xl2-0.2;
    yl1=yl2-0.2;
    
    TLegend* leg = new TLegend(xl1,yl1,xl2,yl2);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetBorderSize(0);
    
    for(unsigned int j=1;j<=NETA;j++)
    {
        TString Name;
        if(lepton_type=="El")
        {
            Name = TString::Format("%.2f",ETA_EL[j-1]);
            Name += "<|#eta|<";
            Name += TString::Format("%.2f",ETA_EL[j]);
        }
        else if(lepton_type=="Mu")
        {
            Name = TString::Format("%.2f",ETA_MU[j-1]);
            Name += "<|#eta|<";
            Name += TString::Format("%.2f",ETA_MU[j]);
        }
        
        leg->AddEntry(h2[j],Name.Data(),"pl");
    }
    
    leg->Draw();
    
    TString NameTemp = lepton_type + "_hEff";
    NameTemp += ".eps";
    c2->Print(NameTemp.Data(),"eps");
    
    for(unsigned int j=1;j<=NETA;j++)
    {
        delete h2[j];
    }
    delete c2;
}

void getFakeEffs()
{
    TFile* fake_eff = new TFile("fake_eff.root","RECREATE");
    TH2D* El_hEff = new TH2D("El_hEff", ";p_{T} [GeV];|#eta|", NPT, PT, NETA_EL, ETA_EL);
    GetEffs("El",NETA_EL,El_hEff);
    fake_eff->cd();
    El_hEff->Write();
    
    TH2D* Mu_hEff = new TH2D("Mu_hEff", ";p_{T} [GeV];|#eta|", NPT, PT, NETA_MU, ETA_MU);
    GetEffs("Mu",NETA_MU,Mu_hEff);
    fake_eff->cd();
    Mu_hEff->Write();
    
    delete fake_eff;
}
