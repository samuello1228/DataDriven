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

const double ETA_EL[] = {0, 1.37, 1.52, 2.47};
const double ETA_MU[] = {0, 1.37, 1.52, 2.4};
const double PT_EL[]  = {25, 35, 45, 120, 200};
const double PT_MU[]  = {25, 30, 45, 120, 200};
const unsigned int NETA_EL = sizeof(ETA_EL) / sizeof(ETA_EL[0]) - 1;
const unsigned int NETA_MU = sizeof(ETA_MU) / sizeof(ETA_MU[0]) - 1;
const unsigned int NPT_EL  = sizeof(PT_EL)  / sizeof(PT_EL[0])  - 1;
const unsigned int NPT_MU  = sizeof(PT_MU)  / sizeof(PT_MU[0])  - 1;

//Peter result
const double Params_Fake_el_nominal[] = {
    0.115837, 0.147943, 0.21007, 0.177491, 0.0621802, 0.0566687, 0.137219, 0.126778, 0.13501, 0.0963788, 0.0909373, 0.121558
};

const double Params_Fake_el_statDOWNdiff[] = {
    0.0181598, 0.0300393, 0.0274507, 0.0438151, 0.029465, 0.0375158, 0.0326709, 0.0248499, 0.0246109, 0.0257812, 0.0136109, 0.036831
};

const double Params_Fake_mu_nominal[] = {
    0.175618, 0.127712, 0.208001, 0.118104, 0.24175, 0.0582174, 0.217779, 0.142331, 0.127783, 0.152324, 0.229239, 0.392727
};

const double Params_Fake_mu_statDOWNdiff[] = {
    0.0376544, 0.0361995, 0.0540197, 0.0446926, 0.110971, 0.0637628, 0.275801, 0.027893, 0.0384687, 0.0513694, 0.120531, 0.538107
};

void GetEffs(TString lepton_type,unsigned int NETA,const double ETA[],unsigned int NPT,const double PT[],TH2D* hEff)
{
    TFile* data = new TFile("../../run/output/fr_mxm/data_fr_mxm.root","READ");
    TH2D *hTight_data = (TH2D*) data->Get(lepton_type+"_hTight");
    TH2D *hLoose_data = (TH2D*) data->Get(lepton_type+"_hLoose");
    
    TFile* mc = new TFile("../../run/output/fr_mxm/mc_fr_mxm.root","READ");
    TH2D *hTight_prompt_mc = (TH2D*) mc->Get(lepton_type+"_prompt_hTight");
    TH2D *hLoose_prompt_mc = (TH2D*) mc->Get(lepton_type+"_prompt_hLoose");
    
    for(unsigned int i=1;i<=NPT;i++)
    {
        cout<<"For "<<PT[i-1]<<" < pt < "<<PT[i]<<" :"<<endl;
        for(unsigned int j=1;j<=NETA;j++)
        {
            cout<<"For "<<ETA[j-1]<<" < |eta| < "<<ETA[j]<<" : "<<endl;
            
            cout<<"N^{data}_{signal}: "<<hTight_data->GetBinContent(i,j)<<" +/- "<<hTight_data->GetBinError(i,j)<<endl;
            cout<<"N^{prompt bkg}_{signal}: "<<hTight_prompt_mc->GetBinContent(i,j)<<" +/- "<<hTight_prompt_mc->GetBinError(i,j)<<endl;
            cout<<"N^{data}_{baseline}: "<<hLoose_data->GetBinContent(i,j)<<" +/- "<<hLoose_data->GetBinError(i,j)<<endl;
            cout<<"N^{prompt bkg}_{baseline}: "<<hLoose_prompt_mc->GetBinContent(i,j)<<" +/- "<<hLoose_prompt_mc->GetBinError(i,j)<<endl;
            
            double n = hTight_data->GetBinContent(i,j) - hTight_prompt_mc->GetBinContent(i,j);
            double ne = sqrt(hTight_data->GetBinError(i,j) * hTight_data->GetBinError(i,j)
                           + hTight_prompt_mc->GetBinError(i,j) * hTight_prompt_mc->GetBinError(i,j));
            double d = hLoose_data->GetBinContent(i,j) - hLoose_prompt_mc->GetBinContent(i,j);
            double de = sqrt(hLoose_data->GetBinError(i,j) * hLoose_data->GetBinError(i,j)
                           + hLoose_prompt_mc->GetBinError(i,j) * hLoose_prompt_mc->GetBinError(i,j));
            
            double eff = n/d;
            double effe = eff * sqrt( (ne/n)*(ne/n) + (de/d)*(de/d) );
            cout<<eff<<" +/- "<<effe<<endl;
        }
        cout<<endl;
    }
    
    //For hEff
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
        
        //fill for pt overflow
        for(unsigned int j=1;j<=NETA;j++)
        {
            hEff->SetBinContent(NPT+1,j,hEff->GetBinContent(NPT,j));
            hEff->SetBinError(NPT+1,j,hEff->GetBinError(NPT,j));
        }
    }
    
    //hEff_eta
    TH1D* hEff_eta = new TH1D(lepton_type+"_hEff_eta", ";|#eta|", NETA, ETA);
    {
        TH2D *hTight_data_eta = (TH2D*) hTight_data->RebinX(NPT,lepton_type+"_hTight_eta");
        TH2D *hLoose_data_eta = (TH2D*) hLoose_data->RebinX(NPT,lepton_type+"_hLoose_eta");
        TH2D *hTight_prompt_mc_eta = (TH2D*) hTight_prompt_mc->RebinX(NPT,lepton_type+"_prompt_hTight_eta");
        TH2D *hLoose_prompt_mc_eta = (TH2D*) hLoose_prompt_mc->RebinX(NPT,lepton_type+"_prompt_hLoose_eta");
        
        //hEff = numerator = hTight_data - hTight_prompt_mc
        hTight_prompt_mc_eta->Scale(-1);
        hTight_data_eta->Add(hTight_prompt_mc_eta);
        
        //denominator = hLoose_data - hLoose_prompt_mc
        hLoose_prompt_mc_eta->Scale(-1);
        hLoose_data_eta->Add(hLoose_prompt_mc_eta);
        
        //hEff = hEff / denominator
        hTight_data_eta->Divide(hLoose_data_eta);
        
        for(unsigned int i=1;i<=NETA;i++)
        {
            cout<<hTight_data_eta->GetBinContent(1,i)<<", "<<hTight_data_eta->GetBinError(1,i)<<endl;
            hEff_eta->SetBinContent(i,hTight_data_eta->GetBinContent(1,i));
            hEff_eta->SetBinError(i,hTight_data_eta->GetBinError(1,i));
        }
        cout<<endl;
        
        delete hTight_data_eta;
        delete hLoose_data_eta;
        delete hTight_prompt_mc_eta;
        delete hLoose_prompt_mc_eta;
    }
    
    //For hEff_pt
    TH1D* hEff_pt = new TH1D(lepton_type+"_hEff_pt", ";p_{T} [GeV]", NPT, PT);
    {
        TH2D *hTight_data_pt = (TH2D*) hTight_data->RebinY(NETA,lepton_type+"_hTight_pt");
        TH2D *hLoose_data_pt = (TH2D*) hLoose_data->RebinY(NETA,lepton_type+"_hLoose_pt");
        TH2D *hTight_prompt_mc_pt = (TH2D*) hTight_prompt_mc->RebinY(NETA,lepton_type+"_prompt_hTight_pt");
        TH2D *hLoose_prompt_mc_pt = (TH2D*) hLoose_prompt_mc->RebinY(NETA,lepton_type+"_prompt_hLoose_pt");
        
        //hEff = numerator = hTight_data - hTight_prompt_mc
        hTight_prompt_mc_pt->Scale(-1);
        hTight_data_pt->Add(hTight_prompt_mc_pt);
        
        //denominator = hLoose_data - hLoose_prompt_mc
        hLoose_prompt_mc_pt->Scale(-1);
        hLoose_data_pt->Add(hLoose_prompt_mc_pt);
        
        //hEff = hEff / denominator
        hTight_data_pt->Divide(hLoose_data_pt);
        
        for(unsigned int i=1;i<=NPT;i++)
        {
            cout<<hTight_data_pt->GetBinContent(i,1)<<", "<<hTight_data_pt->GetBinError(i,1)<<endl;
            hEff_pt->SetBinContent(i,hTight_data_pt->GetBinContent(i,1));
            hEff_pt->SetBinError(i,hTight_data_pt->GetBinError(i,1));
        }
        cout<<endl;
        
        delete hTight_data_pt;
        delete hLoose_data_pt;
        delete hTight_prompt_mc_pt;
        delete hLoose_prompt_mc_pt;
    }
    
    //output result for latex
    {
        double Params_Fake_nominal[NETA*NPT];
        double Params_Fake_statDOWNdiff[NETA*NPT];
        if(lepton_type=="El")
        {
            for(unsigned int i=0;i<NETA*NPT;i++) Params_Fake_nominal[i] = Params_Fake_el_nominal[i];
            for(unsigned int i=0;i<NETA*NPT;i++) Params_Fake_statDOWNdiff[i] = Params_Fake_el_statDOWNdiff[i];
        }
        else if(lepton_type=="Mu")
        {
            for(unsigned int i=0;i<NETA*NPT;i++) Params_Fake_nominal[i] = Params_Fake_mu_nominal[i];
            for(unsigned int i=0;i<NETA*NPT;i++) Params_Fake_statDOWNdiff[i] = Params_Fake_mu_statDOWNdiff[i];
        }
        
        TString PathName = "fake_";
        PathName += lepton_type;
        PathName += ".tex";
        
        ofstream fout;
        fout.open(PathName.Data());
        
        fout<<"\\begin{tabular}{|";
        for(unsigned int i=1;i<=NETA+1;i++) fout<<"c|";
        fout<<"}"<<endl;
        fout<<"\\hline"<<endl;
        
        for(unsigned int i=1;i<=NETA;i++) fout<<" & $"<<ETA[i-1]<<" < |\\eta| < "<<ETA[i]<<"$";
        fout<<" \\\\"<<endl;
        fout<<"\\hline"<<endl;
        
        for(unsigned int i=1;i<=NPT;i++)
        {
            //For my result
            fout<<"$"<<PT[i-1]<<" < p_T < "<<PT[i]<<"$";
            for(unsigned int j=1;j<=NETA;j++)
            {
                fout<<" & \\color{orange} $"<<hEff->GetBinContent(i,j)<<" \\pm "<<hEff->GetBinError(i,j)<<"$";
            }
            fout<<" \\\\"<<endl;
            
            //For Peter result
            for(unsigned int j=1;j<=NETA;j++)
            {
                int index = (j-1) * NPT + (i-1);
                fout<<" & \\color{blue} $"<<Params_Fake_nominal[index]<<" \\pm "<<Params_Fake_statDOWNdiff[index]<<"$";
            }
            fout<<" \\\\"<<endl;
            fout<<"\\hline"<<endl;
        }
        
        fout<<"\\end{tabular}";
        fout.close();
    }
    
    //For h2
    TH1D* h2[NETA+1];
    for(unsigned int j=1;j<=NETA;j++)
    {
        TString Name = "eff";
        Name += TString::Itoa(j,10);
        h2[j] = new TH1D(Name.Data(),Name.Data(),NPT, PT);
    }
    
    for(unsigned int i=1;i<=NPT;i++)
    {
        cout<<"For "<<PT[i-1]<<" < pt < "<<PT[i]<<" :"<<endl;
        for(unsigned int j=1;j<=NETA;j++)
        {
            cout<<"For "<<ETA[j-1]<<" < |eta| < "<<ETA[j]<<" : ";
            cout<<hEff->GetBinContent(i,j)<<" +/- "<<hEff->GetBinError(i,j)<<endl;
            h2[j]->SetBinContent(i,hEff->GetBinContent(i,j));
            h2[j]->SetBinError(i,hEff->GetBinError(i,j));
        }
        cout<<endl;
    }
    
    SetAtlasStyle();
    
    h2[1]->SetLineColor(kBlue);
    h2[2]->SetLineColor(kRed);
    h2[3]->SetLineColor(kGreen);
    
    h2[1]->SetMarkerStyle(20);
    h2[2]->SetMarkerStyle(21);
    h2[3]->SetMarkerStyle(22);
    
    h2[1]->SetMarkerColor(kBlue);
    h2[2]->SetMarkerColor(kRed);
    h2[3]->SetMarkerColor(kGreen);
    
    h2[1]->GetXaxis()->SetTitle("p_{T} [GeV]");
    h2[1]->GetYaxis()->SetRangeUser(-1,2);
    
    if(lepton_type=="El") h2[1]->GetYaxis()->SetTitle("Electron Fake Efficiency");
    else if(lepton_type=="Mu") h2[1]->GetYaxis()->SetTitle("Muon Fake Efficiency");
    
    TCanvas* c2 = new TCanvas();
    c2->cd();
    gStyle->SetOptStat(0);
    
    h2[1]->Draw();
    for(unsigned int j=2;j<=NETA;j++)
    {
        h2[j]->Draw("same");
    }
    
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
        Name = TString::Format("%.2f",ETA[j-1]);
        Name += "<|#eta|<";
        Name += TString::Format("%.2f",ETA[j]);
        
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
    
    //For hEff_eta and hEff_pt
    hEff_eta->SetLineColor(kBlack);
    hEff_eta->SetMarkerStyle(20);
    hEff_eta->SetMarkerColor(kBlack);
    
    hEff_pt->SetLineColor(kBlack);
    hEff_pt->SetMarkerStyle(20);
    hEff_pt->SetMarkerColor(kBlack);
    
    hEff_eta->GetXaxis()->SetTitle("|#eta|");
    hEff_eta->GetYaxis()->SetRangeUser(0,1);
    
    hEff_pt->GetXaxis()->SetTitle("p_{T} [GeV]");
    hEff_pt->GetYaxis()->SetRangeUser(0,1);
    
    if(lepton_type=="El")
    {
        hEff_eta->GetYaxis()->SetTitle("Electron Fake Efficiency");
        hEff_pt->GetYaxis()->SetTitle("Electron Fake Efficiency");
    }
    else if(lepton_type=="Mu")
    {
        hEff_eta->GetYaxis()->SetTitle("Muon Fake Efficiency");
        hEff_pt->GetYaxis()->SetTitle("Muon Fake Efficiency");
    }
    
    c2 = new TCanvas();
    c2->cd();
    hEff_eta->Draw();
    ATLASLabel(0.2,0.88,"Internal");
    NameTemp = lepton_type + "_hEff_eta";
    NameTemp += ".eps";
    c2->Print(NameTemp.Data(),"eps");
    delete hEff_eta;
    delete c2;
    
    c2 = new TCanvas();
    c2->cd();
    hEff_pt->Draw();
    ATLASLabel(0.2,0.88,"Internal");
    NameTemp = lepton_type + "_hEff_pt";
    NameTemp += ".eps";
    c2->Print(NameTemp.Data(),"eps");
    delete hEff_pt;
    delete c2;
    
    delete data;
    delete mc;
}

void getFakeEffs()
{
    TFile* fake_eff = new TFile("fake_eff.root","RECREATE");
    TH2D* El_hEff = new TH2D("El_hEff", ";p_{T} [GeV];|#eta|", NPT_EL, PT_EL, NETA_EL, ETA_EL);
    GetEffs("El", NETA_EL, ETA_EL, NPT_EL, PT_EL, El_hEff);
    fake_eff->cd();
    El_hEff->Write();
    
    TH2D* Mu_hEff = new TH2D("Mu_hEff", ";p_{T} [GeV];|#eta|", NPT_MU, PT_MU, NETA_MU, ETA_MU);
    GetEffs("Mu", NETA_MU, ETA_MU, NPT_MU, PT_MU, Mu_hEff);
    fake_eff->cd();
    Mu_hEff->Write();
    
    delete fake_eff;
}
