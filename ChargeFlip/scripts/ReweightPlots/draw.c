// Author    : mat@ihep.ac.cn
// Date      : May 13, 2018 
// Descrption: It is used to draw some plots of reweighted Z->ee control samples for validation
// Inputs    : root files generated by program ../../util/re_weight.cxx 
// Outputs   : eps files store at ./plots 
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>
#include "Rtypes.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TROOT.h"

#include "AtlasLabels.C"
#include "AtlasStyle.C"

using namespace std;
#define PRINT(x) {std::cout << #x << " = " << x << std::endl;}

void DrawDataMC(TCanvas *c, 
		bool isLogy,
		TH1D *h_data,
		TH1D *h_mc, 
		TGraphAsymmErrors *g_mc,
		TGraphAsymmErrors *g_ratio,
		TLegend *leg,
		string leg_title,
		string x_title,
		string y_title = "Events",
		string z_title = "#frac{#font[52]{N}_{Weighted OS}}{#font[52]{N}_{Observed SS}}",
		string leg_entry_1 = "Observed SS",
		string leg_entry_2 = "Weighted OS",
		string leg_format_1 = "lep",
		string leg_format_2 = "fl");
void SetStyle();
int drawPlots(bool isData);

const double VETAS[] = {0, 0.5, 1, 1.37, 1.52, 1.8, 2.0, 2.47};
const double VPTS[]  = {25, 60, 90, 130, 150, 1000};
const int NETA  = sizeof(VETAS)/sizeof(VETAS[0]) - 1;
const int NPT   = sizeof(VPTS)/sizeof(VPTS[0]) - 1;
const int SIZE  = NETA*NPT;

//double LEG_LEFT_X  = 0.65;
//double LEG_LEFT_Y  = 0.75;
double LEG_LEFT_X  = 0.75;
double LEG_LEFT_Y  = 0.80;
double LEG_RIGHT_X = 0.85;
double LEG_RIGHT_Y = 0.90;
int    LINEWIDTH   = 2;

int draw()
{
	SetStyle();
	//drawPlots(true);
	drawPlots(false);

	return 0;
}

int drawPlots(bool isData)
{
	string file, output;
	if (isData)
	{
		file   = "../../run/checks_weight/data_80_100_20_20.root";
		output = "./plots/data_";
	}
	else 
	{
		file   = "../../run/checks_weight/mc_80_100_20_20.root";
		output = "./plots/mc_";
	}
	TFile *f = new TFile(file.c_str());
	if (!f->IsOpen())
	{
		cerr << "Cannot open file : " << file << endl;
		return -1;
	}
	
	struct VariableData
	{
		VariableData(string title, string name, bool islog = true)
		{
			x_title = title;
			o_name  = name;
			isLogy  = islog;
		}
		TH1D *h_data;
		TH1D *h_mc;
		TH1D *h_mc_down;
		TH1D *h_mc_up;
		TGraphAsymmErrors *g_mc;
		string x_title;
		string o_name;
		bool isLogy;
	};
	
	vector<VariableData> hist;
	hist.clear();
	hist.push_back( VariableData("m_{ll} (GeV)",                      "mll"));
	hist.push_back( VariableData("leading #font[52]{p}_{T} (GeV)",    "pt_1"));
	hist.push_back( VariableData("subleading #font[52]{p}_{T} (GeV)", "pt_2"));
	hist.push_back( VariableData("leading #eta",                      "eta_1", false));
	hist.push_back( VariableData("subleading #eta",                   "eta_2", false));
	hist.push_back( VariableData("leading #phi",                      "phi_1", false));
	hist.push_back( VariableData("subleading #phi",                   "phi_2", false));
	hist.push_back( VariableData("m_{ll} (GeV)",                      "mll_unweighted"));

	for (unsigned int i = 0; i < hist.size(); i++)
	{
		TString hName = "h_" + hist[i].o_name;
		hist[i].h_mc = (TH1D*) f->Get(hName.Data());
		
		TString hName2 = hName + "_down";
		hist[i].h_mc_down = (TH1D*) f->Get(hName2.Data());
		
		hName2 = hName + "_up";
		hist[i].h_mc_up = (TH1D*) f->Get(hName2.Data());
		
		hName2 = hName + "_ss";
		hist[i].h_data = (TH1D*) f->Get(hName2.Data());
		
		if(hist[i].o_name != "mll_unweighted")
		{
		int nbin = hist[i].h_mc->GetNbinsX();
		hist[i].g_mc = new TGraphAsymmErrors(nbin);
		
		TGraphAsymmErrors* g_ratio = new TGraphAsymmErrors(nbin);
		
		for (int j = 0; j < nbin; j++)
		{
			//set point
			double x = hist[i].h_mc->GetBinCenter(j+1);
			double y = hist[i].h_mc->GetBinContent(j+1);
			hist[i].g_mc->SetPoint(j,x,y);
			
			//set point error
			double xerr_down = hist[i].h_mc->GetBinWidth(j+1) /2;
			double xerr_up = hist[i].h_mc->GetBinWidth(j+1) /2;
			
			double y_down = hist[i].h_mc_down->GetBinContent(j+1);
			double y_up = hist[i].h_mc_up->GetBinContent(j+1);
			
			if(0<=y_down && y_down <= y && y <= y_up)
			{
			}
			else if(y_down > y && y > y_up && y_up<0)
			{
				y_down = hist[i].h_mc_up->GetBinContent(j+1);
				y_up = hist[i].h_mc_down->GetBinContent(j+1);
			}
			else
			{
				cout<<"Error!!!!!!!"<<endl;
				cout<<y_down<<", "<<y<<", "<<y_up<<endl;
			}
			
			double yerr_down = y - y_down;
			double yerr_up = y_up - y;
			//cout<<xerr_down<<", "<<x<<", "<<xerr_up<<"; "<<yerr_down<<", "<<y<<", "<<yerr_up<<endl;
			hist[i].g_mc->SetPointError(j,xerr_down,xerr_up,yerr_down,yerr_up);
			
			//set ratio plot
			double y_ss = hist[i].h_data->GetBinContent(j+1);
			double y_ss_error = hist[i].h_data->GetBinError(j+1);
			double y_ss_down = y_ss - y_ss_error;
			double y_ss_up = y_ss + y_ss_error;
			
			if(y_down>0 && y_ss_down>0)
			{
				double ratio = y/y_ss;
				g_ratio->SetPoint(j,x,ratio);
				
				double ratio_down = y_down/y_ss_up;
				double ratio_up = y_up/y_ss_down;
				g_ratio->SetPointError(j,xerr_down,xerr_up, ratio - ratio_down, ratio_up - ratio);
			}
			else
			{
				g_ratio->SetPoint(j,x,-99);
				g_ratio->SetPointError(j,xerr_down,xerr_up,0,0);
			}
		}
		
		TCanvas *c2 = new TCanvas(Form("c%d%d", isData, i), "c", 900, 900);
		TLegend *leg = new TLegend(LEG_LEFT_X, LEG_LEFT_Y, LEG_RIGHT_X, LEG_RIGHT_Y);
		string leg_title;
	   	if (isData)
		{
			leg_title = "Data";
		}
		else 
		{
			leg_title = "MC";
		}
		DrawDataMC(c2, hist[i].isLogy, hist[i].h_data, hist[i].h_mc, hist[i].g_mc, g_ratio, leg, leg_title, hist[i].x_title);
		
		//c2->Print((output + hist[i].o_name + ".eps").c_str(),"eps");
		c2->Print((output + hist[i].o_name + ".pdf").c_str(),"pdf");
		delete g_ratio;
		delete c2;
		}
		else
		{
			hist[i].h_mc->SetMinimum(10);
			hist[i].h_mc->SetMaximum(10000000);
			hist[i].h_mc->SetLineColor(kGreen);
			hist[i].h_data->SetLineColor(kBlue);
			
			hist[i].h_mc->GetYaxis()->SetTitle("Events");
			hist[i].h_mc->GetXaxis()->SetTitle(hist[i].x_title.c_str());
			
			TCanvas *c2 = new TCanvas();
			TPad* pad1 = new TPad("pad1","pad1",0,0,1,1);
			pad1->SetLogy(true);
			
			TLegend *leg = new TLegend(0.65, 0.7, 0.75, 0.8);
			leg->SetTextFont(42);
			leg->SetTextSize(0.04);
			leg->SetHeader("Data");
			leg->AddEntry(hist[i].h_mc,   "Observed OS");
			leg->AddEntry(hist[i].h_data, "Observed SS");
			
			c2->cd();
			pad1->Draw();
			pad1->cd();
			hist[i].h_mc->Draw("hist e2");
			hist[i].h_data->Draw("hist e2 same");
			leg->Draw();
			
			//ATLAS Label
			ATLASLabel(0.55,0.88,"Work in progress");
			TLatex lt1;
			lt1.DrawLatexNDC(0.55,0.83,"#sqrt{#it{s}} = 13 TeV, 36.1 fb^{-1}");
			
			c2->Print((output + hist[i].o_name + ".pdf").c_str(),"pdf");
			
			delete leg;
			delete pad1;
			delete c2;
		}
	}

	return 0;
}

void DrawDataMC(TCanvas *c, 
		bool isLogy,
		TH1D *h_data,
		TH1D *h_mc, 
		TGraphAsymmErrors *g_mc,
		TGraphAsymmErrors *g_ratio,
		TLegend *leg,
		string leg_title,
		string x_title,
		string y_title,
		string z_title,
		string leg_entry_1,
		string leg_entry_2,
		string leg_format_1,
		string leg_format_2)
{
	char pad_name[100];
	sprintf(pad_name, "%s_%d", c->GetName(), 0);

	double size = 0.25;

	TPad *pad_1 = new TPad(pad_name,"This is pad1", 0.0, size, 1, 1);
	sprintf(pad_name, "%s_%d", c->GetName(), 1);
	TPad *pad_2 = new TPad(pad_name,"This is pad2", 0.0, 0.0, 1, size);

	pad_1->SetFillColor(0);
	pad_1->SetLeftMargin(0.23);
	pad_1->SetBottomMargin(0.17);
	pad_2->SetFillColor(0);
	pad_2->SetLeftMargin(0.23);
	pad_2->SetBottomMargin(0.17);
	pad_1->Draw();
	pad_2->Draw();
	if (isLogy)
	{
		pad_1->SetLogy();
	}
	pad_1->cd();

	h_data->SetMarkerSize(0.9);
	h_data->SetMarkerStyle(20);
	h_data->SetMarkerColor( kBlack );
	h_data->SetLineWidth(LINEWIDTH);
	h_data->SetLineColor(kBlack);
	h_data->Draw("E");
	
	double x_label_size = h_data->GetXaxis()->GetLabelSize();
	double y_label_size = h_data->GetYaxis()->GetLabelSize();
	h_data->GetXaxis()->SetTitle( x_title.c_str() );
	h_data->GetXaxis()->CenterTitle();
	h_data->GetYaxis()->SetTitle( y_title.c_str() );
	h_data->GetYaxis()->CenterTitle();
	double max = (h_data->GetMaximum() > h_mc->GetMaximum()) ? h_data->GetMaximum() : h_mc->GetMaximum();
	double min = (h_data->GetMinimum() < h_mc->GetMinimum()) ? h_data->GetMinimum() : h_mc->GetMinimum();
	if (min < 1e-5)
	{
		min = 1.;
	}
	double x_tick_length = h_data->GetXaxis()->GetTickLength();
	double font_size     = h_data->GetYaxis()->GetTitleSize();
	h_data->GetYaxis()->SetRangeUser(0.8 * min, 1.5 * max);
	h_mc->SetLineColor(kBlue);
	h_mc->Draw("HIST, SAME");
	h_data->Draw("E, SAME");
	
	g_mc->SetFillColor(kBlue);
	g_mc->SetFillStyle(3004);
	g_mc->SetLineColor(kBlue);
	g_mc->Draw("e2, SAME");

	pad_1->RedrawAxis(); //force the axis redrawing

	leg->SetTextFont(42);
	//leg->SetTextSize(0.05);
	leg->SetTextSize(0.04);
	//leg->SetHeader(leg_title.c_str(), "c");
	leg->SetHeader(leg_title.c_str());
	leg->AddEntry(h_data, leg_entry_1.c_str(), leg_format_1.c_str());
	leg->AddEntry(g_mc,   leg_entry_2.c_str(), leg_format_2.c_str());
	leg->Draw();
	
	//ATLAS Label
	ATLASLabel(0.01,0.96,"Work in progress");
	TLatex lt1;
	lt1.DrawLatexNDC(0.54,0.955,"#sqrt{#it{s}} = 13 TeV, 36.1 fb^{-1}");
	
	pad_1->Update();

	pad_2->cd();
	pad_2->SetGridy();
	TH1D *h_ratio = (TH1D*)h_mc->Clone(Form("%s+%s", h_data->GetName(), h_mc->GetName()));
	h_ratio->Divide(h_data);

	h_ratio->SetLineColor(kBlue);
	h_ratio->Draw("HIST");
	h_ratio->GetYaxis()->SetTitle( z_title.c_str() );
	h_ratio->GetYaxis()->CenterTitle();
	h_ratio->GetXaxis()->SetTitle(0);
	h_ratio->GetYaxis()->SetTitleSize( font_size * (1-size) / size * 0.8 );
	h_ratio->GetYaxis()->SetTitleOffset( 0.4 );
	h_ratio->GetYaxis()->SetRangeUser(0.0, 3.0);
	h_ratio->GetYaxis()->SetNdivisions(505);
	h_ratio->GetXaxis()->SetLabelSize( x_label_size * (1-size) / size );
	h_ratio->GetYaxis()->SetLabelSize( y_label_size * (1-size) / size );
	h_ratio->GetXaxis()->SetTickLength( x_tick_length * (1-size) / size );
	
	g_ratio->SetFillColor(kBlue);
	g_ratio->SetFillStyle(3004);
	g_ratio->SetLineColor(kBlue);
	g_ratio->Draw("e2, SAME");

	pad_2->Update();
}

void SetStyle()
{
	gStyle->SetOptStat(kFALSE);
	// No Canvas Border
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetCanvasBorderSize(0);
	// White BG
	gStyle->SetCanvasColor(10);
	// Format for axes
	gStyle->SetLabelFont(42,"xyz");
	gStyle->SetLabelSize(0.05,"xyz");
	gStyle->SetLabelOffset(0.01,"xyz");
	gStyle->SetNdivisions(510,"xyz");
	gStyle->SetTitleFont(42,"xyz");
	gStyle->SetTitleColor(1,"xyz");
	gStyle->SetTitleSize(0.07,"xyz");
	gStyle->SetTitleOffset(1.0,"x");
	gStyle->SetTitleOffset(1.0,"y");

	//gStyle->SetTitleOffset(1.15,"xyz");
	// No pad borders
	gStyle->SetPadBorderMode(0);
	gStyle->SetPadBorderSize(0);
	// White BG
	gStyle->SetPadColor(10);
	// Margins for labels etc.
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetPadBottomMargin(0.17);
	gStyle->SetPadTopMargin(0.05);
	// No error bars in x direction
	gStyle->SetErrorX(0);
	// Format legend
	gStyle->SetLegendBorderSize(0);
	gStyle->SetOptTitle(0);
	gStyle->SetEndErrorSize(0);
	gROOT->ForceStyle();
}
