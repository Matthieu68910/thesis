#include "TGraph.h"
#include "TGraphErrors.h"
#include <iostream>
#include <math.h>

/*
	Graph from Adam2020 fig.20
*/
Double_t myFit1(double x){
	double a, b, c, d;

	//cout << "1 selected" << endl;
	a = -0.01289257;
	b = 27.09822;
	c = 1.874107;
	d = 0.9877558;
	
	return d + (a - d) / (1 + pow((x / c), b));
}

Double_t myFit2(double x){
	double a, b, c, d;

	//cout << "2 selected" << endl;
	a = -0.01153673;
	b = 27.92117;
	c = 1.930889;
	d = 1.001919;
	
	return d + (a - d) / (1 + pow((x / c), b));
}

void Graph20(){
	auto c1 = new TCanvas("c1","c1",1920,1080);
	c1->SetTitle("Stub efficiency for 2S mini-module");
	//c1->UseCurrentStyle();
    //gStyle->SetOptStat(0);

    /*TStyle *st1 = new TStyle("st1","my style");
    st1->SetAxisColor(1, "xy");
    st1->SetGridColor(1);
    st1->SetGridStyle(2);
    st1->SetGridWidth(1);
    st1->SetLabelSize(0.03, "xy");
    st1->SetNdivisions(510, "xy");
    st1->SetTitleAlign(-10);
    st1->cd();  // this becomes now the current style gStyle

    gROOT->SetStyle("st1");*/

    gPad->SetGridx(1);
    gPad->SetGridy(1);
    gPad->SetTitle("Stub efficiency for 2S mini-module");

	const Int_t n = 24;
	const Int_t m = 43;
 
	Double_t x1[n] = {	1.32089,
						1.48578,
						1.57853,
						1.71422,
						1.79152,
						1.87568,
						1.93064,
						1.98904,
						2.13847,
						2.42703,
						2.51978,
						2.80491,
						2.97495,
						3.1656,
						4.26316,
						4.66851,
						5.15975,
						6.32257,
						7.25867,
						9.79902,
						10.888,
						12.2483,
						13.0642,
						15.0738};

    Double_t y1[n] = {	0.00162849,
						0.00298557,
						0.00298557,
						0.02605590,
						0.24183200,
						0.49832000,
						0.66659800,
						0.81587700,
						0.97872600,
						0.98279800,
						0.98144100,
						0.98415500,
						0.98279800,
						0.98415500,
						0.98415500,
						0.98279800,
						0.98415500,
						0.98415500,
						0.98415500,
						0.98279800,
						0.98415500,
						0.98415500,
						0.98415500,
						0.98551200};

    Double_t ex1[n] = {0.};
    Double_t ey1[n] = {0.};

    Double_t x2[m] = {	1.00,
						1.10,
						1.20,
						1.30,
						1.40,
						1.50,
						1.60,
						1.63,
						1.66,
						1.69,
						1.70,
						1.72,
						1.75,
						1.78,
						1.80,
						1.81,
						1.84,
						1.87,
						1.90,
						1.93,
						1.96,
						1.99,
						2.00,
						2.02,
						2.05,
						2.08,
						2.10,
						2.11,
						2.14,
						2.17,
						2.20,
						2.30,
						2.40,
						2.50,
						2.60,
						2.70,
						2.80,
						2.90,
						3.00,
						3.10,
						3.20,
						3.30,
						3.40};

    Double_t y2[m] = {	0.003410,
						0.003500,
						0.003573,
						0.003530,
						0.003498,
						0.003422,
						0.003941,
						0.004520,
						0.005858,
						0.007763,
						0.008178,
						0.009758,
						0.013008,
						0.042132,
						0.099101,
						0.131649,
						0.226638,
						0.320932,
						0.410526,
						0.497192,
						0.581813,
						0.662773,
						0.689642,
						0.741712,
						0.817616,
						0.890964,
						0.937125,
						0.955632,
						0.980813,
						0.983966,
						0.985541,
						0.988737,
						0.990732,
						0.991758,
						0.992378,
						0.992363,
						0.992336,
						0.992471,
						0.992546,
						0.992666,
						0.992572,
						0.992495,
						0.992585};
    Double_t ex2[m] = {0.};
    Double_t ey2[m] = {	0.0001751824,
						0.0001058767,
						0.0001737239,
						0.0002390676,
						0.0002249724,
						0.0002168391,
						0.0002315863,
						0.0002106161,
						0.0001951977,
						0.0002878101,
						0.0003769537,
						0.0003746398,
						0.0003657131,
						0.0005516748,
						0.0010768980,
						0.0008963404,
						0.0012533680,
						0.0021316580,
						0.0013478380,
						0.0013506140,
						0.0014890100,
						0.0013906180,
						0.0011904360,
						0.0016616650,
						0.0012097190,
						0.0009347524,
						0.0009852348,
						0.0007244555,
						0.0004788929,
						0.0005969346,
						0.0005259888,
						0.0002960620,
						0.0002843341,
						0.0002788878,
						0.0002643942,
						0.0003160370,
						0.0002290850,
						0.0002347978,
						0.0002987325,
						0.0002991991,
						0.0003590853,
						0.0002544258,
						0.0002663712};

    TGraphErrors *gr1 = new TGraphErrors(n,x1,y1,ex1,ey1);
    gr1->SetName("gr1");
    gr1->SetMarkerColor(2);
    gr1->SetMarkerStyle(20);
    gr1->SetMarkerSize(1.2);

    TGraphErrors *gr2 = new TGraphErrors(m,x2,y2,ex2,ey2);
    gr2->SetName("gr2");
    gr2->SetMarkerColor(4);
    gr2->SetMarkerStyle(20);
    gr2->SetMarkerSize(1.2);

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(gr1);
    mg->Add(gr2);
    mg->SetTitle("Stub efficiency for 2S mini-module");
    mg->Draw("AP");

    TAxis *xaxis = mg->GetXaxis();
    TAxis *yaxis = mg->GetYaxis();
    xaxis->SetTitle("Emulated pT [GeV]");
    xaxis->SetRangeUser(0, 3.5);
    yaxis->SetTitle("Stub efficiency");

    TF1* func1 = new TF1("func1", "([0]/(1+ TMath::Exp(-[1]*(x-[2]))))", x1[0], x1[9]);
    func1->SetParameters(1, 15, 2);
    func1->SetLineColor(kRed);
    func1->SetLineWidth(4);
    func1->SetLineStyle(9);
    //func1->SetParLimits(0,0.01,1);
    //func1->SetParLimits(1,0.01,100);
    //func1->SetParLimits(2,1,3);

    gr1->Fit(func1);

    TF1* func2 = new TF1("func2", "([0]/(1+ TMath::Exp(-[1]*(x-[2]))))", x1[0], x1[9]);
    func2->SetParameters(1, 15, 2);
    func2->SetLineColor(kBlue);
    func2->SetLineWidth(4);
    func2->SetLineStyle(9);
    //func2->SetParLimits(0,0.01,1);
    //func2->SetParLimits(1,0.01,100);
    //func2->SetParLimits(2,1,3);

    gr2->Fit(func2);

    auto legend = new TLegend(0.1,0.7,0.35,0.9);
    legend->SetHeader("Legend","C");
    legend->AddEntry("gr1","Adam et all. 2020","p");
    legend->AddEntry("gr2","Geant4","ep");
    legend->AddEntry(func1,"turn-on p_{T} = 1.88 GeV","l");
    legend->AddEntry(func2,"turn-on p_{T} = 1.94 GeV","l");
    legend->Draw();

    gPad->Modified();

    /*axis->SetAxisColor(Color_t color = 1);
    axis->SetLabelColor(Color_t color = 1);
    axis->SetLabelFont(Style_t font = 62);
    axis->SetLabelOffset(Float_t offset = 0.005);
    axis->SetLabelSize(Float_t size = 0.04);
    axis->SetNdivisions(Int_t n = 510, Bool_t optim = kTRUE);
    axis->SetNoExponent(Bool_t noExponent = kTRUE);
    axis->SetTickLength(Float_t length = 0.03);
    axis->SetTitleOffset(Float_t offset = 1);
    axis->SetTitleSize(Float_t size = 0.02)*/

    /*TF1 *f1 = new TF1("f1","myFit1(x)", 1.4, 2.6);
    f1->SetLineColor(kRed);
    f1->SetLineWidth(1);
    f1->Draw("same");

    TF1 *f2 = new TF1("f2","myFit2(x)", 1.4, 2.6);
    f2->SetLineColor(kBlue);
    f2->SetLineWidth(1);
    f2->Draw("same");*/
}