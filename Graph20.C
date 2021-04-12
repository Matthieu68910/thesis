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
	const Int_t m = 33;
 
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

    Double_t x2[m] = {	1.32089,
						1.48578,
						1.53000,
						1.57853,
						1.65000,
						1.71422,
						1.75000,
						1.79152,
						1.83000,
						1.87568,
						1.90000,
						1.93064,
						1.95500,
						1.98904,
						2.05000,
						2.13847,
						2.25000,
						2.42703,
						2.47000,
						2.51978,
						2.80491,
						2.97495,
						3.16560,
						4.26316,
						4.66851,
						5.15975,
						6.32257,
						7.25867,
						9.79902,
						10.88800,
						12.24830,
						13.06420,
						15.07380};

    Double_t y2[m] = {	0.0016030,
						0.0015846,
						0.0016105,
						0.0019906,
						0.0036492,
						0.0075969,
						0.0114934,
						0.0713456,
						0.1942445,
						0.3366737,
						0.4094951,
						0.4985349,
						0.5670335,
						0.6597470,
						0.8178548,
						0.9804535,
						0.9875308,
						0.9913831,
						0.9919240,
						0.9922130,
						0.9928290,
						0.9926380,
						0.9928450,
						0.9929330,
						0.9930480,
						0.9930800,
						0.9930800,
						0.9931430,
						0.9932380,
						0.9931540,
						0.9932470,
						0.9932640,
						0.9932580};
    Double_t ex2[m] = {0.};
    Double_t ey2[m] = {	0.00014989260,
						0.00003747355,
						0.00003985878,
						0.00003643625,
						0.00008181388,
						0.00013297070,
						0.00010621380,
						0.00026055080,
						0.00032839350,
						0.00048082290,
						0.00040088330,
						0.00036160800,
						0.00031512260,
						0.00044961120,
						0.00028656420,
						0.00011370840,
						0.00014011330,
						0.00006699826,
						0.00007939913,
						0.00008266398,
						0.00036722530,
						0.00024160570,
						0.00013492800,
						0.00022005300,
						0.00038406600,
						0.00020741800,
						0.00020768570,
						0.00030114780,
						0.00019106720,
						0.00030111640,
						0.00020559130,
						0.00022789860,
						0.00021410280};

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
    legend->AddEntry(func2,"turn-on p_{T} = 1.93 GeV","l");
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