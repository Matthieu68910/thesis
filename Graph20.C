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
	auto c1 = new TCanvas("c1","c1",600,500);
    gStyle->SetOptStat(0);

	const Int_t n = 24;
	const Int_t m = 24;
 
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

    Double_t y2[m] = {	0.0016030,
						0.0015846,
						0.0019906,
						0.0075969,
						0.0713456,
						0.3366737,
						0.4985349,
						0.6597470,
						0.9804535,
						0.9913831,
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
						0.00003643625,
						0.00013297070,
						0.00026055080,
						0.00048082290,
						0.00036160800,
						0.00044961120,
						0.00011370840,
						0.00006699826,
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
    mg->Draw("AP");

    /*TF1 *f1 = new TF1("f1","myFit1(x)", 1.4, 2.6);
    f1->SetLineColor(kRed);
    f1->SetLineWidth(1);
    f1->Draw("same");

    TF1 *f2 = new TF1("f2","myFit2(x)", 1.4, 2.6);
    f2->SetLineColor(kBlue);
    f2->SetLineWidth(1);
    f2->Draw("same");*/

    TF1* func1 = new TF1("func", "([0]/(1+ TMath::Exp(-[1]*(x-[2]))))", x1[0], x1[9]);
    func1->SetParameters(1, 15, 2);
    //func1->SetParLimits(0,0.01,1);
    //func1->SetParLimits(1,0.01,100);
    //func1->SetParLimits(2,1,3);

    gr1->Fit(func1);

    TF1* func2 = new TF1("func", "([0]/(1+ TMath::Exp(-[1]*(x-[2]))))", x1[0], x1[9]);
    func2->SetParameters(1, 15, 2);
    //func2->SetParLimits(0,0.01,1);
    //func2->SetParLimits(1,0.01,100);
    //func2->SetParLimits(2,1,3);

    gr2->Fit(func2);
}