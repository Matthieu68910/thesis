#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <random>


void figure20A() {

    const int NBR_BINS = 25;
    const int NBR_BINS1 = 24;

    // array creation
    Double_t x1[NBR_BINS] = {	1.32089E+00,
                                1.48578E+00,
                                1.57853E+00,
                                1.71422E+00,
                                1.79152E+00,
                                1.87568E+00,
                                1.93064E+00,
                                1.98904E+00,
                                2.13847E+00,
                                2.42703E+00,
                                2.51978E+00,
                                2.80491E+00,
                                2.97495E+00,
                                3.16560E+00,
                                4.26316E+00,
                                4.66851E+00,
                                5.15975E+00,
                                6.32257E+00,
                                6.32257E+00,
                                7.25867E+00,
                                9.79902E+00,
                                1.08880E+01,
                                1.22483E+01,
                                1.30642E+01,
                                1.50738E+01
                                }; // in keV !!!
    Double_t y1[NBR_BINS] = {	4.1000000E-03,
                                3.7300000E-03,
                                5.0100000E-03,
                                3.0780000E-02,
                                2.6073000E-01,
                                5.2267000E-01,
                                6.7643000E-01,
                                8.3019000E-01,
                                9.8634000E-01,
                                9.9261000E-01,
                                9.9244000E-01,
                                9.9267000E-01,
                                9.9288000E-01,
                                9.9277000E-01,
                                9.9281000E-01,
                                9.9307000E-01,
                                9.9318000E-01,
                                9.9367000E-01,
                                9.9269000E-01,
                                9.9314000E-01,
                                9.9323000E-01,
                                9.9373000E-01,
                                9.9379000E-01,
                                9.9316000E-01,
                                9.9362000E-01
                                };

    Double_t x2[NBR_BINS1] = {   1.32089, 
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
                                15.0738}; // in keV !!!}; // in keV !!!
    Double_t y2[NBR_BINS1] = {	0.00162849,
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


    Double_t y5[NBR_BINS] = {   0.};

    for (int i = 0; i < NBR_BINS; ++i)
    {
    	y5[i] = (y1[i] - y2[i]) / y2[i];
    }


    //****************** Create Histo ************************************//
    auto c1 = new TCanvas("c1","c1",1000,400);
	c1->SetTitle("Figure 1: threshold variation");
	gStyle->SetOptStat(0);
	gPad->SetGridx(1);
    gPad->SetGridy(1);

    gPad->SetTitle("Figure 1: threshold variation");

    TGraph *gr1 = new TGraph(NBR_BINS,x1,y5);
    gr1->SetName("gr1");
    gr1->SetMarkerColor(kRed+2);
    gr1->SetMarkerStyle(20);
    gr1->SetMarkerSize(1.);

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(gr1);
    mg->SetTitle("");
    mg->Draw("AP");

    TAxis *xaxis = mg->GetXaxis();
    TAxis *yaxis = mg->GetYaxis();
    xaxis->SetLabelFont(42);
    xaxis->SetLabelSize(0.04);
    xaxis->SetTitle("pT [GeV]");
    xaxis->SetTitleFont(22);
    xaxis->SetTitleSize(0.05);
    xaxis->SetTitleOffset(0.95);
    xaxis->SetRangeUser(1.0, 3.5);
    yaxis->SetLabelFont(42);
    yaxis->SetLabelSize(0.04);
    yaxis->SetTitle("#varepsilon_{rel} = #frac{x_{sim} #minus x_{exp}}{x_{exp}}");
    yaxis->SetTitleFont(22);
    yaxis->SetTitleSize(0.05);
    yaxis->SetTitleOffset(0.9);

    TF1* f1 = new TF1("f1", "x", 1.0, 3.5); // 3.5
    TGaxis* A1 = new TGaxis(1.0, yaxis->GetXmax(), 3.5, yaxis->GetXmax(), "f1", 510, "-");
    A1->SetLabelFont(42);
    A1->SetLabelSize(0.04);
    A1->SetTitle("Angle d'incidence [deg]");
    A1->SetTitleFont(22);
    A1->SetTitleSize(0.05);
    A1->SetTitleOffset(0.95);
    A1->ChangeLabel(1, -1, -1, -1, -1, -1, "19.6");
    A1->ChangeLabel(2, -1, -1, -1, -1, -1, "13.1");
    A1->ChangeLabel(3, -1, -1, -1, -1, -1, "9.8");
    A1->ChangeLabel(4, -1, -1, -1, -1, -1, "7.8");
    A1->ChangeLabel(5, -1, -1, -1, -1, -1, "6.5");
    A1->ChangeLabel(6, -1, -1, -1, -1, -1, "5.6");
    A1->ChangeLabel(7, -1, -1, -1, -1, -1, "");
    A1->Draw("SAME");

    TLine *line = new TLine(1.0,0.,3.5,0.);
    line->SetLineWidth(1.5);
    line->SetLineColor(12);
    line->Draw("SAME");

    auto legend = new TLegend(0.65,0.85,0.85,0.75);
    legend->AddEntry("gr1","Geant4","p");
    legend->Draw();

    gPad->Modified();
    //*********************** 

    c1->SaveAs("figure20A.pdf");
    // Close file when finished
    //f.Close();
}