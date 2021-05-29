#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <random>


void figure14A() {

    const int NBR_BINS = 14;

    // array creation
    Double_t x1[NBR_BINS] = {	0.69750,
                                1.46438,
                                2.23125,
                                2.95800,
                                5.21850,
                                5.98538,
                                7.47874,
                                9.69881,
                                11.23268,
                                16.48005,
                                19.46704,
                                20.97015,
                                22.45403,
                                23.94750}; // in keV !!!
    Double_t y1[NBR_BINS] = {	5.3042E+01,
                                1.4471E+01,
                                3.7583E+00,
                                1.5217E+00,
                                1.0055E+00,
                                1.0043E+00,
                                1.0043E+00,
                                1.0051E+00,
                                1.0055E+00,
                                1.0018E+00,
                                9.1113E-01,
                                7.8413E-01,
                                6.3020E-01,
                                4.8454E-01};

    Double_t x2[NBR_BINS] = {   0.70200,
                                1.46175,
                                2.22150,
                                2.95538,
                                5.20875,
                                5.96850,
                                7.46201,
                                9.71539,
                                11.20890,
                                16.44934,
                                19.46258,
                                20.95609,
                                22.44960,
                                23.96933}; // in keV !!!
    Double_t y2[NBR_BINS] = {	5.2740E+01,
                                1.4549E+01,
                                3.8212E+00,
                                1.5309E+00,
                                1.0165E+00,
                                1.0153E+00,
                                1.0159E+00,
                                1.0168E+00,
                                1.0171E+00,
                                1.0116E+00,
                                9.2189E-01,
                                7.9677E-01,
                                6.4258E-01,
                                4.9255E-01};

    // adam D0
    Double_t x3[NBR_BINS] = {   0.69750,
                                1.46438,
                                2.23125,
                                2.95800,
                                5.21850,
                                5.98538,
                                7.47874,
                                9.69881,
                                11.23268,
                                16.48005,
                                19.46704,
                                20.97015,
                                22.45403,
                                23.94750}; // in keV !!!
    Double_t y3[NBR_BINS] = {   34.313600,
                                21.360400,
                                3.745310,
                                1.221800,
                                1.025700,
                                1.024490,
                                1.023320,
                                1.021920,
                                1.014540,
                                0.912504,
                                0.692076,
                                0.559103,
                                0.443103,
                                0.345991};

    // adam D1
    Double_t x4[NBR_BINS] = {   0.73800,
                                1.50488,
                                2.23125,
                                2.99813,
                                5.25863,
                                5.98538,
                                7.47874,
                                9.73916,
                                11.23268,
                                16.49768,
                                19.50739,
                                21.00090,
                                22.49438,
                                23.98785}; // in keV !!!
    Double_t y4[NBR_BINS] = {   34.737000,
                                18.986900,
                                2.858890,
                                1.150090,
                                1.041510,
                                1.041510,
                                1.041510,
                                1.041510,
                                1.032160,
                                0.934708,
                                0.706798,
                                0.574432,
                                0.450319,
                                0.356219};

    Double_t y5[NBR_BINS] = {   0.};

    Double_t y6[NBR_BINS] = {   0.};

    for (int i = 0; i < NBR_BINS; ++i)
    {
    	x1[i] *= 1000.;
    	x2[i] *= 1000.;
    	y5[i] = (y1[i] - y3[i]) / y3[i];
    	y6[i] = (y2[i] - y4[i]) / y4[i];
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

    TGraph *gr2 = new TGraph(NBR_BINS,x2,y6);
    gr2->SetName("gr2");
    gr2->SetMarkerColor(kRed+2);
    gr2->SetMarkerStyle(25);
    gr2->SetMarkerSize(1.);

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(gr1);
    mg->Add(gr2);
    mg->SetTitle("");
    mg->Draw("AP");

    TAxis *xaxis = mg->GetXaxis();
    TAxis *yaxis = mg->GetYaxis();
    xaxis->SetMaxDigits(3);
    xaxis->SetLabelFont(42);
    xaxis->SetLabelSize(0.04);
    xaxis->SetTitle("Seuil [e]");
    xaxis->SetTitleFont(22);
    xaxis->SetTitleSize(0.05);
    xaxis->SetTitleOffset(0.95);
    xaxis->SetRangeUser(0., 25000);
    yaxis->SetLabelFont(42);
    yaxis->SetLabelSize(0.04);
    yaxis->SetTitle("Erreur relative");
    yaxis->SetTitleFont(22);
    yaxis->SetTitleSize(0.05);
    yaxis->SetTitleOffset(0.5);
    yaxis->SetRangeUser(-0.5, 0.5);

    c1->RedrawAxis();

    TF1* f1 = new TF1("f1", "x", 0, 25000);
    TGaxis* A1 = new TGaxis(0., 0.5, 25000, 0.5, "f1", 510, "-");
    A1->SetLabelFont(42);
    A1->SetLabelSize(0.04);
    A1->SetTitle("Seuil [keV]");
    A1->SetTitleFont(22);
    A1->SetTitleSize(0.05);
    A1->SetTitleOffset(0.95);
    A1->ChangeLabel(1, -1, -1, -1, -1, -1, "0.0");
    A1->ChangeLabel(2, -1, -1, -1, -1, -1, "18.1");
    A1->ChangeLabel(3, -1, -1, -1, -1, -1, "36.2");
    A1->ChangeLabel(4, -1, -1, -1, -1, -1, "54.3");
    A1->ChangeLabel(5, -1, -1, -1, -1, -1, "72.4");
    A1->ChangeLabel(6, -1, -1, -1, -1, -1, "90.5");
    //A1->ChangeLabel(7, -1, -1, -1, -1, -1, "1.4");
    A1->Draw("SAME");

    TLine *line = new TLine(xaxis->GetXmin(),0.,xaxis->GetXmax(),0.);
    line->SetLineWidth(1.5);
    line->SetLineColor(12);
    line->Draw("SAME");

    auto legend = new TLegend(0.65,0.35,0.9,0.1);
    legend->AddEntry("gr1","Geant4 - top sensor","p");
    legend->AddEntry("gr2","Geant4 - bottom sensor","p");
    legend->Draw();

    gPad->Modified();
    //*********************** 

    c1->SaveAs("figure14A.pdf");
    // Close file when finished
    //f.Close();
}