#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <random>


void figure13A() {

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
    Double_t y1[NBR_BINS] = {	9.7275E+01,
								1.7652E+01,
								2.3284E+00,
								1.1142E+00,
								1.0581E+00,
								1.0563E+00,
								1.0522E+00,
								1.0460E+00,
								1.0415E+00,
								1.0264E+00,
								9.3548E-01,
								8.0852E-01,
								6.5354E-01,
								5.0348E-01};

    Double_t x2[NBR_BINS] = {   0.73800,
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
    Double_t y2[NBR_BINS] = {	9.0413E+01,
								1.5879E+01,
								2.3654E+00,
								1.1443E+00,
								1.0944E+00,
								1.0916E+00,
								1.0850E+00,
								1.0733E+00,
								1.0660E+00,
								1.0423E+00,
								9.4739E-01,
								8.2097E-01,
								6.6396E-01,
								5.1355E-01};

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
    Double_t y3[NBR_BINS] = {   38.080700,
                                20.937800,
                                4.368410,
                                1.542830,
                                1.157780,
                                1.144010,
                                1.116960,
                                1.090560,
                                1.077590,
                                0.944719,
                                0.717475,
                                0.581773,
                                0.466412,
                                0.362798};

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
    Double_t y4[NBR_BINS] = {   41.406900,
                                19.959500,
                                3.438860,
                                1.470740,
                                1.200090,
                                1.185810,
                                1.144010,
                                1.116960,
                                1.103680,
                                0.972574,
                                0.734848,
                                0.599618,
                                0.472025,
                                0.376055};

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
    gr1->SetMarkerStyle(50);
    gr1->SetMarkerSize(1.2);

    TGraph *gr2 = new TGraph(NBR_BINS,x2,y6);
    gr2->SetName("gr2");
    gr2->SetMarkerColor(kBlue+2);
    gr2->SetMarkerStyle(50);
    gr2->SetMarkerSize(1.2);

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(gr1);
    mg->Add(gr2);
    mg->SetTitle("");
    mg->Draw("AP");

    TAxis *xaxis = mg->GetXaxis();
    TAxis *yaxis = mg->GetYaxis();
    xaxis->SetLabelFont(42);
    xaxis->SetLabelSize(0.04);
    xaxis->SetTitle("Seuil [e-]");
    xaxis->SetTitleFont(22);
    xaxis->SetTitleSize(0.05);
    xaxis->SetTitleOffset(0.95);
    xaxis->SetRangeUser(0., 25000);
    yaxis->SetLabelFont(42);
    yaxis->SetLabelSize(0.04);
    yaxis->SetTitle("Erreur relative");
    yaxis->SetTitleFont(22);
    yaxis->SetTitleSize(0.05);
    yaxis->SetTitleOffset(0.9);
    yaxis->SetRangeUser(-0.6, 0.6);

    c1->RedrawAxis();

    TF1* f1 = new TF1("f1", "x", 0, 25000);
    TGaxis* A1 = new TGaxis(0., yaxis->GetXmax(), 25000, yaxis->GetXmax(), "f1", 510, "-");
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

    c1->RedrawAxis();

    auto legend = new TLegend(0.65,0.35,0.9,0.1);
    legend->AddEntry("gr1","Geant4 - top sensor","p");
    legend->AddEntry("gr2","Geant4 - bottom sensor","p");
    legend->Draw();

    gPad->Modified();
    //*********************** 

    c1->SaveAs("figure13A.pdf");
    // Close file when finished
    //f.Close();
}