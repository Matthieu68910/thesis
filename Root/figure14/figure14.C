#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <random>

std::default_random_engine generator;
std::uniform_real_distribution<> distribution1(0.0, 1.0);

bool CBC2(
    const vector<double> &strip_A, 
    const vector<double> &strip_B, 
    vector<double> &res_A, // [0:4]-> 1-5 strip wide clusters,[5]-> number of clusters (tot), [6]-> mean cluster width (tot), [7:8]-> accepted, [9] -> nbr hits
    vector<double> &res_B,
    const int MAX_CLUSTER_WIDTH = 3,
    const int CLUSTER_WINDOW = 5,
    double THRESHOLD = 5.1975
    //double kill_value = 0.04331
    ){

    const int NBR_STRIP = strip_A.size();

    double noise = 0.;
    // Clusters in sensor A
    std::vector<double> clus_pos_A;
    std::vector<double> clus_size_A;
    bool inside = false;
    int size = 0;
    int nbr_hits_A = 0;
    // Loop on sensor A strips
    for (int i = 0; i < NBR_STRIP; ++i)
    {
        // noise parameters determination
        if (distribution1(generator) >= 0.5)
        {
            std::normal_distribution<double> dist(1.36, 0.06);
            noise = abs(dist(generator)) * 0.375;
        } else
        {
            std::normal_distribution<double> dist(2.38, 0.6);
            noise = abs(dist(generator)) * 0.375;
        }
        // noise value deternmination
        std::normal_distribution<double> dist1(0., noise);
        // noise creation
        double strip_energy = (strip_A[i] / 0.00362) + abs(dist1(generator));
        //if(distribution1(generator) < kill_value){strip_energy = 0.;}
        if (strip_energy < THRESHOLD && !inside)        
        {} else if (strip_energy < THRESHOLD && inside)
        {
            clus_size_A.push_back(size);
            if(size <= 5) res_A.at(size - 1) += 1;
            clus_pos_A.push_back(floor((i - 1) - (size / 2) + 0.5));
            size = 0;
            inside = false;
        } else if (strip_energy >= THRESHOLD && !inside)
        {
            nbr_hits_A += 1;
            size = 1;
            inside = true;
            if (i == (NBR_STRIP - 1))
            {
                clus_size_A.push_back(1);
                res_A.at(0) += 1;
                clus_pos_A.push_back(i);
                size = 0;
                inside = false;
            }
        } else if (strip_energy >= THRESHOLD && inside)
        {
            nbr_hits_A += 1;
            size += 1;
            if (i == (NBR_STRIP - 1))
            {
                clus_size_A.push_back(size);
                if(size <= 5) res_A.at(size - 1) += 1;
                clus_pos_A.push_back(floor(i - (size / 2) + 0.5));
                size = 0;
                inside = false;
            }
        }
    }

    // Clusters in sensor B
    std::vector<double> clus_pos_B;
    std::vector<double> clus_size_B;
    inside = false;
    size = 0;
    int nbr_hits_B = 0;
    // Loop on sensor B strips
    for (int i = 0; i < NBR_STRIP; ++i)
    {
        // noise parameters determination
        if (distribution1(generator) >= 0.5)
        {
            std::normal_distribution<double> dist(1.36, 0.06);
            noise = abs(dist(generator)) * 0.375;
        } else
        {
            std::normal_distribution<double> dist(2.38, 0.6);
            noise = abs(dist(generator)) * 0.375;
        }
        // noise value deternmination
        std::normal_distribution<double> dist1(0., noise);
        // noise creation
        double strip_energy = (strip_B[i] / 0.00362) + abs(dist1(generator));
        //if(distribution1(generator) < kill_value){strip_energy = 0.;}
        if (strip_energy < THRESHOLD && !inside)        
        {} else if (strip_energy < THRESHOLD && inside)
        {
            clus_size_B.push_back(size);
            if(size <= 5) res_B.at(size - 1) += 1;
            clus_pos_B.push_back(floor((i - 1) - (size / 2) + 0.5));
            size = 0;
            inside = false;
        } else if (strip_energy >= THRESHOLD && !inside)
        {
            nbr_hits_B += 1;
            size = 1;
            inside = true;
            if (i == (NBR_STRIP - 1))
            {
                clus_size_B.push_back(1);
                res_B.at(0) += 1;
                clus_pos_B.push_back(i);
                size = 0;
                inside = false;
            }
        } else if (strip_energy >= THRESHOLD && inside)
        {
            nbr_hits_B += 1;
            size += 1;
            if (i == (NBR_STRIP - 1))
            {
                clus_size_B.push_back(size);
                if(size <= 5) res_B.at(size - 1) += 1;
                clus_pos_B.push_back(floor(i - (size / 2) + 0.5));
                size = 0;
                inside = false;
            }
        }
    }

    // fill results (tot)
    if(clus_pos_A.size() != 0){res_A.at(5) = clus_pos_A.size();}else{res_A.at(5) = 0;}
    if(clus_pos_B.size() != 0){res_B.at(5) = clus_pos_B.size();}else{res_B.at(5) = 0;}
    if(clus_size_A.size() != 0){res_A.at(6) = std::accumulate(clus_size_A.begin(), clus_size_A.end(), 0.0) / clus_size_A.size();}else{res_A.at(6) = 0;}
    if(clus_size_B.size() != 0){res_B.at(6) = std::accumulate(clus_size_B.begin(), clus_size_B.end(), 0.0) / clus_size_B.size();}else{res_B.at(6) = 0;}

    // Clean clusters in both sensors depending on MAX_CLUSTER_WIDTH
    // for A
    for (int i = clus_size_A.size() - 1; i >= 0; --i)
    {
        if (clus_size_A.at(i) > MAX_CLUSTER_WIDTH)
        {   
            clus_size_A.erase(clus_size_A.begin() + i);
            clus_pos_A.erase(clus_pos_A.begin() + i);
        }
    }
    // for B
    for (int i = clus_size_B.size() - 1; i >= 0; --i)
    {
        if (clus_size_B.at(i) > MAX_CLUSTER_WIDTH)
        {   
            clus_size_B.erase(clus_size_B.begin() + i);
            clus_pos_B.erase(clus_pos_B.begin() + i);
        }
    }

    // fill results (accepted)
    res_A.at(7) = clus_pos_A.size();
    res_B.at(7) = clus_pos_B.size();
    res_A.at(8) = std::accumulate(clus_size_A.begin(), clus_size_A.end(), 0.0) / clus_size_A.size();
    res_B.at(8) = std::accumulate(clus_size_B.begin(), clus_size_B.end(), 0.0) / clus_size_B.size();
    res_A.at(9) = nbr_hits_A;
    res_B.at(9) = nbr_hits_B;

    // return false if no possible stub
    if (clus_pos_A.size() == 0 || clus_pos_B.size() == 0)
    {
        return false;
    }

    // stub finding logic
    for (int i = 0; i < clus_pos_A.size(); ++i)
    {
        for (int j = 0; j < clus_pos_B.size(); ++j)
        {
            if (abs(clus_pos_A.at(i) - clus_pos_B.at(j)) <= CLUSTER_WINDOW)
            {
                return true;
            }
        }
    }
    return false;
}

void figure14() {
	// data for Adam2020

	// open file
    TFile *f = TFile::Open("/media/matthieu/ssd1/Geant4/Data/Data_figure13-16/data100k_pions.root", "read");

    //************* variable ***************//
    const int NBR_STRIP = 254;

    const int MAX_CLUSTER_WIDTH = 3;
    const int CLUSTER_WINDOW = 5;

    const int NBR_BINS = 14;
    const int NBR_BINS_IRR = 16;

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
    Double_t y1[NBR_BINS] = {0.};
    Double_t ex1[NBR_BINS] = {0.};
    Double_t ey1[NBR_BINS] = {0.};

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
    Double_t y2[NBR_BINS] = {0.};
    Double_t ex2[NBR_BINS] = {0.};
    Double_t ey2[NBR_BINS] = {0.};

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
    Double_t ex3[NBR_BINS] = {0.};
    Double_t ey3[NBR_BINS] = {0.};

    // adam D1
    Double_t x4[NBR_BINS] = {   0.70200,
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
    Double_t ex4[NBR_BINS] = {0.};
    Double_t ey4[NBR_BINS] = {0.};

    // adam D0_IRR
    Double_t x5[NBR_BINS_IRR] = {   0.70275,
                                    1.47000,
                                    2.20538,
                                    2.97263,
                                    4.47488,
                                    5.21025,
                                    5.97713,
                                    8.95009,
                                    9.71726,
                                    11.21966,
                                    14.22450,
                                    19.46693,
                                    19.47071,
                                    20.96933,
                                    22.47173,
                                    23.97413}; // in keV !!!
    Double_t y5[NBR_BINS_IRR] = {   25.7860000,
									5.5754500,
									1.3866300,
									0.8155560,
									0.7533410,
									0.7484100,
									0.7257880,
									0.6025730,
									0.5778540,
									0.5144990,
									0.3679600,
									0.1341110,
									0.1313330,
									0.1017400,
									0.0804747,
									0.0628790};
    Double_t ex5[NBR_BINS_IRR] = {0.};
    Double_t ey5[NBR_BINS_IRR] = {0.};

    // adam D1_IRR
    Double_t x6[NBR_BINS_IRR] = {   0.69750,
                                    1.44750,
                                    2.20538,
                                    2.94750,
                                    4.44750,
                                    5.19750,
                                    5.94750,
                                    8.94750,
                                    9.69750,
                                    11.19750,
                                    14.19750,
                                    19.44750,
                                    19.46693,
                                    20.97015,
                                    22.44750,
                                    23.94750}; // in keV !!!
    Double_t y6[NBR_BINS_IRR] = {   25.3367000,
									5.0908500,
									1.2586000,
									0.8018850,
									0.7596600,
									0.7528420,
									0.7327520,
									0.6063620,
									0.5848840,
									0.5155280,
									0.3760190,
									0.1394780,
									0.1369850,
									0.1064240,
									0.0849478,
									0.0671969};
    Double_t ex6[NBR_BINS_IRR] = {0.};
    Double_t ey6[NBR_BINS_IRR] = {0.};


    //************ Get the TTree ******************************************
    auto data = f->Get<TTree>("data");

    // Get the number of entries in TTree
    const int ENTRIES = data->GetEntries() / 10;
    //cout << std::scientific << "Number of entries: " << ENTRIES << endl;

    //**************** Set BranchAddress for datas recovery ***************
    // for strips
    std::vector<double> strip_A(NBR_STRIP, 0);
    std::vector<double> strip_B(NBR_STRIP, 0);

    for (int i = 0; i < NBR_STRIP; ++i)
    {
        string strip_name = "s" + std::to_string(i);
        char const *pchar = strip_name.c_str();
        data->SetBranchAddress(pchar, &strip_B.at(i));
    }
    for (int i = NBR_STRIP; i < (2*NBR_STRIP); ++i)
    {
        string strip_name = "s" + std::to_string(i);
        char const *pchar = strip_name.c_str();
        data->SetBranchAddress(pchar, &strip_A.at(i-NBR_STRIP));
    }
    // for other datas (pT)

    //****************** Create Histo ************************************//
    auto c1 = new TCanvas("c1","c1", 1000, 600);
	c1->SetTitle("Figure 1: threshold variation");
	gStyle->SetOptStat(0);
	gPad->SetGridx(1);
    gPad->SetGridy(1);
    gPad->SetLogy();

    gPad->SetTitle("Figure 1: threshold variation");

    for (int j = 0; j < NBR_BINS; ++j)
    {
	    // Create vectors for stat
	    std::vector<double> nbr_hitsA(100, 0);
	    std::vector<double> nbr_hitsB(100, 0);

	    //****************** Main loop over all entries **********************//
	    int count_loop = 0;
	    int loop_number = 0;
	    for (int k = 0; k < ENTRIES; k++)
	    {
	        // fill variables with datas from entry i
	        data->GetEntry(k);

	        std::vector<double> res_A1(10, 0);
	        std::vector<double> res_B1(10, 0);
            std::vector<double> res_A2(10, 0);
            std::vector<double> res_B2(10, 0);

	        bool stub1 = CBC2(strip_A, strip_B, res_A1, res_B1, MAX_CLUSTER_WIDTH, CLUSTER_WINDOW, x1[j]);
            bool stub2 = CBC2(strip_A, strip_B, res_A2, res_B2, MAX_CLUSTER_WIDTH, CLUSTER_WINDOW, x2[j]);

	        nbr_hitsA.at(loop_number) += (double) res_A1.at(5) / (ENTRIES / 100);

	        nbr_hitsB.at(loop_number) += (double) res_B2.at(5) / (ENTRIES / 100);

	        count_loop += 1;
	        if (count_loop == ENTRIES / 100)
	        {
	        	count_loop = 0;
	        	loop_number += 1;
	        }
	    }
	    double variance1 = 0., deviation1 = 0., average1 = 0.;
		average1 = std::accumulate(nbr_hitsA.begin(), nbr_hitsA.end(), 0.0) / 100;
		for (int i = 0; i < 100; ++i){variance1 += pow((nbr_hitsA.at(i) - average1), 2);}
		variance1 /= 99;
		deviation1 = sqrt(variance1);
		

		y1[j] = average1;
		ey1[j] = deviation1;

		double variance2 = 0., deviation2 = 0., average2 = 0.;
		average2 = std::accumulate(nbr_hitsB.begin(), nbr_hitsB.end(), 0.0) / 100;
		for (int i = 0; i < 100; ++i){variance2 += pow((nbr_hitsB.at(i) - average2), 2);}
		variance2 /= 99;
		deviation2 = sqrt(variance2);
		cout << std::scientific << x1[j] << "\t" << average1 << "\t" << deviation1 << "\t" << average2 << "\t" << deviation2 << endl;

		y2[j] = average2;
		ey2[j] = deviation2;
	}
    for (int i = 0; i < NBR_BINS; ++i)
    {
        x1[i] *= 1000.;
        x2[i] *= 1000.;
        x3[i] *= 1000.;
        x4[i] *= 1000.;
    }
    for (int i = 0; i < NBR_BINS_IRR; ++i)
    {
        x5[i] *= 1000.;
        x6[i] *= 1000.;
    }

    TGraphErrors *gr1 = new TGraphErrors(NBR_BINS,x1,y1,ex1,ey1);
    gr1->SetName("gr1");
    gr1->SetMarkerColor(kRed+2);
    gr1->SetMarkerStyle(20);
    gr1->SetMarkerSize(1.0);

    TGraphErrors *gr2 = new TGraphErrors(NBR_BINS,x2,y2,ex2,ey2);
    gr2->SetName("gr2");
    gr2->SetMarkerColor(kRed+2);
    gr2->SetMarkerStyle(25);
    gr2->SetMarkerSize(1.0);

    TGraphErrors *gr3 = new TGraphErrors(NBR_BINS,x3,y3,ex3,ey3);
    gr3->SetName("gr3");
    gr3->SetMarkerColor(12);
    gr3->SetMarkerStyle(20);
    gr3->SetMarkerSize(1.0);

    TGraphErrors *gr4 = new TGraphErrors(NBR_BINS,x4,y4,ex4,ey4);
    gr4->SetName("gr4");
    gr4->SetMarkerColor(12);
    gr4->SetMarkerStyle(25);
    gr4->SetMarkerSize(1.0);

    TGraphErrors *gr5 = new TGraphErrors(NBR_BINS_IRR,x5,y5,ex5,ey5);
    gr5->SetName("gr5");
    gr5->SetMarkerColor(12);
    gr5->SetMarkerStyle(22);
    gr5->SetMarkerSize(1.0);

    TGraphErrors *gr6 = new TGraphErrors(NBR_BINS_IRR,x6,y6,ex6,ey6);
    gr6->SetName("gr6");
    gr6->SetMarkerColor(12);
    gr6->SetMarkerStyle(26);
    gr6->SetMarkerSize(1.0);

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(gr3);
    mg->Add(gr4);
    //mg->Add(gr5);
    //mg->Add(gr6);
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
    yaxis->SetTitle("Nombre moyen de cluster/evenement");
    yaxis->SetTitleFont(22);
    yaxis->SetTitleSize(0.05);
    yaxis->SetTitleOffset(0.9);

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

    auto legend = new TLegend(0.55,0.9,0.9,0.7);
    legend->AddEntry("gr1","Geant4 - top sensor","ep");
    legend->AddEntry("gr2","Geant4 - bottom sensor","ep");
    legend->AddEntry("gr3","Adam et al. - top sensor","p");
    legend->AddEntry("gr4","Adam et al. - bottom sensor","p");
    //legend->AddEntry("gr5","Adam et al. - top sensor irr.","p");
    //legend->AddEntry("gr6","Adam et al. - bottom sensor irr.","p");
    legend->Draw();

    gPad->Modified();
    //*********************** 

    c1->SaveAs("figure14_386.pdf");
    // Close file when finished
    //f.Close();
}