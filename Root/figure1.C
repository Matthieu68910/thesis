#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <random>

std::default_random_engine generator;
std::normal_distribution<double> distribution(0.00362, 0.00220994475); // µ = 1000 e-, s = 800 e-

bool CBC2(
	const vector<double> &strip_A, 
	const vector<double> &strip_B, 
	vector<double> &res_A, // [0:4]-> 1-5 strip wide clusters,[5]-> number of clusters (tot), [6]-> mean cluster width (tot), [7:8]-> accepted
	vector<double> &res_B,
	const int MAX_CLUSTER_WIDTH = 3,
	const int CLUSTER_WINDOW = 5,
	double THRESHOLD = 0.0148
	){

	const int NBR_STRIP = strip_A.size();

	// Clusters in sensor A
	std::vector<double> clus_pos_A;
	std::vector<double> clus_size_A;
	bool inside = false;
	int size = 0;
	// Loop on sensor A strips
	for (int i = 0; i < NBR_STRIP; ++i)
    {
    	double strip_energy = abs(strip_A[i] + distribution(generator));
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
	// Loop on sensor B strips
	for (int i = 0; i < NBR_STRIP; ++i)
    {
    	double strip_energy = abs(strip_B[i] + distribution(generator));
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

void figure1() {
	// data for Adam2020

	// open file
    TFile *f = TFile::Open("/media/matthieu/ssd1/Geant4/Data/DataSet_1/data100k.root", "read");

    //************* variable ***************//
    const int NBR_STRIP = 254;

    const int MAX_CLUSTER_WIDTH = 3;
    const int CLUSTER_WINDOW = 5;

    const int NBR_BINS = 18;
    //const double LOWER_SIDE = 0.5;
    //const double INCREMENT = 0.5;

    // array creation
    Double_t x1[NBR_BINS] = {	2.0,
							    2.2,
							    2.4,
							    2.6,
							    2.8,
							    3.0,
							    3.5,
							    4.0,
							    4.5,
							    5.0,
							    7.5,
							    10.0,
							    15.0,
							    20.0,
							    25.0,
							    30.0,
							    35.0,
							    40.0}; // in keV !!!
    Double_t y1[NBR_BINS] = {0.};
    Double_t ex1[NBR_BINS] = {0.};
    Double_t ey1[NBR_BINS] = {0.};

    Double_t x2[NBR_BINS] = {0.};
    Double_t y2[NBR_BINS] = {0.};
    Double_t ex2[NBR_BINS] = {0.};
    Double_t ey2[NBR_BINS] = {0.};

    for (int i = 0; i < NBR_BINS; ++i)
    {
    	x1[i] = x2[i] = x1[i] * 3.62;
    }


    //************ Get the TTree ******************************************
    auto data = f->Get<TTree>("data");

    // Get the number of entries in TTree
    const int ENTRIES = data->GetEntries();
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
    auto c1 = new TCanvas("c1","c1",1920,1080);
	c1->SetTitle("Figure 1: threshold variation");
	gStyle->SetOptStat(0);
	gPad->SetGridx(1);
    gPad->SetGridy(1);
    gPad->SetLogy();

    gPad->SetTitle("Figure 1: threshold variation");

    for (int j = 0; j < NBR_BINS; ++j)
    {
	    // Create vectors for stat
	    std::vector<double> nbr_clusterA(10, 0);
	    std::vector<double> nbr_clusterB(10, 0);

	    //****************** Main loop over all entries **********************//
	    int count_loop = 0;
	    int loop_number = 0;
	    for (int k = 0; k < ENTRIES; k++)
	    {
	        // fill variables with datas from entry i
	        data->GetEntry(k);

	        std::vector<double> res_A(9, 0);
	        std::vector<double> res_B(9, 0);

	        bool stub1 = CBC2(strip_A, strip_B, res_A, res_B, MAX_CLUSTER_WIDTH, CLUSTER_WINDOW, x1[j]/1000);

	        nbr_clusterA.at(loop_number) += (double)(res_A.at(5) * res_A.at(6)) / (ENTRIES / 10);

	        nbr_clusterB.at(loop_number) += (double)(res_B.at(5) * res_B.at(6)) / (ENTRIES / 10);

	        count_loop += 1;
	        if (count_loop == ENTRIES / 10)
	        {
	        	count_loop = 0;
	        	loop_number += 1;
	        }
	    }
	    double variance1 = 0., deviation1 = 0., average1 = 0.;
		average1 = std::accumulate(nbr_clusterA.begin(), nbr_clusterA.end(), 0.0) / 10;
		for (int i = 0; i < 10; ++i){variance1 += pow((nbr_clusterA.at(i) - average1), 2);}
		variance1 /= 9;
		deviation1 = sqrt(variance1);
		cout << std::scientific << "File " << j << " mean cluster number A:\t" << average1 << "\t" << deviation1 << endl;

		y1[j] = average1;
		ey1[j] = deviation1;

		double variance2 = 0., deviation2 = 0., average2 = 0.;
		average2 = std::accumulate(nbr_clusterB.begin(), nbr_clusterB.end(), 0.0) / 10;
		for (int i = 0; i < 10; ++i){variance2 += pow((nbr_clusterB.at(i) - average2), 2);}
		variance2 /= 9;
		deviation2 = sqrt(variance2);
		cout << std::scientific << "\t" << j << " mean cluster number B:\t" << average2 << "\t" << deviation2 << endl;

		y2[j] = average2;
		ey2[j] = deviation2;
	}

    TGraphErrors *gr1 = new TGraphErrors(NBR_BINS,x1,y1,ex1,ey1);
    gr1->SetName("gr1");
    gr1->SetMarkerColor(kRed+2);
    gr1->SetMarkerStyle(20);
    gr1->SetMarkerSize(1);

    TGraphErrors *gr2 = new TGraphErrors(NBR_BINS,x2,y2,ex2,ey2);
    gr2->SetName("gr2");
    gr2->SetMarkerColor(kRed+2);
    gr2->SetMarkerStyle(25);
    gr2->SetMarkerSize(1.2);

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(gr1);
    mg->Add(gr2);
    mg->SetTitle("Figure 1: threshold variation");
    mg->Draw("AP");

    TAxis *xaxis = mg->GetXaxis();
    TAxis *yaxis = mg->GetYaxis();
    xaxis->SetTitle("Threshold [keV]");
    //xaxis->Set(25, 1.0, 3.5);
    xaxis->SetRangeUser(0., 150.);

    yaxis->SetTitle("Mean number of hits per event");
    yaxis->SetRangeUser(0.05, 30);

    c1->RedrawAxis();

    auto legend = new TLegend(0.7,0.9,0.9,0.75);
    legend->AddEntry("gr1","Top sensor","p");
    legend->AddEntry("gr2","Bottom sensor","p");
    legend->Draw();

    gPad->Modified();
    //*********************** 

    // Close file when finished
    //f.Close();
}