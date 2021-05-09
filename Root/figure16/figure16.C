#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <random>

std::default_random_engine generator;
std::normal_distribution<double> distribution(0., 0.002891475); // µ = 0 e-, s = 2.13*375*3.62
std::uniform_real_distribution<> distribution1(0.0, 1.0);

bool CBC2(
    const double &x0,
	const vector<double> &strip_A, 
	const vector<double> &strip_B, 
	vector<double> &res_A, // [0:4]-> 1-5 strip wide clusters,[5]-> number of clusters (tot), [6]-> mean cluster width (tot), [7:8]-> accepted, [9] -> match
	vector<double> &res_B,
	const int MAX_CLUSTER_WIDTH = 3,
	const int CLUSTER_WINDOW = 5,
	double THRESHOLD = 0.0222,
    double kill_value = 0.04331
	){

	const int NBR_STRIP = strip_A.size();

	// Clusters in sensor A
	std::vector<double> clus_pos_A;
	std::vector<double> clus_size_A;
	bool inside = false;
	int size = 0;
    int nbr_hits_A = 0;
    bool match = false;
	// Loop on sensor A strips
	for (int i = 0; i < NBR_STRIP; ++i)
    {
    	double strip_energy = strip_A[i] + abs(distribution(generator));
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
    	double strip_energy = strip_B[i] + abs(distribution(generator));
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

    for (int b = 0; b < clus_pos_A.size(); ++b)
    {
        double pos_sensor = (clus_pos_A.at(b) - (NBR_STRIP / 2)) * 0.09;
        //cout << "pos_sensor = " << pos_sensor << endl;
        //cout << "x0 = " << x0 << endl;
        if(abs(x0 - pos_sensor) <= 0.103923){match = true;}
        //cout << "match : " << match << endl;
    }

    if(match){res_A.at(9) = 1;}else{res_A.at(9) = 0;}

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

void figure16() {
	// data for Adam2020

	// open file
    TFile *f = TFile::Open("/media/matthieu/ssd1/Geant4/Data/Data_figure13-16/data100k_pions.root", "read");

    //************* variable ***************//
    const int NBR_STRIP = 254;

    const int MAX_CLUSTER_WIDTH = 3;
    const int CLUSTER_WINDOW = 5;

    const int NBR_BINS = 14;
    const int NBR_BINS1 = 16;

    // G4
    Double_t x1[NBR_BINS] = {	2.715,
                                5.430,
                                8.145,
                                10.860,
                                19.005,
                                21.720,
                                27.150,
                                35.295,
                                40.725,
                                59.730,
                                70.590,
                                76.020,
                                81.450,
                                86.880}; // in keV !!!
    Double_t y1[NBR_BINS] = {0.};
    Double_t ex1[NBR_BINS] = {0.};
    Double_t ey1[NBR_BINS] = {0.};

    // Adam
    Double_t x2[NBR_BINS] = {   2.715,
                                5.430,
                                8.145,
                                10.860,
                                19.005,
                                21.720,
                                27.150,
                                35.295,
                                40.725,
                                59.730,
                                70.590,
                                76.020,
                                81.450,
                                86.880}; // in keV !!!
    Double_t y2[NBR_BINS] = {   0.60058,
                                0.82031,
                                0.91971,
                                0.96250,
                                0.98769,
                                0.98777,
                                0.99033,
                                0.98455,
                                0.97596,
                                0.78219,
                                0.45542,
                                0.29548,
                                0.18861,
                                0.12298};
    Double_t ex2[NBR_BINS] = {  0.034517,
                                0.048871,
                                0.029415,
                                0.012339,
                                0.002080,
                                0.002063,
                                0.001988,
                                0.002542,
                                0.009636,
                                0.093510,
                                0.110278,
                                0.078436,
                                0.057550,
                                0.038654};
    Double_t ey2[NBR_BINS] = {  0.034517,
                                0.048871,
                                0.029415,
                                0.012339,
                                0.002080,
                                0.002063,
                                0.001988,
                                0.002542,
                                0.009636,
                                0.093510,
                                0.110278,
                                0.078436,
                                0.057550,
                                0.038654};

    // Adam IRR
    Double_t x3[NBR_BINS1] = {   2.715,
                                5.430,
                                8.145,
                                10.860,
                                16.290,
                                19.005,
                                21.720,
                                32.692,
                                35.295,
                                40.725,
                                51.585,
                                70.590,
                                70.590,
                                76.020,
                                81.450,
                                86.880}; // in keV !!!
    Double_t y3[NBR_BINS1] = {  0.87545,
                                0.96947,
                                0.97383,
                                0.97342,
                                0.96847,
                                0.93369,
                                0.91345,
                                0.74739,
                                0.71188,
                                0.62001,
                                0.43720,
                                0.15504,
                                0.15474,
                                0.11836,
                                0.08795,
                                0.06827};
    Double_t ex3[NBR_BINS1] = { 0.010911,
                                0.005027,
                                0.006544,
                                0.008014,
                                0.014938,
                                0.058135,
                                0.053471,
                                0.085289,
                                0.084985,
                                0.095667,
                                0.114933,
                                0.054585,
                                0.055896,
                                0.038783,
                                0.025682,
                                0.026984};
    Double_t ey3[NBR_BINS1] = { 0.010911,
                                0.005027,
                                0.006544,
                                0.008014,
                                0.014938,
                                0.058135,
                                0.053471,
                                0.085289,
                                0.084985,
                                0.095667,
                                0.114933,
                                0.054585,
                                0.055896,
                                0.038783,
                                0.025682,
                                0.026984};


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
    double x0 = 0.;
    data->SetBranchAddress("x0", &x0);

    //****************** Create Histo ************************************//
    auto c1 = new TCanvas("c1","c1",1920,1080);
	c1->SetTitle("Figure 1: threshold variation");
	gStyle->SetOptStat(0);
	gPad->SetGridx(1);
    gPad->SetGridy(1);
    //gPad->SetLogy();

    gPad->SetTitle("Figure 1: threshold variation");

    for (int j = 0; j < NBR_BINS; ++j)
    {
	    // Create vectors for stat
	    std::vector<double> nbr_hitsA(10, 0);
	    //std::vector<double> nbr_hitsB(10, 0);

	    //****************** Main loop over all entries **********************//
	    int count_loop = 0;
	    int loop_number = 0;
	    for (int k = 0; k < ENTRIES; k++)
	    {
	        // fill variables with datas from entry i
	        data->GetEntry(k);

	        std::vector<double> res_A1(10, 0);
	        std::vector<double> res_B1(10, 0);
            //std::vector<double> res_A2(10, 0);
            //std::vector<double> res_B2(10, 0);

	        bool stub1 = CBC2(x0, strip_A, strip_B, res_A1, res_B1, MAX_CLUSTER_WIDTH, CLUSTER_WINDOW, x1[j]/1000);
            //bool stub2 = CBC2(strip_A, strip_B, res_A2, res_B2, MAX_CLUSTER_WIDTH, CLUSTER_WINDOW, x2[j]/1000);

	        if(res_A1.at(9) == 1)
            {
                nbr_hitsA.at(loop_number) += (double) 1.0 / (ENTRIES / 10);
                //cout << "hint" << endl;
            }

	        //nbr_hitsB.at(loop_number) += (double) res_B2.at(5) / (ENTRIES / 10);

	        count_loop += 1;
	        if (count_loop == ENTRIES / 10)
	        {
	        	count_loop = 0;
	        	loop_number += 1;
                cout << loop_number << endl;
	        }
	    }
	    double variance1 = 0., deviation1 = 0., average1 = 0.;
		average1 = std::accumulate(nbr_hitsA.begin(), nbr_hitsA.end(), 0.0) / 10;
		for (int i = 0; i < 10; ++i){variance1 += pow((nbr_hitsA.at(i) - average1), 2);}
		variance1 /= 9;
		deviation1 = sqrt(variance1);
		

		y1[j] = average1;
		ey1[j] = deviation1;

		/*double variance2 = 0., deviation2 = 0., average2 = 0.;
		average2 = std::accumulate(nbr_hitsB.begin(), nbr_hitsB.end(), 0.0) / 10;
		for (int i = 0; i < 10; ++i){variance2 += pow((nbr_hitsB.at(i) - average2), 2);}
		variance2 /= 9;
		deviation2 = sqrt(variance2);
		cout << std::scientific << x1[j] << "\t" << average1 << "\t" << deviation1 << "\t" << average2 << "\t" << deviation2 << endl;

		y2[j] = average2;
		ey2[j] = deviation2;*/
	}

    TGraphErrors *gr1 = new TGraphErrors(NBR_BINS,x1,y1,ex1,ey1);
    gr1->SetName("gr1");
    gr1->SetMarkerColor(kRed+2);
    gr1->SetMarkerStyle(20);
    gr1->SetMarkerSize(1.0);

    TGraphErrors *gr2 = new TGraphErrors(NBR_BINS,x2,y2,ex2,ey2);
    gr2->SetName("gr2");
    gr2->SetMarkerColor(kBlue+2);
    gr2->SetMarkerStyle(20);
    gr2->SetMarkerSize(1.0);

    TGraphErrors *gr3 = new TGraphErrors(NBR_BINS1,x3,y3,ex3,ey3);
    gr3->SetName("gr3");
    gr3->SetMarkerColor(kGreen+3);
    gr3->SetMarkerStyle(20);
    gr3->SetMarkerSize(1.0);

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(gr1);
    mg->Add(gr2);
    mg->Add(gr3);
    mg->SetTitle("Figure 16: threshold variation");
    mg->Draw("AP");

    TAxis *xaxis = mg->GetXaxis();
    TAxis *yaxis = mg->GetYaxis();
    xaxis->SetTitle("Threshold [keV]");
    //xaxis->Set(25, 1.0, 3.5);
    xaxis->SetRangeUser(0., 150.);

    yaxis->SetTitle("Cluster efficiency");
    yaxis->SetRangeUser(0., 1.2);

    c1->RedrawAxis();

    auto legend = new TLegend(0.7,0.9,0.9,0.75);
    legend->AddEntry("gr1","Geant4","ep");
    legend->AddEntry("gr2","Adam et al. - non-irr.","ep");
    legend->AddEntry("gr3","Adam et al. - irradiated","ep");
    legend->Draw();

    gPad->Modified();
    //*********************** 

    // Close file when finished
    //f.Close();
}