#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <random>

std::default_random_engine generator;
std::normal_distribution<double> distribution(0., 0.795); // µ = 0 e-, s = 2.13*375*3.62
std::uniform_real_distribution<> distribution1(0.0, 1.0);

bool CBC2(
    const double &x0,
	const vector<double> &strip_A, 
	const vector<double> &strip_B, 
	vector<double> &res_A, // [0:4]-> 1-5 strip wide clusters,[5]-> number of clusters (tot), [6]-> mean cluster width (tot), [7:8]-> accepted, [9] -> match
	vector<double> &res_B,
	const int MAX_CLUSTER_WIDTH = 3,
	const int CLUSTER_WINDOW = 5,
	double THRESHOLD = 5.1975,
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
    	double strip_energy = (strip_A[i] / 0.00362) + abs(distribution(generator));
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
    	double strip_energy = (strip_B[i] / 0.00362) + abs(distribution(generator));
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
        if(abs(x0 - pos_sensor) <= 0.09){match = true;}
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
    Double_t x1[NBR_BINS] = {	0.6975,
                                1.4475,
                                2.1975,
                                2.9475,
                                5.1975,
                                5.9475,
                                7.4475,
                                9.6975,
                                11.1975,
                                16.4475,
                                19.4475,
                                20.9475,
                                22.4475,
                                23.9475}; // in keV !!!
    Double_t y1[NBR_BINS] = {0.};
    Double_t ex1[NBR_BINS] = {0.};
    Double_t ey1[NBR_BINS] = {0.};

    // Adam
    Double_t x2[NBR_BINS] = {   0.6975,
                                1.4475,
                                2.1975,
                                2.9475,
                                5.1975,
                                5.9475,
                                7.4475,
                                9.6975,
                                11.1975,
                                16.4475,
                                19.4475,
                                20.9475,
                                22.4475,
                                23.9475}; // in keV !!!
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
    Double_t x3[NBR_BINS1] = {  0.6975,
                                1.4475,
                                2.1975,
                                2.9475,
                                4.4475,
                                5.1975,
                                5.9475,
                                8.9783,
                                9.6975,
                                11.1975,
                                14.1975,
                                19.4475,
                                19.4475,
                                20.9475,
                                22.4475,
                                23.9475}; // in keV !!!
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
    auto c1 = new TCanvas("c1","c1",1000,600);
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

	        bool stub1 = CBC2(x0, strip_A, strip_B, res_A1, res_B1, MAX_CLUSTER_WIDTH, CLUSTER_WINDOW, x1[j]);
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

    for (int i = 0; i < NBR_BINS; ++i)
    {
        x1[i] *= 1000.;
        x2[i] *= 1000.;
    }
    for (int i = 0; i < NBR_BINS1; ++i)
    {
        x3[i] *= 1000.;
    }

    TGraphErrors *gr1 = new TGraphErrors(NBR_BINS,x1,y1,ex1,ey1);
    gr1->SetName("gr1");
    gr1->SetMarkerColor(kRed+2);
    gr1->SetMarkerStyle(20);
    gr1->SetMarkerSize(1.0);

    TGraphErrors *gr2 = new TGraphErrors(NBR_BINS,x2,y2,ex2,ey2);
    gr2->SetName("gr2");
    gr2->SetMarkerColor(1);
    gr2->SetMarkerStyle(21);
    gr2->SetMarkerSize(1.0);

    TGraphErrors *gr3 = new TGraphErrors(NBR_BINS1,x3,y3,ex3,ey3);
    gr3->SetName("gr3");
    gr3->SetMarkerColor(1);
    gr3->SetMarkerStyle(22);
    gr3->SetMarkerSize(1.0);

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(gr1);
    mg->Add(gr2);
    mg->Add(gr3);
    mg->SetTitle("");
    mg->Draw("AP");

    TAxis *xaxis = mg->GetXaxis();
    TAxis *yaxis = mg->GetYaxis();
    xaxis->SetMaxDigits(3);
    xaxis->SetLabelFont(42);
    xaxis->SetLabelSize(0.04);
    xaxis->SetTitle("Seuil [e-]");
    xaxis->SetTitleFont(22);
    xaxis->SetTitleSize(0.05);
    xaxis->SetTitleOffset(0.95);
    xaxis->SetRangeUser(0., 25000);
    yaxis->SetLabelFont(42);
    yaxis->SetLabelSize(0.04);
    yaxis->SetTitle("Cluster efficiency");
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

    auto legend = new TLegend(0.1,0.4,0.4,0.1);
    legend->AddEntry("gr1","Geant4","ep");
    legend->AddEntry("gr2","Adam et al. - non-irr.","ep");
    legend->AddEntry("gr3","Adam et al. - irr.","ep");
    legend->Draw();

    gPad->Modified();
    //*********************** 
    c1->SaveAs("figure16.pdf");
    // Close file when finished
    //f.Close();
}