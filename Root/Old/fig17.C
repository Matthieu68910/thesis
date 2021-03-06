#include "TFile.h"
#include "TTree.h"
#include <iostream>

bool CBC2(
	const vector<double> &strip_A, 
	const vector<double> &strip_B, 
	vector<double> &res_A, // [0:4]-> 1-5 strip wide clusters,[5]-> number of clusters (tot), [6]-> mean cluster width (tot), [7:8]-> accepted
	vector<double> &res_B,
	const int MAX_CLUSTER_WIDTH = 3,
	const int CLUSTER_WINDOW = 5,
	const double THRESHOLD = 0.0144
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
        if (strip_A[i] < THRESHOLD && !inside)        
        {} else if (strip_A[i] < THRESHOLD && inside)
        {
	        clus_size_A.push_back(size);
	        if(size <= 5) res_A.at(size - 1) += 1;
	        clus_pos_A.push_back(floor((i - 1) - (size / 2) + 0.5));
	        size = 0;
	        inside = false;
        } else if (strip_A[i] >= THRESHOLD && !inside)
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
        } else if (strip_A[i] >= THRESHOLD && inside)
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
        if (strip_B[i] < THRESHOLD && !inside)        
        {} else if (strip_B[i] < THRESHOLD && inside)
        {
	        clus_size_B.push_back(size);
	        if(size <= 5) res_B.at(size - 1) += 1;
	        clus_pos_B.push_back(floor((i - 1) - (size / 2) + 0.5));
	        size = 0;
	        inside = false;
        } else if (strip_B[i] >= THRESHOLD && !inside)
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
        } else if (strip_B[i] >= THRESHOLD && inside)
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

void CopyFile(int& n){
    string cmd = ".!rootcp -r /media/matthieu/ssd1/Geant4/Data/fig17/test_" + std::to_string(n);
    cmd += ".root /media/matthieu/ssd1/Geant4/Data/fig17/test_";
    cmd += std::to_string(n);
    cmd += "_new.root";      
    char const *pchar = cmd.c_str();
    gROOT->ProcessLine(pchar);
}

void fig17() {
	const Int_t n = 25; // adam2020
 
	Double_t x1[n] = {	1.300,
						1.500,
						1.600,
						1.800,
						2.000,
						2.700,
						3.100,
						3.100,
						3.800,
						4.168,
						4.600,
						6.200,
						6.600,
						7.000,
						7.800,
						8.100,
						9.200,
						9.900,
						10.200,
						10.500,
						11.000,
						11.500,
						12.500,
						13.300,
						15.000};
    Double_t y1[n] = {	1.11851,
						1.12435,
						1.12561,
						1.12936,
						1.13301,
						1.14403,
						1.15330,
						1.15349,
						1.17478,
						1.18285,
						1.19106,
						1.22892,
						1.24545,
						1.25452,
						1.27972,
						1.29204,
						1.32176,
						1.33246,
						1.34754,
						1.35009,
						1.37482,
						1.38456,
						1.41346,
						1.44631,
						1.50565};
	Double_t ex1[n] = {0.};
    Double_t ey1[n] = {0.};
    					
	const Int_t k = 31; // simulation
 
	Double_t x2[k] = {	0.0,
						0.5,
						1.0,
						1.5,
						2.0,
						2.5,
						3.0,
						3.5,
						4.0,
						4.5,
						5.0,
						5.5,
						6.0,
						6.5,
						7.0,
						7.5,
						8.0,
						8.5,
						9.0,
						9.5,
						10.0,
						10.5,
						11.0,
						11.5,
						12.0,
						12.5,
						13.0,
						13.5,
						14.0,
						14.5,
						15.0};
    Double_t y2[k] = {0.};
    Double_t ex2[k] = {0.};
    Double_t ey2[k] = {0.};

    //****************** Create Histo ************************************//
    auto c1 = new TCanvas("c1","c1",1920,1080);
	c1->SetTitle("Mean cluster width for 2S mini-module");
	gStyle->SetOptStat(0);
	gPad->SetGridx(1);
    gPad->SetGridy(1);

    gPad->SetTitle("Mean cluster width for 2S mini-module");

    // loop on all files
    for (int j = 0; j < k; ++j)
    {
        //CopyFile(j);

        string file_path = "/media/matthieu/ssd1/Geant4/Data/fig17/test_";
        file_path += std::to_string(j);
        file_path += ".root";
        char const *pchar = file_path.c_str();
        TFile f(pchar, "read");
        //cout << "File " << j << " opened:" << endl;

	    //************* variable ***************//
	    const int NBR_STRIP = 101;
	    const int MAX_CLUSTER_WIDTH = 3;
	    const int CLUSTER_WINDOW = 5;
	    const double THRESHOLD = 0.017; // MeV -> = 6 * (1000 * 3.6 keV)
	    
	    Int_t nbrCAT, nbrCA, nbrCBT, nbrCB; // for A and B detectors
	    Double_t mCWAT, mCWBT, mCWA, mCWB;
	    Int_t CA1, CA2, CA3, CA4, CA5;
	    Int_t CB1, CB2, CB3, CB4, CB5;
	    Double_t momentum;

	    Bool_t corr_stub = false; // for general

	    //************ Get the TTree ******************************************
	    auto data = f.Get<TTree>("data");

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
	    data->SetBranchAddress("momentum", &momentum);

	    //****************** Main loop over all entries **********************//
	    int count_loop = 0;
	    bool stop = true;
	    double mean_cluster_width = 0.;
	    std::vector<double> mClusWidth;
	    for (int k = 0; k < ENTRIES; k++)
	    {
	        // fill variables with datas from entry i
	        data->GetEntry(k);

	        std::vector<double> res_A(9, 0);
	        std::vector<double> res_B(9, 0);

	        bool stub = CBC2(strip_A, strip_B, res_A, res_B, MAX_CLUSTER_WIDTH, CLUSTER_WINDOW, THRESHOLD);

	        CA1 = (int)res_A.at(0); // for A
	        CA2 = (int)res_A.at(1);
	        CA3 = (int)res_A.at(2);
	        CA4 = (int)res_A.at(3);
	        CA5 = (int)res_A.at(4);
	        nbrCAT = (int)res_A.at(5);
	        mCWAT = (double)res_A.at(6);
	        nbrCA = (int)res_A.at(7);
	        mCWA = (double)res_A.at(8);

	        CB1 = (int)res_B.at(0); // for B
	        CB2 = (int)res_B.at(1);
	        CB3 = (int)res_B.at(2);
	        CB4 = (int)res_B.at(3);
	        CB5 = (int)res_B.at(4);
	        nbrCBT = (int)res_B.at(5);
	        mCWBT = (double)res_B.at(6);
	        nbrCB = (int)res_B.at(7);
	        mCWB = (double)res_B.at(8);

	        corr_stub = stub; // for general

	        // fig 17 (mean cluster width)
	        //if(nbrCAT != 0 && nbrCBT != 0){mean_cluster_width += (double) (mCWAT*nbrCAT + mCWBT*nbrCBT) / (nbrCAT + nbrCBT);}
	        if(nbrCAT != 0){mean_cluster_width += (double) mCWAT;}
	        if (isnan(mean_cluster_width) && stop)
	        {
	        	cout << "at " << k << " see: " << mean_cluster_width << endl;
	        	cout << mCWAT << nbrCAT << mCWBT << nbrCBT << endl;
	        	stop = false;
	        }
	        if (!isnan(mean_cluster_width) && !stop)
	        {
	        	stop = true;
	        }
	        count_loop += 1;
	        if (count_loop == ENTRIES /10)
	        {
	            count_loop = 0;
	            mClusWidth.push_back((double) mean_cluster_width / (ENTRIES / 10));
	            //cout << mean_cluster_width << endl;
	            mean_cluster_width = 0.;
	        }
	    }
	    //********************* fig 17 computation and printing **********************************//
	    // for A
	    double variance, deviation, average;
	    average = std::accumulate(mClusWidth.begin(), mClusWidth.end(), 0.0) / 10;
	    for (int i = 0; i < 10; ++i){variance += pow((mClusWidth.at(i) - average), 2);}
	    variance /= 9;
	    deviation = sqrt(variance);
	    cout << std::scientific << "File " << j << " mean cluster width:\t" << average << "\t" << deviation << endl;

	    //*********************** 
	    y2[j] = average;
	    ey2[j] = deviation;

	    // Close file when finished
	    f.Close();
	}   
	// Fill graphs
	TGraphErrors *gr1 = new TGraphErrors(n,x1,y1,ex1,ey1);
    gr1->SetName("gr1");
    gr1->SetMarkerColor(2);
    gr1->SetMarkerStyle(20);
    gr1->SetMarkerSize(1.2);

    TGraphErrors *gr2 = new TGraphErrors(k,x2,y2,ex2,ey2);
    gr2->SetName("gr2");
    gr2->SetMarkerColor(4);
    gr2->SetMarkerStyle(20);
    gr2->SetMarkerSize(1.2);

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(gr1);
    mg->Add(gr2);
    mg->SetTitle("Mean cluster width for 2S mini-module");
    mg->Draw("AP");

    TAxis *xaxis = mg->GetXaxis();
    TAxis *yaxis = mg->GetYaxis();
    xaxis->SetTitle("Incident angle [deg]");
    xaxis->SetRangeUser(0, 15);
    yaxis->SetTitle("Mean cluster width [strip]");

    auto legend = new TLegend(0.1,0.7,0.35,0.9);
    legend->SetHeader("Legend","C");
    legend->AddEntry("gr1","Adam et all. 2020","p");
    legend->AddEntry("gr2","Geant4","ep");
    legend->Draw();

    gPad->Modified();
}