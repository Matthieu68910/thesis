#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <random>

std::default_random_engine generator;
std::normal_distribution<double> distribution(0., 0.002891475); // µ = 0 e-, s = 2.13*375*3.62
std::uniform_real_distribution<> distribution1(0.0, 1.0);

bool CBC2(
    const double &x0,
    const double &theta_i,
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
        double pos_real = x0 + (tan(theta_i) * 0.135);
        //cout << "pos_sensor = " << pos_sensor << endl;
        //cout << "pos_real = " << pos_real << endl;
        if(abs(pos_real - pos_sensor) <= (0.103923 / cos(theta_i))){match = true;}
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

void figure19() {
	// data for Adam2020

	// open file
    TFile *f = TFile::Open("/media/matthieu/ssd1/Geant4/Data/Data_figure13-16/data100k_pions.root", "read");

    const Int_t n = 25; // adam2020 non-irradiated
 
    Double_t x2[n] = {  1.300,
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
    Double_t y2[n] = {  0.9945519,
                        0.9949264,
                        0.9948899,
                        0.9950949,
                        0.9949529,
                        0.9952242,
                        0.9952743,
                        0.9952294,
                        0.9954266,
                        0.9953091,
                        0.9958590,
                        0.9953932,
                        0.9949760,
                        0.9950563,
                        0.9965728,
                        0.9956231,
                        0.9963466,
                        0.9959599,
                        0.9958317,
                        0.9965371,
                        0.9966218,
                        0.9962193,
                        0.9964897,
                        0.9957233,
                        0.9966631};
    Double_t ex2[n] = { 0.0014096,
                        0.0013615,
                        0.0009794,
                        0.0012729,
                        0.0014898,
                        0.0011438,
                        0.0015762,
                        0.0015474,
                        0.0009673,
                        0.0013345,
                        0.0011523,
                        0.0010456,
                        0.0013778,
                        0.0010189,
                        0.0014869,
                        0.0005323,
                        0.0016299,
                        0.0009557,
                        0.0009313,
                        0.0012692,
                        0.0007909,
                        0.0006010,
                        0.0008025,
                        0.0013576,
                        0.0013455};
    Double_t ey2[n] = { 0.0014096,
                        0.0013615,
                        0.0009794,
                        0.0012729,
                        0.0014898,
                        0.0011438,
                        0.0015762,
                        0.0015474,
                        0.0009673,
                        0.0013345,
                        0.0011523,
                        0.0010456,
                        0.0013778,
                        0.0010189,
                        0.0014869,
                        0.0005323,
                        0.0016299,
                        0.0009557,
                        0.0009313,
                        0.0012692,
                        0.0007909,
                        0.0006010,
                        0.0008025,
                        0.0013576,
                        0.0013455};

    const Int_t m = 17; // adam2020 irradiated
 
    Double_t x3[m] = {  5.0,
                        6.2,
                        6.9,
                        7.0,
                        7.8,
                        8.1,
                        8.5,
                        9.0,
                        9.1,
                        9.1,
                        9.2,
                        9.5,
                        9.9,
                        10.2,
                        10.5,
                        11.5,
                        12.5};
    Double_t y3[m] = {  0.9880492,
                        0.9793732,
                        0.9920439,
                        0.9890069,
                        0.9822211,
                        0.9811468,
                        0.9803337,
                        0.9768101,
                        0.9787623,
                        0.9761969,
                        0.9770920,
                        0.9737097,
                        0.9759281,
                        0.9761464,
                        0.9750419,
                        0.9759553,
                        0.9660823};
    Double_t ex3[m] = { 0.0069608,
                        0.0104303,
                        0.0027381,
                        0.0055620,
                        0.0043319,
                        0.0053377,
                        0.0057922,
                        0.0047982,
                        0.0036445,
                        0.0033799,
                        0.0085658,
                        0.0058824,
                        0.0022672,
                        0.0082825,
                        0.0044714,
                        0.0070758,
                        0.0024858};
    Double_t ey3[m] = { 0.0069608,
                        0.0104303,
                        0.0027381,
                        0.0055620,
                        0.0043319,
                        0.0053377,
                        0.0057922,
                        0.0047982,
                        0.0036445,
                        0.0033799,
                        0.0085658,
                        0.0058824,
                        0.0022672,
                        0.0082825,
                        0.0044714,
                        0.0070758,
                        0.0024858};
                        
    const Int_t index = 25; // simulation
 
    Double_t x1[index] = {  1.300,
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
    Double_t y1[index] = {0.};
    Double_t ex1[index] = {0.};
    Double_t ey1[index] = {0.};

    //****************** Create Histo ************************************//
    auto c1 = new TCanvas("c1","c1",1920,1080);
	c1->SetTitle("Figure 1: threshold variation");
	gStyle->SetOptStat(0);
	gPad->SetGridx(1);
    gPad->SetGridy(1);
    //gPad->SetLogy();

    gPad->SetTitle("Figure 1: threshold variation");

    // loop on all files
    for (int j = 0; j < n; ++j)
    {
        //CopyFile(j);

        string file_path = "/media/matthieu/ssd1/Geant4/Data/Data_figure17-19/data_";
        file_path += std::to_string(j);
        file_path += ".root";
        char const *pchar = file_path.c_str();
        TFile f(pchar, "read");
        //cout << "File " << j << " opened:" << endl;

        //************* variable ***************//
        const int NBR_STRIP = 254;
        const int MAX_CLUSTER_WIDTH = 5;
        const int CLUSTER_WINDOW = 5;
        const double THRESHOLD = 0.019005; // (14 * 375 * 3.62) / 1000000

        Int_t nbrCAT, nbrCA, nbrCBT, nbrCB; // for A and B detectors
        Double_t mCWAT, mCWBT, mCWA, mCWB;
        Int_t CA1, CA2, CA3, CA4, CA5;
        Int_t CB1, CB2, CB3, CB4, CB5;
        Double_t momentum;

        //************ Get the TTree ******************************************
        auto data = f.Get<TTree>("data");

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
        double theta_i = 0.;
        data->SetBranchAddress("x0", &x0);
        data->SetBranchAddress("theta_i", &theta_i);

        // Create vectors for stat
        std::vector<double> nbr_hitsA(10, 0);
        std::vector<double> nbr_eventA(10, 0);

        //****************** Main loop over all entries **********************//
        int count_loop = 0;
        int loop_number = 0;
        for (int k = 0; k < ENTRIES; k++)
        {
            // fill variables with datas from entry i
            data->GetEntry(k);

            std::vector<double> res_A1(10, 0);
            std::vector<double> res_B1(10, 0);

            bool stub1 = CBC2(x0, theta_i, strip_A, strip_B, res_A1, res_B1, MAX_CLUSTER_WIDTH, CLUSTER_WINDOW, THRESHOLD);
            //bool stub2 = CBC2(strip_A, strip_B, res_A2, res_B2, MAX_CLUSTER_WIDTH, CLUSTER_WINDOW, x2[j]/1000);

            if(res_A1.at(5) == 1 && res_B1.at(5) == 1)
            {
                if(res_A1.at(9) == 1){nbr_hitsA.at(loop_number) += 1.0;}
                nbr_eventA.at(loop_number) += 1.0;
                //cout << "hint" << endl;
            }

            //nbr_hitsB.at(loop_number) += (double) res_B2.at(5) / (ENTRIES / 10);

            count_loop += 1;
            if (count_loop == ENTRIES / 10)
            {
                count_loop = 0;
                nbr_hitsA.at(loop_number) /= nbr_eventA.at(loop_number);
                loop_number += 1;
            }
        }
        double variance1 = 0., deviation1 = 0., average1 = 0.;
        average1 = std::accumulate(nbr_hitsA.begin(), nbr_hitsA.end(), 0.0) / 10;
        for (int i = 0; i < 10; ++i){variance1 += pow((nbr_hitsA.at(i) - average1), 2);}
        variance1 /= 9;
        deviation1 = sqrt(variance1);

        y1[j] = average1;
        ey1[j] = deviation1;

        // Close file when finished
        f.Close();
    }  

    TGraphErrors *gr1 = new TGraphErrors(index,x1,y1,ex1,ey1);
    gr1->SetName("gr1");
    gr1->SetMarkerColor(kRed+2);
    gr1->SetMarkerStyle(20);
    gr1->SetMarkerSize(1.0);

    TGraphErrors *gr2 = new TGraphErrors(n,x2,y2,ex2,ey2);
    gr2->SetName("gr2");
    gr2->SetMarkerColor(kBlue+2);
    gr2->SetMarkerStyle(20);
    gr2->SetMarkerSize(1.0);

    TGraphErrors *gr3 = new TGraphErrors(m,x3,y3,ex3,ey3);
    gr3->SetName("gr3");
    gr3->SetMarkerColor(kGreen+3);
    gr3->SetMarkerStyle(20);
    gr3->SetMarkerSize(1.0);

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(gr1);
    mg->Add(gr2);
    mg->Add(gr3);
    mg->SetTitle("Figure 19: Cluster efficiency");
    mg->Draw("AP");

    TAxis *xaxis = mg->GetXaxis();
    TAxis *yaxis = mg->GetYaxis();
    xaxis->SetTitle("Angle [deg]");
    //xaxis->Set(25, 1.0, 3.5);
    xaxis->SetRangeUser(0., 16.5);

    yaxis->SetTitle("Cluster efficiency");
    yaxis->SetRangeUser(0.92, 1.01);

    c1->RedrawAxis();

    auto legend = new TLegend(0.7,0.9,0.9,0.75);
    legend->AddEntry("gr1","Geant4","ep");
    legend->AddEntry("gr2","Adam et al. - non-irr.","ep");
    legend->AddEntry("gr3","Adam et al. - irradiated","ep");
    legend->Draw();

    gPad->Modified();
}