#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <random>

std::default_random_engine generator;
std::normal_distribution<double> distribution(0., 0.002891475); // µ = 0 e-, s = 2.13*375*3.62
std::uniform_real_distribution<> distribution1(0.0, 1.0);

bool CBC2(
	const vector<double> &strip_A, 
	const vector<double> &strip_B, 
	vector<double> &res_A, // [0:4]-> 1-5 strip wide clusters,[5]-> number of clusters (tot), [6]-> mean cluster width (tot), [7:8]-> accepted, [9] -> nbr hits
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

void figure13() {
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
    Double_t x1[NBR_BINS] = {	2.715000,
                                5.491088,
                                8.267175,
                                10.898010,
                                19.081020,
                                21.857108,
                                27.263080,
                                35.299751,
                                40.852334,
                                59.847831,
                                70.660726,
                                76.101993,
                                81.473621,
                                86.880000}; // in keV !!!
    Double_t y1[NBR_BINS] = {0.};
    Double_t ex1[NBR_BINS] = {0.};
    Double_t ey1[NBR_BINS] = {0.};

    Double_t x2[NBR_BINS] = {   2.86161000,
                                5.63769750,
                                8.26717500,
                                11.04326250,
                                19.22627250,
                                21.85710750,
                                27.26307975,
                                35.44581825,
                                40.85233350,
                                59.91163350,
                                70.80679275,
                                76.21330800,
                                81.61968750,
                                87.02606700}; // in keV !!!
    Double_t y2[NBR_BINS] = {0.};
    Double_t ex2[NBR_BINS] = {0.};
    Double_t ey2[NBR_BINS] = {0.};

    // adam D0
    Double_t x3[NBR_BINS] = {   2.715000,
                                5.491088,
                                8.267175,
                                10.898010,
                                19.081020,
                                21.857108,
                                27.263080,
                                35.299751,
                                40.852334,
                                59.847831,
                                70.660726,
                                76.101993,
                                81.473621,
                                86.880000}; // in keV !!!
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
    Double_t ex3[NBR_BINS] = {0.};
    Double_t ey3[NBR_BINS] = {0.};

    // adam D1
    Double_t x4[NBR_BINS] = {   2.86161000,
                                5.63769750,
                                8.26717500,
                                11.04326250,
                                19.22627250,
                                21.85710750,
                                27.26307975,
                                35.44581825,
                                40.85233350,
                                59.91163350,
                                70.80679275,
                                76.21330800,
                                81.61968750,
                                87.02606700}; // in keV !!!
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
    Double_t ex4[NBR_BINS] = {0.};
    Double_t ey4[NBR_BINS] = {0.};

    // adam D0_IRR
    Double_t x5[NBR_BINS_IRR] = {   2.73400500,
                                    5.51145000,
                                    8.17350750,
                                    10.95095250,
                                    16.38909750,
                                    19.05115500,
                                    21.82724250,
                                    32.58936675,
                                    35.36654025,
                                    40.80522825,
                                    51.68274000,
                                    70.66031850,
                                    70.67402925,
                                    76.09900650,
                                    81.53769450,
                                    86.97638250}; // in keV !!!
    Double_t y5[NBR_BINS_IRR] = {   30.1215000,
                                    6.7070400,
                                    1.9182700,
                                    1.1759900,
                                    0.9366170,
                                    0.8848140,
                                    0.8264180,
                                    0.6433920,
                                    0.6147630,
                                    0.5424320,
                                    0.3899660,
                                    0.1416410,
                                    0.1412730,
                                    0.1090240,
                                    0.0858497,
                                    0.0668366};
    Double_t ex5[NBR_BINS_IRR] = {0.};
    Double_t ey5[NBR_BINS_IRR] = {0.};

    // adam D1_IRR
    Double_t x6[NBR_BINS_IRR] = {   2.7150000,
                                    5.4300000,
                                    8.1735075,
                                    10.8600000,
                                    16.2900000,
                                    19.0050000,
                                    21.7200000,
                                    32.5800000,
                                    35.2950000,
                                    40.7250000,
                                    51.5850000,
                                    70.5900000,
                                    70.6603185,
                                    76.1019930,
                                    81.4500000,
                                    86.8800000}; // in keV !!!
    Double_t y6[NBR_BINS_IRR] = {   29.2426000,
                                    6.1077500,
                                    1.7713900,
                                    1.1762000,
                                    0.9612890,
                                    0.8955480,
                                    0.8407000,
                                    0.6526200,
                                    0.6244690,
                                    0.5483190,
                                    0.3952350,
                                    0.1473400,
                                    0.1482370,
                                    0.1148180,
                                    0.0918478,
                                    0.0727227};
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
	    std::vector<double> nbr_hitsA(10, 0);
	    std::vector<double> nbr_hitsB(10, 0);

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

	        bool stub1 = CBC2(strip_A, strip_B, res_A1, res_B1, MAX_CLUSTER_WIDTH, CLUSTER_WINDOW, x1[j]/1000);
            bool stub2 = CBC2(strip_A, strip_B, res_A2, res_B2, MAX_CLUSTER_WIDTH, CLUSTER_WINDOW, x2[j]/1000);

	        nbr_hitsA.at(loop_number) += (double) res_A1.at(9) / (ENTRIES / 10);

	        nbr_hitsB.at(loop_number) += (double) res_B2.at(9) / (ENTRIES / 10);

	        count_loop += 1;
	        if (count_loop == ENTRIES / 10)
	        {
	        	count_loop = 0;
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

		double variance2 = 0., deviation2 = 0., average2 = 0.;
		average2 = std::accumulate(nbr_hitsB.begin(), nbr_hitsB.end(), 0.0) / 10;
		for (int i = 0; i < 10; ++i){variance2 += pow((nbr_hitsB.at(i) - average2), 2);}
		variance2 /= 9;
		deviation2 = sqrt(variance2);
		cout << std::scientific << x1[j] << "\t" << average1 << "\t" << deviation1 << "\t" << average2 << "\t" << deviation2 << endl;

		y2[j] = average2;
		ey2[j] = deviation2;
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
    gr3->SetMarkerColor(kBlue+3);
    gr3->SetMarkerStyle(20);
    gr3->SetMarkerSize(1.0);

    TGraphErrors *gr4 = new TGraphErrors(NBR_BINS,x4,y4,ex4,ey4);
    gr4->SetName("gr4");
    gr4->SetMarkerColor(kBlue+3);
    gr4->SetMarkerStyle(25);
    gr4->SetMarkerSize(1.0);

    TGraphErrors *gr5 = new TGraphErrors(NBR_BINS_IRR,x5,y5,ex5,ey5);
    gr5->SetName("gr5");
    gr5->SetMarkerColor(kGreen+3);
    gr5->SetMarkerStyle(20);
    gr5->SetMarkerSize(1.0);

    TGraphErrors *gr6 = new TGraphErrors(NBR_BINS_IRR,x6,y6,ex6,ey6);
    gr6->SetName("gr6");
    gr6->SetMarkerColor(kGreen+3);
    gr6->SetMarkerStyle(25);
    gr6->SetMarkerSize(1.0);

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(gr1);
    mg->Add(gr2);
    mg->Add(gr3);
    mg->Add(gr4);
    mg->Add(gr5);
    mg->Add(gr6);
    mg->SetTitle("Figure 13: threshold variation");
    mg->Draw("AP");

    TAxis *xaxis = mg->GetXaxis();
    TAxis *yaxis = mg->GetYaxis();
    xaxis->SetTitle("Threshold [keV]");
    //xaxis->Set(25, 1.0, 3.5);
    xaxis->SetRangeUser(0., 150.);

    yaxis->SetTitle("Mean number of hits per event");
    yaxis->SetRangeUser(0.05, 100);

    c1->RedrawAxis();

    auto legend = new TLegend(0.7,0.9,0.9,0.75);
    legend->AddEntry("gr1","Geant4 - top sensor","ep");
    legend->AddEntry("gr2","Geant4 - bottom sensor","ep");
    legend->AddEntry("gr3","Adam et al. - top sensor","p");
    legend->AddEntry("gr4","Adam et al. - bottom sensor","p");
    legend->AddEntry("gr5","Adam et al. - top sensor irr.","p");
    legend->AddEntry("gr6","Adam et al. - bottom sensor irr.","p");
    legend->Draw();

    gPad->Modified();
    //*********************** 

    // Close file when finished
    //f.Close();
}