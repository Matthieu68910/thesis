#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <random>

std::default_random_engine generator;
std::normal_distribution<double> distribution(0.00362, 0.795); // µ = 1000 e-, s = 800 e-

bool CBC2(
	const vector<double> &strip_A, 
	const vector<double> &strip_B, 
	vector<double> &res_A, // [0:4]-> 1-5 strip wide clusters,[5]-> number of clusters (tot), [6]-> mean cluster width (tot), [7:8]-> accepted
	vector<double> &res_B,
	const int MAX_CLUSTER_WIDTH = 3,
	const int CLUSTER_WINDOW = 5,
	const double THRESHOLD = 5.1975
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
    	double strip_energy = (strip_A[i] / 0.00362) + abs(distribution(generator));
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
    	double strip_energy = (strip_B[i] / 0.00362) + abs(distribution(generator));
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

void SaveData(
	const int &k,
	Double_t x[],
	Double_t y1[],
	Double_t ey1[],
	Double_t y2[],
	Double_t ey2[],
	Double_t y3[],
	Double_t ey3[]
	){

	// open file
	ofstream myfile;
    myfile.open ("figure18_data.txt");
    myfile << "x\ty1\tey1\ty2\tey2\ty3\tey3\n";
    for (int i = 0; i < k; ++i)
    {
    	myfile 	<< std::scientific 
    			<< x[i] 
    			<< "\t" 
    			<< y1[i] 
    			<< "\t" 
    			<< ey1[i] 
    			<< "\t" 
    			<< y2[i] 
    			<< "\t" 
    			<< ey2[i]
    			<< "\t" 
    			<< y3[i] 
    			<< "\t" 
    			<< ey3[i]
    			<< "\n";
    }
    myfile.close();

    return;
}

void figure18() {
	const Int_t n = 25; // adam2020 1 strip
 
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
    Double_t y1[n] = {	0.89879,
						0.89353,
						0.89353,
						0.88827,
						0.88564,
						0.87250,
						0.86461,
						0.86461,
						0.84621,
						0.83832,
						0.82780,
						0.79100,
						0.77523,
						0.76471,
						0.74105,
						0.73053,
						0.69899,
						0.68847,
						0.67533,
						0.67007,
						0.64641,
						0.63852,
						0.60960,
						0.57543,
						0.52022};
	Double_t ex1[n] = {0.};
    Double_t ey1[n] = {0.};

    const Int_t m = 25; // adam2020 2 strips
 
	Double_t x2[m] = {	1.300,
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
    Double_t y2[m] = {	0.09302,
						0.09683,
						0.09743,
						0.10252,
						0.10343,
						0.11812,
						0.12727,
						0.12568,
						0.14430,
						0.15040,
						0.16068,
						0.19949,
						0.21547,
						0.22468,
						0.24904,
						0.25910,
						0.28931,
						0.30226,
						0.31401,
						0.31838,
						0.34249,
						0.35111,
						0.37945,
						0.41475,
						0.46736};
	Double_t ex2[m] = {0.};
    Double_t ey2[m] = {0.};

    const Int_t p = 25; // adam2020 >2 strips
 
	Double_t x3[p] = {	1.300,
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
    Double_t y3[p] = {	0.00958,
						0.01047,
						0.01064,
						0.00995,
						0.01143,
						0.01017,
						0.01009,
						0.01283,
						0.01128,
						0.01159,
						0.01154,
						0.01117,
						0.01159,
						0.01150,
						0.01207,
						0.01296,
						0.01246,
						0.01138,
						0.01265,
						0.01273,
						0.01276,
						0.01251,
						0.01285,
						0.01233,
						0.01436};
	Double_t ex3[p] = {0.};
    Double_t ey3[p] = {0.};
    					
	const Int_t k = 25; // Simulation // 1 strip
 
	Double_t x4[k] = {	1.300,
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
    Double_t y4[k] = {0.};
	Double_t ex4[k] = {0.};
    Double_t ey4[k] = {0.};

    Double_t x5[k] = {	1.300,
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
						15.000}; // 2 strips
    Double_t y5[k] = {0.};
	Double_t ex5[k] = {0.};
    Double_t ey5[k] = {0.};

    Double_t x6[k] = {	1.300,
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
						15.000}; // >2 strips
    Double_t y6[k] = {0.};
	Double_t ex6[k] = {0.};
    Double_t ey6[k] = {0.};

    //****************** Create Histo ************************************//
    auto c1 = new TCanvas("c1","c1",1000,600);
	c1->SetTitle("Mean cluster width for 2S mini-module");
	gStyle->SetOptStat(0);
	gPad->SetGridx(1);
    gPad->SetGridy(1);
    gPad->SetLogy();

    gPad->SetTitle("Mean cluster width for 2S mini-module");

    // loop on all files
    for (int j = 0; j < k; ++j)
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
	    const int MAX_CLUSTER_WIDTH = 3;
	    const int CLUSTER_WINDOW = 5;
	    const double THRESHOLD = 5.1975; // MeV -> = 6 * (1000 * 3.6 keV)
	    
	    Int_t nbrCAT, nbrCA, nbrCBT, nbrCB; // for A and B detectors
	    Double_t mCWAT, mCWBT, mCWA, mCWB;
	    Int_t CA1, CA2, CA3, CA4, CA5;
	    Int_t CB1, CB2, CB3, CB4, CB5;
	    Double_t momentum;

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
	    double nbr_cluster_1 = 0.;
	    double nbr_cluster_2 = 0.;
	    double total_clusters = 0.;
	    std::vector<double> fracC1;
	    std::vector<double> fracC2;
	    std::vector<double> fracC3;
	    for (int k = 0; k < ENTRIES; k++)
	    {
	        // fill variables with datas from entry i
	        data->GetEntry(k);

	        std::vector<double> res_A(9, 0);
	        std::vector<double> res_B(9, 0);

	        bool stub = CBC2(strip_A, strip_B, res_A, res_B, MAX_CLUSTER_WIDTH, CLUSTER_WINDOW, THRESHOLD);

	        if(!isnan(res_A.at(0))){nbr_cluster_1 += res_A.at(0);}
	        if(!isnan(res_A.at(1))){nbr_cluster_2 += res_A.at(1);}
	        if(!isnan(res_A.at(5))){total_clusters += res_A.at(5);}

	        count_loop += 1;
	        if (count_loop == ENTRIES /10)
	        {
	            count_loop = 0;
	            fracC1.push_back((double) nbr_cluster_1 / total_clusters);
	            fracC2.push_back((double) nbr_cluster_2 / total_clusters);
	            fracC3.push_back((double) 1 - ((nbr_cluster_1 + nbr_cluster_2) / total_clusters));
	            //cout << mean_cluster_width << endl;
	            nbr_cluster_1 = nbr_cluster_2 = total_clusters = 0.;
	        }
	    }
	    //********************* fig 17 computation and printing **********************************//
	    // for A
	    double variance, deviation, average;
	    average = std::accumulate(fracC1.begin(), fracC1.end(), 0.0) / 10;
	    for (int i = 0; i < 10; ++i){variance += pow((fracC1.at(i) - average), 2);}
	    variance /= 9;
	    deviation = sqrt(variance);
	    cout << std::scientific << "File " << j << " CW1:\t" << average << "\t" << deviation << endl;
	    y4[j] = average;
	    ey4[j] = deviation;

	    average = std::accumulate(fracC2.begin(), fracC2.end(), 0.0) / 10;
	    for (int i = 0; i < 10; ++i){variance += pow((fracC2.at(i) - average), 2);}
	    variance /= 9;
	    deviation = sqrt(variance);
	    cout << std::scientific << "File " << j << " CW2:\t" << average << "\t" << deviation << endl;
	    y5[j] = average;
	    ey5[j] = deviation;

	    average = std::accumulate(fracC3.begin(), fracC3.end(), 0.0) / 10;
	    for (int i = 0; i < 10; ++i){variance += pow((fracC3.at(i) - average), 2);}
	    variance /= 9;
	    deviation = sqrt(variance);
	    cout << std::scientific << "File " << j << " CW3p:\t" << average << "\t" << deviation << endl;
	    y6[j] = average;
	    ey6[j] = deviation;

	    // Close file when finished
	    f.Close();
	}   
	// Fill graphs
	TGraphErrors *gr1 = new TGraphErrors(n,x1,y1,ex1,ey1); // Adam 1
    gr1->SetName("gr1");
    gr1->SetMarkerColor(12);
    gr1->SetMarkerStyle(24);
    gr1->SetMarkerSize(1.);

    TGraphErrors *gr2 = new TGraphErrors(m,x2,y2,ex2,ey2); // Adam 2
    gr2->SetName("gr2");
    gr2->SetMarkerColor(12);
    gr2->SetMarkerStyle(25);
    gr2->SetMarkerSize(1.);

    TGraphErrors *gr3 = new TGraphErrors(p,x3,y3,ex3,ey3); // Adam >2
    gr3->SetName("gr3");
    gr3->SetMarkerColor(12);
    gr3->SetMarkerStyle(26);
    gr3->SetMarkerSize(1.);

    TGraphErrors *gr4 = new TGraphErrors(k,x4,y4,ex4,ey4); // Adam 1
    gr4->SetName("gr4");
    gr4->SetMarkerColor(kRed+2);
    gr4->SetMarkerStyle(20);
    gr4->SetMarkerSize(1.);

    TGraphErrors *gr5 = new TGraphErrors(k,x5,y5,ex5,ey5); // Adam 2
    gr5->SetName("gr5");
    gr5->SetMarkerColor(kBlue+2);
    gr5->SetMarkerStyle(21);
    gr5->SetMarkerSize(1.);

    TGraphErrors *gr6 = new TGraphErrors(k,x6,y6,ex6,ey6); // Adam >2
    gr6->SetName("gr6");
    gr6->SetMarkerColor(kGreen+2);
    gr6->SetMarkerStyle(22);
    gr6->SetMarkerSize(1.);

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(gr1);
    mg->Add(gr2);
    mg->Add(gr3);
    mg->Add(gr4);
    mg->Add(gr5);
    mg->Add(gr6);
    mg->SetTitle("");
    mg->Draw("AP");

    TAxis *xaxis = mg->GetXaxis();
    TAxis *yaxis = mg->GetYaxis();
    xaxis->SetLabelFont(42);
	xaxis->SetLabelSize(0.04);
    xaxis->SetTitle("Angle d'incidence [deg]");
    xaxis->SetTitleFont(22);
	xaxis->SetTitleSize(0.05);
	xaxis->SetTitleOffset(0.95);
    xaxis->SetRangeUser(0.5, 15);
    yaxis->SetLabelFont(42);
	yaxis->SetLabelSize(0.04);
    yaxis->SetTitle("Fraction");
    yaxis->SetTitleFont(22);
	yaxis->SetTitleSize(0.05);
	yaxis->SetTitleOffset(0.9);
	//yaxis->SetRangeUser(0.002, 1.1);

    TF1* f1 = new TF1("f1", "x", 0.5, 15);
    TGaxis* A1 = new TGaxis(0.5, yaxis->GetXmax(), 15.0, yaxis->GetXmax(), "f1", 510, "-");
    A1->SetLabelFont(42);
	A1->SetLabelSize(0.04);
	A1->SetTitle("pT [GeV]");
	A1->SetTitleFont(22);
	A1->SetTitleSize(0.05);
	A1->SetTitleOffset(0.95);
    A1->ChangeLabel(1, -1, -1, -1, -1, -1, "9.8");
    A1->ChangeLabel(2, -1, -1, -1, -1, -1, "4.9");
    A1->ChangeLabel(3, -1, -1, -1, -1, -1, "3.3");
    A1->ChangeLabel(4, -1, -1, -1, -1, -1, "2.5");
    A1->ChangeLabel(5, -1, -1, -1, -1, -1, "2.0");
    A1->ChangeLabel(6, -1, -1, -1, -1, -1, "1.6");
    A1->ChangeLabel(7, -1, -1, -1, -1, -1, "1.4");
	A1->Draw("SAME");

    auto legend = new TLegend(0.65,0.45,0.9,0.7);
    legend->AddEntry("gr1","Adam et al. - 1 strip","p");
    legend->AddEntry("gr2","Adam et al. - 2 strips","p");
    legend->AddEntry("gr3","Adam et al. - >2 strips","p");
    legend->AddEntry("gr4","Geant4 - 1 strip","ep");
    legend->AddEntry("gr5","Geant4 - 2 strips","ep");
    legend->AddEntry("gr6","Geant4 - >2 strips","ep");
    legend->Draw();

    gPad->Modified();

    c1->SaveAs("figure18.pdf");

    SaveData(k, x4, y4, ey4, y5, ey5, y6, ey6);
}