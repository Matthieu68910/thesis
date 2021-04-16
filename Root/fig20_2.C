#include "TFile.h"
#include "TTree.h"
#include <iostream>

bool CBC2(
	const vector<double> &strip_A, 
	const vector<double> &strip_B, 
	vector<double> &res_A, // [0:4]-> 1-5 strip wide clusters,[5]-> number of clusters (tot), [6]-> mean cluster width (tot), [7:8]-> accepted
	vector<double> &res_B,
	const int MAX_CLUSTER_WIDTH = 2,
	const int CLUSTER_WINDOW = 5,
	const double THRESHOLD = 0.0216
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
    string cmd = ".!rootcp -r /media/matthieu/ssd1/Geant4/Data/test_" + std::to_string(n);
    cmd += ".root /media/matthieu/ssd1/Geant4/Data/test_";
    cmd += std::to_string(n);
    cmd += "_new.root";      
    char const *pchar = cmd.c_str();
    gROOT->ProcessLine(pchar);
}

void fig20_2() {
	// data for Adam2020
	const Int_t n = 24;
 
	Double_t x1[n] = {	1.32089,
						1.48578,
						1.57853,
						1.71422,
						1.79152,
						1.87568,
						1.93064,
						1.98904,
						2.13847,
						2.42703,
						2.51978,
						2.80491,
						2.97495,
						3.1656,
						4.26316,
						4.66851,
						5.15975,
						6.32257,
						7.25867,
						9.79902,
						10.888,
						12.2483,
						13.0642,
						15.0738};

    Double_t y1[n] = {	0.00162849,
						0.00298557,
						0.00298557,
						0.02605590,
						0.24183200,
						0.49832000,
						0.66659800,
						0.81587700,
						0.97872600,
						0.98279800,
						0.98144100,
						0.98415500,
						0.98279800,
						0.98415500,
						0.98415500,
						0.98279800,
						0.98415500,
						0.98415500,
						0.98415500,
						0.98279800,
						0.98415500,
						0.98415500,
						0.98415500,
						0.98551200};

    Double_t ex1[n] = {0.};
    Double_t ey1[n] = {0.};

	// open file
    TFile *f = TFile::Open("/media/matthieu/ssd1/Geant4/Data/data100M.root", "read");

    //************* variable ***************//
    const int NBR_STRIP = 101;

    const int MAX_CLUSTER_WIDTH1 = 1;
    const int CLUSTER_WINDOW1 =5;
    const double THRESHOLD1 = 0.0216; // MeV -> = 4 * (1000 * 3.6 keV)

    const int MAX_CLUSTER_WIDTH2 = 2;
    const int CLUSTER_WINDOW2 =5;
    const double THRESHOLD2 = 0.0216; // MeV -> = 4 * (1000 * 3.6 keV)

    const int MAX_CLUSTER_WIDTH3 = 3;
    const int CLUSTER_WINDOW3 =5;
    const double THRESHOLD3 = 0.0216; // MeV -> = 4 * (1000 * 3.6 keV)

    const int NBR_BINS = 500;
    
    /*Int_t nbrCAT, nbrCA, nbrCBT, nbrCB; // for A and B detectors
    Double_t mCWAT, mCWBT, mCWA, mCWB;
    Int_t CA1, CA2, CA3, CA4, CA5;
    Int_t CB1, CB2, CB3, CB4, CB5;*/
    Double_t momentum;

    Bool_t corr_stub = false; // for general

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
    data->SetBranchAddress("momentum", &momentum);

    //****************** Create Histo ************************************//
    auto c1 = new TCanvas("c1","c1",1920,1080);
	c1->SetTitle("Stub efficiency for 2S mini-module");
	gStyle->SetOptStat(0);
	gPad->SetGridx(1);
    gPad->SetGridy(1);

    TH1* h1 = new TH1D("h1", "Stub efficiency for 2S mini-module", NBR_BINS, 1.0, 3.5);
    h1->SetName("h1");
    TH1* h2 = new TH1D("h2", "Stub efficiency 2", NBR_BINS, 1.0, 3.5);
    h2->SetName("h2");
    TH1* h3 = new TH1D("h3", "Stub efficiency 3", NBR_BINS, 1.0, 3.5);
    h3->SetName("h3");

    gPad->SetTitle("Stub efficiency for 2S mini-module");

    // Create vectors for bin content
    std::vector<double> nbr_stub1(NBR_BINS+2, 0);
    std::vector<double> nbr_event1(NBR_BINS+2, 0);
    std::vector<double> nbr_stub2(NBR_BINS+2, 0);
    std::vector<double> nbr_event2(NBR_BINS+2, 0);
    std::vector<double> nbr_stub3(NBR_BINS+2, 0);
    std::vector<double> nbr_event3(NBR_BINS+2, 0);

    //****************** Main loop over all entries **********************//
    int count_loop = 0;
    for (int k = 0; k < ENTRIES; k++)
    {
        // fill variables with datas from entry i
        data->GetEntry(k);

        std::vector<double> res_A(9, 0);
        std::vector<double> res_B(9, 0);

        bool stub1 = CBC2(strip_A, strip_B, res_A, res_B, MAX_CLUSTER_WIDTH1, CLUSTER_WINDOW1, THRESHOLD1);
        bool stub2 = CBC2(strip_A, strip_B, res_A, res_B, MAX_CLUSTER_WIDTH2, CLUSTER_WINDOW2, THRESHOLD2);
        bool stub3 = CBC2(strip_A, strip_B, res_A, res_B, MAX_CLUSTER_WIDTH3, CLUSTER_WINDOW3, THRESHOLD3);


        /*CA1 = (int)res_A.at(0); // for A
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
        mCWB = (double)res_B.at(8);*/

        //corr_stub = stub; // for general
        if(stub1){nbr_stub1.at(h1->FindBin(momentum)) += 1.;}
        nbr_event1.at(h1->FindBin(momentum)) += 1.;

        if(stub2){nbr_stub2.at(h2->FindBin(momentum)) += 1.;}
        nbr_event2.at(h2->FindBin(momentum)) += 1.;

        if(stub3){nbr_stub3.at(h3->FindBin(momentum)) += 1.;}
        nbr_event3.at(h3->FindBin(momentum)) += 1.;

        count_loop += 1;
        if (count_loop == ENTRIES / 10)
        {
        	count_loop = 0;
        	cout << "*" << endl;
        }
    }
    double error = 0., content = 0.;
    for (int i = 1; i <= NBR_BINS; ++i)
    {
    	content = nbr_stub1.at(i) / nbr_event1.at(i);
    	h1->SetBinContent(i, content);
    	//error = sqrt(nbr_event1.at(i)) / nbr_event1.at(i); //* content;
    	//h1->SetBinError(i, error);
    	content = nbr_stub2.at(i) / nbr_event2.at(i);
    	h2->SetBinContent(i, content);
    	content = nbr_stub3.at(i) / nbr_event3.at(i);
    	h3->SetBinContent(i, content);
    }
    
    h1->SetLineColor(kGreen+2);
    h1->SetLineWidth(1);
    //h1->SetLineStyle(9);
    h1->Draw("SAME");

    h2->SetLineColor(kBlue+2);
    h2->SetLineWidth(1);
    //h1->SetLineStyle(9);
    h2->Draw("SAME");

    h3->SetLineColor(kMagenta+2);
    h3->SetLineWidth(1);
    //h1->SetLineStyle(9);
    h3->Draw("SAME");

    TGraphErrors *gr1 = new TGraphErrors(n,x1,y1,ex1,ey1);
    gr1->SetName("gr1");
    gr1->SetMarkerColor(kRed+2);
    gr1->SetMarkerStyle(20);
    gr1->SetMarkerSize(1.2);
    gr1->Draw("SAME P");

    TF1* func1 = new TF1("func1", "([0]/(1+ TMath::Exp(-[1]*(x-[2]))))", x1[0], x1[9]);
    func1->SetParameters(1, 15, 2);
    func1->SetLineColor(kRed+2);
    func1->SetLineWidth(1);
    func1->SetLineStyle(7);
    gr1->Fit(func1);

    TAxis *xaxis = h1->GetXaxis();
    TAxis *yaxis = h1->GetYaxis();
    xaxis->SetTitle("Emulated pT [GeV]");
    //xaxis->Set(25, 1.0, 3.5);
    xaxis->SetRangeUser(1.0, 3.5);

    yaxis->SetTitle("Stub efficiency");
    yaxis->SetRangeUser(0., 1.1);

    c1->RedrawAxis();

    auto legend = new TLegend(0.65,0.1,0.9,0.35);
    legend->AddEntry("gr1","Adam et al. 2020","p");
    legend->AddEntry("h1","Geant4: max cluster width = 1","l");
    legend->AddEntry("h2","Geant4: max cluster width = 2","l");
    legend->AddEntry("h3","Geant4: max cluster width = 3","l");
    legend->Draw();

    gPad->Modified();
    //*********************** 

    // Close file when finished
    //f.Close();
}