#include "TFile.h"
#include "TTree.h"
#include <iostream>

std::default_random_engine generator;
std::uniform_real_distribution<> distribution1(0.0, 1.0);

bool CBC2(
    const vector<double> &strip_A, // vector for sensor A strips' data
    const vector<double> &strip_B, // vector for sensor B strips' data
    vector<double> &res_A, // [0:4]-> 1-5 strip wide clusters,[5]-> number of clusters (tot), [6]-> mean cluster width (tot), [7:8]-> accepted
    vector<double> &res_B,
    const int MAX_CLUSTER_WIDTH = 3,
    const int CLUSTER_WINDOW = 5,
    const double THRESHOLD = 5.1975 // (119.86 - 106) * 375 e- = 5197.5 e- = 5.1975 ke-
    ){

    const int NBR_STRIP = strip_A.size(); // get number of strips

    double noise = 0.;
    // Clusters in sensor A
    std::vector<double> clus_pos_A; // create vector for clusters position
    std::vector<double> clus_size_A; // create vector for clusters size
    bool inside = false; // is inside a cluster?
    int size = 0; // size of the cluster we are in
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
        double strip_energy = (strip_A[i] / 0.00362) + abs(dist1(generator)); // change MeV in ke-, and apply noise
        if (strip_energy < THRESHOLD && !inside)       
        {} else if (strip_energy < THRESHOLD && inside)
        {
            clus_size_A.push_back(size);
            if(size <= 5) res_A.at(size - 1) += 1; // fill stats for cluster size
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


void figure20_5() {
	const Int_t n = 24; // data for Adam2020
 
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

    const Int_t m = 16; // // Adam2020 irradiated
 
    Double_t x2[m] = {  1.57853,
                        1.71422,
                        1.87568,
                        1.93064,
                        1.98904,
                        2.07149,
                        2.13847,
                        2.1608,
                        2.18485,
                        2.31367,
                        2.42703,
                        2.51978,
                        2.80491,
                        2.84613,
                        3.1656,
                        3.92307};
    Double_t y2[m] = {  0.001628,
                        0.001628,
                        0.007057,
                        0.024699,
                        0.097981,
                        0.305615,
                        0.486106,
                        0.547175,
                        0.598744,
                        0.876945,
                        0.947514,
                        0.955656,
                        0.963799,
                        0.96787,
                        0.957013,
                        0.965156};
    Double_t ex2[m] = {0.};
    Double_t ey2[m] = {0.};

    //************* variable ***************//
    const int NBR_STRIP = 254;

    const int MAX_CLUSTER_WIDTH = 3;
    const int CLUSTER_WINDOW =5;
    const double THRESHOLD = 5.1975; // MeV -> = 3 * (1000 * 3.62 keV)

    const int NBR_BINS = 200;

    Double_t momentum;

    //****************** Create Histo ************************************//
    auto c1 = new TCanvas("c1","c1",1000,600);
    c1->SetTitle("");
    gStyle->SetOptStat(0);
    gPad->SetGridx(1);
    gPad->SetGridy(1);

    TH1* h1 = new TH1D("h1", "", NBR_BINS, 1.0, 3.5);
    h1->SetName("h1");
    TH1* h2 = new TH1D("h2", "", NBR_BINS, 1.0, 3.5);
    h2->SetName("h2");
    TH1* h3 = new TH1D("h3", "", NBR_BINS, 1.0, 3.5);
    h3->SetName("h3");

    gPad->SetTitle("");

    // Create vectors for bin content
    std::vector<double> nbr_stub1(NBR_BINS+2, 0);
    std::vector<double> nbr_event1(NBR_BINS+2, 0);
    std::vector<double> nbr_stub2(NBR_BINS+2, 0);
    std::vector<double> nbr_event2(NBR_BINS+2, 0);
    std::vector<double> nbr_stub3(NBR_BINS+2, 0);
    std::vector<double> nbr_event3(NBR_BINS+2, 0);

    /********************************************************************************************/
    // open file
    TFile *f = TFile::Open("/media/matthieu/ssd1/Geant4/Data/Data_figure20-1/data1M.root", "read");
    auto data = f->Get<TTree>("data");
    const int ENTRIES = data->GetEntries();
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
    data->SetBranchAddress("momentum", &momentum);
    //****************** Main loop over all entries **********************//
    int count_loop = 0;
    double percentage = 0.;
    for (int k = 0; k < ENTRIES; k++)
    {
        data->GetEntry(k);
        std::vector<double> res_A(9, 0);
        std::vector<double> res_B(9, 0);

        bool stub = CBC2(strip_A, strip_B, res_A, res_B, MAX_CLUSTER_WIDTH, CLUSTER_WINDOW, THRESHOLD);

        if(stub){nbr_stub1.at(h1->FindBin(momentum)) += 1.;}
        nbr_event1.at(h1->FindBin(momentum)) += 1.;

        count_loop += 1;
        if (count_loop == ENTRIES / 100)
        {
        	count_loop = 0;
        	percentage += 0.333333333;
        	cout << percentage << " %" << endl;
        }
    }
    double error = 0., content = 0.;
    for (int i = 1; i <= NBR_BINS; ++i)
    {
    	content = nbr_stub1.at(i) / nbr_event1.at(i);
    	h1->SetBinContent(i, content);
    }
    f->Close();
    /**********************************************************************************************/
    // open file
    TFile *f2 = TFile::Open("/media/matthieu/ssd1/Geant4/Data/Data_figure20-1/data1M-260.root", "read");
    auto data2 = f2->Get<TTree>("data");
    const int ENTRIES2 = data2->GetEntries();
    std::vector<double> strip_A2(NBR_STRIP, 0);
    std::vector<double> strip_B2(NBR_STRIP, 0);
    for (int i = 0; i < NBR_STRIP; ++i)
    {
        string strip_name = "s" + std::to_string(i);
        char const *pchar = strip_name.c_str();
        data2->SetBranchAddress(pchar, &strip_B2.at(i));
    }
    for (int i = NBR_STRIP; i < (2*NBR_STRIP); ++i)
    {
        string strip_name = "s" + std::to_string(i);
        char const *pchar = strip_name.c_str();
        data2->SetBranchAddress(pchar, &strip_A2.at(i-NBR_STRIP));
    }
    data2->SetBranchAddress("momentum", &momentum);
    //****************** Main loop over all entries **********************//
    count_loop = 0;
    for (int k = 0; k < ENTRIES2; k++)
    {
        data2->GetEntry(k);
        std::vector<double> res_A2(9, 0);
        std::vector<double> res_B2(9, 0);

        bool stub = CBC2(strip_A2, strip_B2, res_A2, res_B2, MAX_CLUSTER_WIDTH, CLUSTER_WINDOW, THRESHOLD);

        if(stub){nbr_stub2.at(h2->FindBin(momentum)) += 1.;}
        nbr_event2.at(h2->FindBin(momentum)) += 1.;

        count_loop += 1;
        if (count_loop == ENTRIES2 / 100)
        {
            count_loop = 0;
            percentage += 0.333333333;
            cout << percentage << " %" << endl;
        }
    }
    error = 0.;
    content = 0.;
    for (int i = 1; i <= NBR_BINS; ++i)
    {
        content = nbr_stub2.at(i) / nbr_event2.at(i);
        h2->SetBinContent(i, content);
    }
    f2->Close();
    /************************************************************************************************************/
    // open file
    TFile *f3 = TFile::Open("/media/matthieu/ssd1/Geant4/Data/Data_figure20-1/data1M-250.root", "read");
    auto data3 = f3->Get<TTree>("data");
    const int ENTRIES3 = data3->GetEntries();
    std::vector<double> strip_A3(NBR_STRIP, 0);
    std::vector<double> strip_B3(NBR_STRIP, 0);
    for (int i = 0; i < NBR_STRIP; ++i)
    {
        string strip_name = "s" + std::to_string(i);
        char const *pchar = strip_name.c_str();
        data3->SetBranchAddress(pchar, &strip_B3.at(i));
    }
    for (int i = NBR_STRIP; i < (2*NBR_STRIP); ++i)
    {
        string strip_name = "s" + std::to_string(i);
        char const *pchar = strip_name.c_str();
        data3->SetBranchAddress(pchar, &strip_A3.at(i-NBR_STRIP));
    }
    data3->SetBranchAddress("momentum", &momentum);
    //****************** Main loop over all entries **********************//
    count_loop = 0;
    for (int k = 0; k < ENTRIES3; k++)
    {
        data3->GetEntry(k);
        std::vector<double> res_A3(9, 0);
        std::vector<double> res_B3(9, 0);

        bool stub = CBC2(strip_A3, strip_B3, res_A3, res_B3, MAX_CLUSTER_WIDTH, CLUSTER_WINDOW, THRESHOLD);

        if(stub){nbr_stub3.at(h3->FindBin(momentum)) += 1.;}
        nbr_event3.at(h3->FindBin(momentum)) += 1.;

        count_loop += 1;
        if (count_loop == ENTRIES3 / 100)
        {
            count_loop = 0;
            percentage += 0.333333333;
            cout << percentage << " %" << endl;
        }
    }
    error = 0.;
    content = 0.;
    for (int i = 1; i <= NBR_BINS; ++i)
    {
        content = nbr_stub3.at(i) / nbr_event3.at(i);
        h3->SetBinContent(i, content);
    }
    f3->Close();
    /*******************************************************************************************************/
    
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

    TGraphErrors *gr1 = new TGraphErrors(n,x1,y1,ex1,ey1); // non-irradiated
    gr1->SetName("gr1");
    gr1->SetMarkerColor(12);
    gr1->SetMarkerStyle(24);
    gr1->SetMarkerSize(0.7);
    gr1->Draw("SAME P");

    TGraphErrors *gr2 = new TGraphErrors(m,x2,y2,ex2,ey2); // irradiated
    gr2->SetName("gr2");
    gr2->SetMarkerColor(12);
    gr2->SetMarkerStyle(25);
    gr2->SetMarkerSize(0.7);
    gr2->Draw("SAME P");

    TF1* func1 = new TF1("func1", "([0]/(1+ TMath::Exp(-[1]*(x-[2]))))", x1[0], x1[n-1]);
    func1->SetParameters(1, 15, 2);
    func1->SetLineColor(13);
    func1->SetLineWidth(1);
    func1->SetLineStyle(7);
    gr1->Fit(func1);

    TF1* func2 = new TF1("func2", "([0]/(1+ TMath::Exp(-[1]*(x-[2]))))", x2[0], x2[m-1]);
    func2->SetParameters(1, 15, 2);
    func2->SetLineColor(13);
    func2->SetLineWidth(1);
    func2->SetLineStyle(7);
    gr2->Fit(func2);

    TAxis *xaxis = h1->GetXaxis();
    TAxis *yaxis = h1->GetYaxis();
    xaxis->SetLabelFont(42);
    xaxis->SetLabelSize(0.04);
    xaxis->SetTitle("pT [GeV]");
    xaxis->SetTitleFont(22);
    xaxis->SetTitleSize(0.05);
    xaxis->SetTitleOffset(0.95);
    xaxis->SetRangeUser(1.0, 3.5);
    yaxis->SetLabelFont(42);
    yaxis->SetLabelSize(0.04);
    yaxis->SetTitle("Stub efficiency");
    yaxis->SetTitleFont(22);
    yaxis->SetTitleSize(0.05);
    yaxis->SetTitleOffset(0.9);
    yaxis->SetRangeUser(0., 1.05);

    TF1* f1 = new TF1("f1", "x", 1.0, 3.5); // 3.5
    TGaxis* A1 = new TGaxis(1.0, 1.05, 3.5, 1.05, "f1", 510, "-");
    A1->SetLabelFont(42);
    A1->SetLabelSize(0.04);
    A1->SetTitle("Angle d'incidence [deg]");
    A1->SetTitleFont(22);
    A1->SetTitleSize(0.05);
    A1->SetTitleOffset(0.95);
    A1->ChangeLabel(1, -1, -1, -1, -1, -1, "19.6");
    A1->ChangeLabel(2, -1, -1, -1, -1, -1, "13.1");
    A1->ChangeLabel(3, -1, -1, -1, -1, -1, "9.8");
    A1->ChangeLabel(4, -1, -1, -1, -1, -1, "7.8");
    A1->ChangeLabel(5, -1, -1, -1, -1, -1, "6.5");
    A1->ChangeLabel(6, -1, -1, -1, -1, -1, "5.6");
    A1->ChangeLabel(7, -1, -1, -1, -1, -1, "");
    A1->Draw("SAME");

    c1->RedrawAxis();

    auto legend = new TLegend(0.65,0.1,0.9,0.35);
    legend->AddEntry("gr1","Adam et al. 2020 - non-irradiated","p");
    legend->AddEntry("gr2","Adam et al. 2020 - irradiated","p");
    legend->AddEntry("h1","Geant4: W = 270 um","l");
    legend->AddEntry("h2","Geant4: W = 260 um","l");
    legend->AddEntry("h3","Geant4: W = 250 um","l");
    legend->Draw();

    gPad->Modified();
    //*********************** 
    c1->SaveAs("figure20_5.pdf");
    // Close file when finished
    //f.Close();
}