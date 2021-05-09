#include "TFile.h"
#include "TTree.h"
#include <iostream>

std::default_random_engine generator;
std::normal_distribution<double> distribution(0.0, 0.795); // µ = 1000 e-, s = 800 e-

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
        double strip_energy = (strip_A[i] / 0.00362) + distribution(generator);
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
        double strip_energy = (strip_B[i] / 0.00362) + distribution(generator);
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
	Double_t y[],
	Double_t ey[],
	Double_t p1,
	Double_t p2,
	Double_t p3,
	Double_t t1,
	Double_t t2,
	Double_t t3,
	Double_t s1,
	Double_t s2,
	Double_t s3
	){

	// open file
	ofstream myfile;
    myfile.open ("figure20_0_data.txt");
    myfile << "x\ty\tey\n";
    for (int i = 0; i < k; ++i)
    {
    	myfile << std::scientific << x[i] << "\t" << y[i] << "\t" << ey[i] << "\n";
    }
    myfile << std::scientific << "serie\t" << "param0\t" << "param1\t" << "param2\n";
    myfile << std::scientific << "Adam\t" << p1 << "\t" << t1 << "\t" << s1 << "\n";
    myfile << std::scientific << "Adam-irr\t" << p2 << "\t" << t2 << "\t" << s2 << "\n";
    myfile << std::scientific << "Geant4\t" << p3 << "\t" << t3 << "\t" << s3 << "\n";
    myfile.close();

    return;
}

void figure20_0() {
	// data for Adam2020
	const Int_t n = 24; // Adam2020 non-irradiated
 
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

    const Int_t p = 41; // // Simulation
 
    Double_t x3[p] = {  1.20000,
                        1.30000,
                        1.32089,
                        1.40000,
                        1.48578,
                        1.50000,
                        1.57853,
                        1.60000,
                        1.70000,
                        1.71422,
                        1.75000,
                        1.79152,
                        1.80000,
                        1.85000,
                        1.87568,
                        1.90000,
                        1.93064,
                        1.95000,
                        1.98904,
                        2.00000,
                        2.05000,
                        2.10000,
                        2.13847,
                        2.15000,
                        2.20000,
                        2.25000,
                        2.30000,
                        2.40000,
                        2.42703,
                        2.50000,
                        2.51978,
                        2.60000,
                        2.70000,
                        2.80000,
                        2.80491,
                        2.90000,
                        2.97495,
                        3.00000,
                        3.10000,
                        3.16560,
                        3.20000};
    Double_t y3[p] = {0.};
    Double_t ex3[p] = {0.};
    Double_t ey3[p] = {0.};

	//****************** Create Histo ************************************//
    auto c1 = new TCanvas("c1","c1",1000,600);
    c1->SetTitle("Cluster efficiency for 2S mini-module");
    gStyle->SetOptStat(0);
    gPad->SetGridx(1);
    gPad->SetGridy(1);

    gPad->SetTitle("Cluster efficiency for 2S mini-module");

    // loop on all files
    for (int j = 0; j < p; ++j)
    {
        //CopyFile(j);

        string file_path = "/media/matthieu/ssd1/Geant4/Data/Data_figure20/data_";
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

        //****************** Main loop over all entries **********************//
        int count_loop = 0;
        bool stop = true;
        double nbr_stub = 0.;
        std::vector<double> efficiencies;
        for (int k = 0; k < ENTRIES; k++)
        {
            // fill variables with datas from entry i
            data->GetEntry(k);

            std::vector<double> res_A(9, 0);
            std::vector<double> res_B(9, 0);

            bool stub = CBC2(strip_A, strip_B, res_A, res_B, MAX_CLUSTER_WIDTH, CLUSTER_WINDOW, THRESHOLD);

            if(stub){nbr_stub += 1.;}

            count_loop += 1;
            if (count_loop == ENTRIES /10)
            {
                count_loop = 0;
                efficiencies.push_back((double) nbr_stub / (ENTRIES / 10));
                //cout << mean_cluster_width << endl;
                nbr_stub = 0.;
            }
        }
        //********************* fig 17 computation and printing **********************************//
        // for A
        double variance = 0., deviation = 0., average = 0.;
        average = std::accumulate(efficiencies.begin(), efficiencies.end(), 0.0) / 10;
        for (int i = 0; i < 10; ++i){variance += pow((efficiencies.at(i) - average), 2);}
        variance /= 9;
        deviation = sqrt(variance);
        cout << std::scientific << "File " << j << " stub efficiency:\t" << average << "\t+/- " << deviation << endl;

        //*********************** 
        y3[j] = average;
        ey3[j] = deviation;

        // Close file when finished
        f.Close();
    }   
    // Fill graphs
    TGraphErrors *gr1 = new TGraphErrors(n,x1,y1,ex1,ey1); // non-irradiated
    gr1->SetName("gr1");
    gr1->SetMarkerColor(13);
    gr1->SetMarkerStyle(24);
    gr1->SetMarkerSize(0.7);

    TGraphErrors *gr2 = new TGraphErrors(m,x2,y2,ex2,ey2); // irradiated
    gr2->SetName("gr2");
    gr2->SetMarkerColor(13);
    gr2->SetMarkerStyle(25);
    gr2->SetMarkerSize(0.7);

    TGraphErrors *gr3 = new TGraphErrors(p,x3,y3,ex3,ey3); // simulation
    gr3->SetName("gr3");
    gr3->SetMarkerColor(kRed+2);
    gr3->SetMarkerStyle(20);
    gr3->SetMarkerSize(0.7);

    TF1* func1 = new TF1("func1", "(0.5*[0]*(1+ TMath::Erf((x-[1])/[2])))", x1[0], x1[n-1]);
    func1->SetParameters(0.98, 1.85, 0.10);
    func1->SetLineColor(13);
    func1->SetLineWidth(1);
    func1->SetLineStyle(7);
    TFitResultPtr r1 = gr1->Fit(func1, "S");

    TF1* func2 = new TF1("func2", "(0.5*[0]*(1+ TMath::Erf((x-[1])/[2])))", x2[0], x2[m-1]);
    func2->SetParameters(0.95, 2.15, 0.10);
    func2->SetLineColor(13);
    func2->SetLineWidth(1);
    func2->SetLineStyle(7);
    TFitResultPtr r2 = gr2->Fit(func2, "S");

    TF1* func3 = new TF1("func3", "(0.5*[0]*(1+ TMath::Erf((x-[1])/[2])))", x3[0], x3[p-1]);
    func3->SetParameters(0.99, 1.85, 0.10);
    func3->SetLineColor(kRed+2);
    func3->SetLineWidth(1);
    func3->SetLineStyle(7);
    TFitResultPtr r3 = gr3->Fit(func3, "S");

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(gr1);
    mg->Add(gr2);
    mg->Add(gr3);
    mg->SetTitle("");
    mg->Draw("AP");

    Double_t plateau1   = r1->Value(0);
    Double_t turn1   = r1->Value(1);
    Double_t sigma1   = r1->Value(2);

    Double_t plateau2   = r2->Value(0);
    Double_t turn2   = r2->Value(1);
    Double_t sigma2   = r2->Value(2);

    Double_t plateau3   = r3->Value(0);
    Double_t turn3   = r3->Value(1);
    Double_t sigma3   = r3->Value(2);

    TAxis *xaxis = mg->GetXaxis();
    TAxis *yaxis = mg->GetYaxis();
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

    TF1* f1 = new TF1("f1", "x", 1.0, 3.5);
    TGaxis* A1 = new TGaxis(1.0, yaxis->GetXmax(), 3.5, yaxis->GetXmax(), "f1", 510, "-");
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

    auto legend = new TLegend(0.6,0.1,0.9,0.4);
    legend->AddEntry("gr1","Adam et al. - non-irr.","p");
    legend->AddEntry("gr2","Adam et al. - irr.","p");
    legend->AddEntry("gr3","Geant4 - 2.63 mm","ep");
    legend->Draw();

    gPad->Modified();

    c1->SaveAs("figure20_0.pdf");

    SaveData(p, x3, y3, ey3, plateau1, plateau2, plateau3, turn1, turn2, turn3, sigma1, sigma2, sigma3);
}