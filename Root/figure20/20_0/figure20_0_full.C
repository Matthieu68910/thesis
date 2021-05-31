#include "TFile.h"
#include "TTree.h"
#include <iostream>

std::default_random_engine generator;
std::uniform_real_distribution<> distribution1(0.0, 1.0);

bool CBC3(
    const vector<double> &strip_A, 
    const vector<double> &strip_B, 
    vector<double> &res,                    // global results for sensors (stubs, bend)
                                                // [0]  -> nbr stubs
                                                // [1]  -> mean bend information
    vector<double> &res_A,                  // results for sensor A
                                                // [0:4]-> 1-5 strip wide clusters
                                                // [5]  -> number of clusters (tot)
                                                // [6]  -> mean cluster width (tot)
                                                // [7:8]-> idem accepted
                                                // [9]  -> nbr hits
    vector<double> &res_B,                  // results for sensor B
                                                // [0:4]-> 1-5 strip wide clusters
                                                // [5]  -> number of clusters (tot)
                                                // [6]  -> mean cluster width (tot)
                                                // [7:8]-> idem accepted
                                                // [9]  -> nbr hits
    const int MAX_CLUSTER_WIDTH = 4,        // max 4
    const double CLUSTER_WINDOW = 7.,       // max +- 7 en 1/2 strip
    const double THRESHOLD = 0,             // default to 0 (in killo-electrons)   
    double NOISE = 1.,                      // noise in killo-elctrons
    double PAIR_CREATION_ENERGY = 0.00362   // pair creation energy in keV            
    ){

    // find number of strips
    const int NBR_STRIP = strip_A.size();

    // noise distribution deternmination
    std::normal_distribution<double> dist1(0., NOISE);

    // Clusters in sensor A
    std::vector<double> clus_pos_A;
    std::vector<double> clus_size_A;
    bool inside = false;
    int size = 0;
    int nbr_hits_A = 0;

    // Loop on sensor A strips
    for (int i = 0; i < NBR_STRIP; ++i)
    {
        // convert strip signal to killo-electrons and add noise 
        double strip_energy = (strip_A[i] / PAIR_CREATION_ENERGY) + abs(dist1(generator));

        // cluster finding logic
        if (strip_energy < THRESHOLD && !inside)        
        {
            // do nothing
        } 
        else if (strip_energy < THRESHOLD && inside)
        {
            clus_size_A.push_back(size);
            if(size <= 5) res_A.at(size - 1) += 1;
            clus_pos_A.push_back((double) (i - 1) - ((double) size / 2) + 0.5);
            size = 0;
            inside = false;
        } 
        else if (strip_energy >= THRESHOLD && !inside)
        {
            nbr_hits_A += 1;
            size = 1;
            inside = true;
            if (i == (NBR_STRIP - 1))
            {
                clus_size_A.push_back(1);
                res_A.at(0) += 1;
                clus_pos_A.push_back((double) i);
                size = 0;
                inside = false;
            }
        } 
        else if (strip_energy >= THRESHOLD && inside)
        {
            nbr_hits_A += 1;
            size += 1;
            if (i == (NBR_STRIP - 1))
            {
                clus_size_A.push_back(size);
                if(size <= 5) res_A.at(size - 1) += 1;
                clus_pos_A.push_back((double) i - ((double) size / 2) + 0.5);
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
        // convert strip signal to killo-electrons and add noise 
        double strip_energy = (strip_B[i] / PAIR_CREATION_ENERGY) + abs(dist1(generator));

        // cluster finding logic
        if (strip_energy < THRESHOLD && !inside)        
        {
            // do nothing
        } 
        else if (strip_energy < THRESHOLD && inside)
        {
            clus_size_B.push_back(size);
            if(size <= 5) res_B.at(size - 1) += 1;
            clus_pos_B.push_back((double) (i - 1) - ((double) size / 2) + 0.5);
            size = 0;
            inside = false;
        } 
        else if (strip_energy >= THRESHOLD && !inside)
        {
            nbr_hits_B += 1;
            size = 1;
            inside = true;
            if (i == (NBR_STRIP - 1))
            {
                clus_size_B.push_back(1);
                res_B.at(0) += 1;
                clus_pos_B.push_back((double) i);
                size = 0;
                inside = false;
            }
        } 
        else if (strip_energy >= THRESHOLD && inside)
        {
            nbr_hits_B += 1;
            size += 1;
            if (i == (NBR_STRIP - 1))
            {
                clus_size_B.push_back(size);
                if(size <= 5) res_B.at(size - 1) += 1;
                clus_pos_B.push_back((double) i - ((double) size / 2) + 0.5);
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
        res.at(0) = 0.;
        res.at(1) = 0.;
        return false;
    }

    // stub finding logic
    std::vector<double> stub_pos;
    std::vector<double> stub_bend;

    // loop over all identified clusters
    for (int i = 0; i < clus_pos_A.size(); ++i) // for each seed in sensor A ....
    {
        for (int j = 0; j < clus_pos_B.size(); ++j) // .... look for correlation in sensor B.
        {
            if (abs(clus_pos_A.at(i) - clus_pos_B.at(j)) <= CLUSTER_WINDOW)  // correlation only if inside window
            {
                stub_pos.push_back(clus_pos_A.at(i));
                stub_bend.push_back(abs((double) clus_pos_A.at(i) - clus_pos_B.at(j)));
            }
        }
    }
    if (stub_pos.size() != 0)   
    {
        res.at(0) = stub_pos.size();
        res.at(1) = std::accumulate(stub_bend.begin(), stub_bend.end(), 0.0) / stub_pos.size();
        return true;
    } 
    else
    {
        res.at(0) = 0.;
        res.at(1) = 0.;
        return false;
    }
}

bool CBC2(
    const vector<double> &strip_A, // vector for sensor A strips' data
    const vector<double> &strip_B, // vector for sensor B strips' data
    vector<double> &res_A, // [0:4]-> 1-5 strip wide clusters,[5]-> number of clusters (tot), [6]-> mean cluster width (tot), [7:8]-> accepted
    vector<double> &res_B,
    const int MAX_CLUSTER_WIDTH = 3,
    const int CLUSTER_WINDOW = 5,
    const double THRESHOLD = 5.1975, // (119.86 - 106) * 375 e- = 5197.5 e- = 5.1975 ke-             // default to 0 (in killo-electrons)   
    double NOISE = 1.,                      // noise in killo-elctrons
    double PAIR_CREATION_ENERGY = 0.00362   // pair creation energy in keV
    ){

    const int NBR_STRIP = strip_A.size(); // get number of strips

    // noise distribution deternmination
    std::normal_distribution<double> dist1(0., NOISE);

    // Clusters in sensor A
    std::vector<double> clus_pos_A; // create vector for clusters position
    std::vector<double> clus_size_A; // create vector for clusters size
    bool inside = false; // is inside a cluster?
    int size = 0; // size of the cluster we are in
    // Loop on sensor A strips
    for (int i = 0; i < NBR_STRIP; ++i)
    {
        // convert strip signal to killo-electrons and add noise 
        double strip_energy = (strip_A[i] / PAIR_CREATION_ENERGY) + abs(dist1(generator));

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
        // convert strip signal to killo-electrons and add noise 
        double strip_energy = (strip_B[i] / PAIR_CREATION_ENERGY) + abs(dist1(generator));

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
	Double_t x1[],
	Double_t y1[],
	Double_t ey1[],
    const int &l,
    Double_t x2[],
    Double_t y2[],
    Double_t ey2[],
	Double_t p1,
	Double_t p2,
	Double_t p3,
    Double_t p4,
	Double_t t1,
	Double_t t2,
	Double_t t3,
    Double_t t4,
	Double_t s1,
	Double_t s2,
	Double_t s3,
    Double_t s4
	){

	// open file
	ofstream myfile;
    myfile.open ("figure20_0-full-txt");
    myfile << "x1\ty1\tey1\tx2\ty2\tey2\n";
    for (int i = 0; i < k; ++i)
    {
    	myfile << std::scientific << x1[i] << "\t" << y1[i] << "\t" << ey1[i] << "\t" << x2[i] << "\t" << y2[i] << "\t" << ey2[i] << "\n";
    }
    myfile << std::scientific << "serie\t" << "param0\t" << "param1\t" << "param2\n";
    myfile << std::scientific << "Adam\t" << p1 << "\t" << t1 << "\t" << s1 << "\n";
    myfile << std::scientific << "Adam-irr\t" << p2 << "\t" << t2 << "\t" << s2 << "\n";
    myfile << std::scientific << "Grand\t" << p3 << "\t" << t3 << "\t" << s3 << "\n";
    myfile << std::scientific << "Mini\t" << p4 << "\t" << t4 << "\t" << s4 << "\n";
    myfile.close();

    return;
}

void figure20_0_full() {
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

    const Int_t p = 25; // // Grand module
 
    Double_t x3[p] = {  1.32089E+00,
                        1.48578E+00,
                        1.57853E+00,
                        1.71422E+00,
                        1.79152E+00,
                        1.87568E+00,
                        1.93064E+00,
                        1.98904E+00,
                        2.13847E+00,
                        2.42703E+00,
                        2.51978E+00,
                        2.80491E+00,
                        2.97495E+00,
                        3.16560E+00,
                        4.26316E+00,
                        4.66851E+00,
                        5.15975E+00,
                        6.32257E+00,
                        6.32257E+00,
                        7.25867E+00,
                        9.79902E+00,
                        1.08880E+01,
                        1.22483E+01,
                        1.30642E+01,
                        1.50738E+01};
    Double_t y3[p] = {  8.2200E-03,
                        8.6300E-03,
                        9.6500E-03,
                        2.1400E-02,
                        1.2312E-01,
                        4.6088E-01,
                        6.6087E-01,
                        7.8988E-01,
                        9.7566E-01,
                        9.8778E-01,
                        9.8829E-01,
                        9.8960E-01,
                        9.8984E-01,
                        9.9010E-01,
                        9.9007E-01,
                        9.9031E-01,
                        9.9056E-01,
                        9.9079E-01,
                        9.9071E-01,
                        9.9063E-01,
                        9.9120E-01,
                        9.9151E-01,
                        9.9088E-01,
                        9.9093E-01,
                        9.9158E-01};
    Double_t ex3[p] = { 0.};
    Double_t ey3[p] = { 2.99353E-03,
                        2.84855E-03,
                        3.04636E-03,
                        4.55050E-03,
                        8.91178E-03,
                        1.60014E-02,
                        1.47557E-02,
                        1.35246E-02,
                        4.99943E-03,
                        3.52360E-03,
                        3.57120E-03,
                        3.68179E-03,
                        3.09029E-03,
                        3.50613E-03,
                        2.69775E-03,
                        3.20572E-03,
                        3.09878E-03,
                        2.99594E-03,
                        3.13112E-03,
                        3.14195E-03,
                        3.12694E-03,
                        2.86919E-03,
                        2.69035E-03,
                        3.02600E-03,
                        2.70868E-03};

    const Int_t q = 25; // mini module
 
    Double_t x4[q] = {  1.32089E+00,
                        1.48578E+00,
                        1.57853E+00,
                        1.71422E+00,
                        1.79152E+00,
                        1.87568E+00,
                        1.93064E+00,
                        1.98904E+00,
                        2.13847E+00,
                        2.42703E+00,
                        2.51978E+00,
                        2.80491E+00,
                        2.97495E+00,
                        3.16560E+00,
                        4.26316E+00,
                        4.66851E+00,
                        5.15975E+00,
                        6.32257E+00,
                        6.32257E+00,
                        7.25867E+00,
                        9.79902E+00,
                        1.08880E+01,
                        1.22483E+01,
                        1.30642E+01,
                        1.50738E+01};
    Double_t y4[q] = {  4.10000E-03,
                        3.73000E-03,
                        5.01000E-03,
                        3.07800E-02,
                        2.60730E-01,
                        5.22670E-01,
                        6.76430E-01,
                        8.30190E-01,
                        9.86340E-01,
                        9.92610E-01,
                        9.92440E-01,
                        9.92670E-01,
                        9.92880E-01,
                        9.92770E-01,
                        9.92810E-01,
                        9.93070E-01,
                        9.93180E-01,
                        9.93670E-01,
                        9.92690E-01,
                        9.93140E-01,
                        9.93230E-01,
                        9.93730E-01,
                        9.93790E-01,
                        9.93160E-01,
                        9.93620E-01};
    Double_t ex4[q] = { 0.};
    Double_t ey4[q] = { 1.83952E-03,
                        1.94290E-03,
                        2.21790E-03,
                        5.87338E-03,
                        1.33590E-02,
                        1.59659E-02,
                        1.45132E-02,
                        1.16997E-02,
                        3.87225E-03,
                        2.60107E-03,
                        3.05611E-03,
                        2.55468E-03,
                        2.40067E-03,
                        2.53004E-03,
                        2.51739E-03,
                        2.73125E-03,
                        2.50002E-03,
                        2.48655E-03,
                        2.31680E-03,
                        2.80339E-03,
                        2.46124E-03,
                        2.26459E-03,
                        2.41333E-03,
                        2.43178E-03,
                        2.62767E-03};

	//****************** Create Histo ************************************//
    auto c1 = new TCanvas("c1","c1",1000,600);
    c1->SetTitle("Cluster efficiency for 2S mini-module");
    gStyle->SetOptStat(0);
    gPad->SetGridx(1);
    gPad->SetGridy(1);

    gPad->SetTitle("Cluster efficiency for 2S mini-module");

    // Fill graphs
    TGraphErrors *gr1 = new TGraphErrors(n,x1,y1,ex1,ey1); // non-irradiated
    gr1->SetName("gr1");
    gr1->SetMarkerColor(12);
    gr1->SetMarkerStyle(24);
    gr1->SetMarkerSize(0.7);

    TGraphErrors *gr2 = new TGraphErrors(m,x2,y2,ex2,ey2); // irradiated
    gr2->SetName("gr2");
    gr2->SetMarkerColor(12);
    gr2->SetMarkerStyle(25);
    gr2->SetMarkerSize(0.7);

    TGraphErrors *gr3 = new TGraphErrors(p,x3,y3,ex3,ey3); // Maxi
    gr3->SetName("gr3");
    gr3->SetMarkerColor(kBlue+2);
    gr3->SetMarkerStyle(20);
    gr3->SetMarkerSize(0.7);

    TGraphErrors *gr4 = new TGraphErrors(q,x4,y4,ex4,ey4); // Mini
    gr4->SetName("gr4");
    gr4->SetMarkerColor(kRed+2);
    gr4->SetMarkerStyle(20);
    gr4->SetMarkerSize(0.7);

    TF1* func1 = new TF1("func1", "(0.5*[0]*(1+ TMath::Erf((x-[1])/[2])))", x1[0], x1[n-1]);
    func1->SetParameters(0.98, 1.85, 0.10);
    func1->SetLineColor(12);
    func1->SetLineWidth(1);
    func1->SetLineStyle(7);
    TFitResultPtr r1 = gr1->Fit(func1, "S");

    TF1* func2 = new TF1("func2", "(0.5*[0]*(1+ TMath::Erf((x-[1])/[2])))", x2[0], x2[m-1]);
    func2->SetParameters(0.95, 2.15, 0.10);
    func2->SetLineColor(12);
    func2->SetLineWidth(1);
    func2->SetLineStyle(7);
    TFitResultPtr r2 = gr2->Fit(func2, "S");

    TF1* func3 = new TF1("func3", "(0.5*[0]*(1+ TMath::Erf((x-[1])/[2])))", x3[0], x3[p-1]);
    func3->SetParameters(0.99, 1.85, 0.10);
    func3->SetLineColor(kBlue+2);
    func3->SetLineWidth(1);
    func3->SetLineStyle(7);
    TFitResultPtr r3 = gr3->Fit(func3, "S");

    TF1* func4 = new TF1("func4", "(0.5*[0]*(1+ TMath::Erf((x-[1])/[2])))", x4[0], x4[p-1]);
    func4->SetParameters(0.99, 1.85, 0.10);
    func4->SetLineColor(kRed+2);
    func4->SetLineWidth(1);
    func4->SetLineStyle(7);
    TFitResultPtr r4 = gr4->Fit(func4, "S");

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(gr1);
    mg->Add(gr2);
    mg->Add(gr4);
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

    Double_t plateau4   = r4->Value(0);
    Double_t turn4   = r4->Value(1);
    Double_t sigma4   = r4->Value(2);

    TAxis *xaxis = mg->GetXaxis();
    TAxis *yaxis = mg->GetYaxis();
    xaxis->SetLabelFont(42);
	xaxis->SetLabelSize(0.04);
    xaxis->SetTitle("pT [GeV]");
    xaxis->SetTitleFont(22);
	xaxis->SetTitleSize(0.05);
	xaxis->SetTitleOffset(0.95);
    xaxis->SetRangeUser(1.3, 2.6);
    yaxis->SetLabelFont(42);
	yaxis->SetLabelSize(0.04);
    yaxis->SetTitle("Stub efficiency");
    yaxis->SetTitleFont(22);
	yaxis->SetTitleSize(0.05);
	yaxis->SetTitleOffset(0.9);

    TF1* f1 = new TF1("f1", "x", 1.3, 2.6); // 3.5
    TGaxis* A1 = new TGaxis(1.3, yaxis->GetXmax(), 2.6, yaxis->GetXmax(), "f1", 510, "-");
    A1->SetLabelFont(42);
	A1->SetLabelSize(0.04);
	A1->SetTitle("Angle d'incidence [deg]");
	A1->SetTitleFont(22);
	A1->SetTitleSize(0.05);
	A1->SetTitleOffset(0.95);
    A1->ChangeLabel(1, -1, -1, -1, -1, -1, "14.1");
    A1->ChangeLabel(2, -1, -1, -1, -1, -1, "12.3");
    A1->ChangeLabel(3, -1, -1, -1, -1, -1, "11.0");
    A1->ChangeLabel(4, -1, -1, -1, -1, -1, "9.8");
    A1->ChangeLabel(5, -1, -1, -1, -1, -1, "8.9");
    A1->ChangeLabel(6, -1, -1, -1, -1, -1, "8.2");
    A1->ChangeLabel(7, -1, -1, -1, -1, -1, "7.6");
	A1->Draw("SAME");

    auto legend = new TLegend(0.1,0.9,0.4,0.6);
    legend->AddEntry("gr1","Adam et al. - mini  module - non-irr.","p");
    legend->AddEntry("gr2","Adam et al. - mini  module - irr.","p");
    legend->AddEntry("gr4","Geant4      - mini  module","ep");
    legend->AddEntry("gr3","Geant4 + By - grand module","ep");
    legend->Draw();

    gPad->Modified();

    c1->SaveAs("figure20_0-full.pdf");

    SaveData(p, x3, y3, ey3, q, x4, y4, ey4, plateau1, plateau2, plateau3, plateau4, turn1, turn2, turn3, turn4, sigma1, sigma2, sigma3, sigma4);
}