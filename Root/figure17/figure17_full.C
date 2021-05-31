#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <random>

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
    const double CLUSTER_WINDOW = 7.,           // max +- 7 en 1/2 strip
    const double THRESHOLD = 0,                   // default to 0 (in killo-electrons)   
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

	double noise = 0.;
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
    	double strip_energy = (strip_A[i] / 0.00362) + abs(dist1(generator));
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
	Double_t ey[]
	){

	// open file
	ofstream myfile;
    myfile.open ("figure17-full_data.txt");
    myfile << "x\ty\tey\n";
    for (int i = 0; i < k; ++i)
    {
    	myfile << std::scientific << x[i] << "\t" << y[i] << "\t" << ey[i] << "\n";
    }
    myfile.close();

    return;
}

void figure17_full() {
	const Int_t n = 25; // mini non-irradiated
 
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
    Double_t y1[n] = {	1.107392E+00,
                        1.117410E+00,
                        1.118199E+00,
                        1.119322E+00,
                        1.124621E+00,
                        1.141700E+00,
                        1.153216E+00,
                        1.154313E+00,
                        1.171031E+00,
                        1.180927E+00,
                        1.193439E+00,
                        1.237726E+00,
                        1.249167E+00,
                        1.263586E+00,
                        1.286432E+00,
                        1.294512E+00,
                        1.325811E+00,
                        1.350611E+00,
                        1.357975E+00,
                        1.368332E+00,
                        1.385643E+00,
                        1.401662E+00,
                        1.433885E+00,
                        1.459057E+00,
                        1.513336E+00};
	Double_t ex1[n] = { 1.255698E-02,
                        1.291455E-02,
                        1.372397E-02,
                        1.209704E-02,
                        1.359613E-02,
                        1.365208E-02,
                        1.356115E-02,
                        1.409360E-02,
                        1.427797E-02,
                        1.238374E-02,
                        1.348783E-02,
                        1.620648E-02,
                        1.584535E-02,
                        1.505500E-02,
                        1.586560E-02,
                        1.559485E-02,
                        1.620876E-02,
                        1.929637E-02,
                        1.850430E-02,
                        1.770109E-02,
                        1.760270E-02,
                        1.874299E-02,
                        1.795925E-02,
                        1.832159E-02,
                        1.960423E-02};
    Double_t ey1[n] = { 1.255698E-02,
                        1.291455E-02,
                        1.372397E-02,
                        1.209704E-02,
                        1.359613E-02,
                        1.365208E-02,
                        1.356115E-02,
                        1.409360E-02,
                        1.427797E-02,
                        1.238374E-02,
                        1.348783E-02,
                        1.620648E-02,
                        1.584535E-02,
                        1.505500E-02,
                        1.586560E-02,
                        1.559485E-02,
                        1.620876E-02,
                        1.929637E-02,
                        1.850430E-02,
                        1.770109E-02,
                        1.760270E-02,
                        1.874299E-02,
                        1.795925E-02,
                        1.832159E-02,
                        1.960423E-02};

    const Int_t m = 35; // adam2020 irradiated
 
	Double_t x3[m] = {	1.50,
						3.50,
						4.00,
						4.50,
						5.00,
						6.20,
						6.60,
						6.60,
						6.90,
						7.00,
						7.20,
						7.50,
						7.80,
						8.00,
						8.10,
						8.50,
						8.90,
						9.00,
						9.10,
						9.10,
						9.20,
						9.50,
						9.80,
						9.90,
						10.00,
						10.20,
						10.30,
						10.50,
						11.20,
						11.50,
						11.70,
						12.00,
						12.50,
						13.00,
						13.50};
    Double_t y3[m] = {	1.3039,
						1.30538,
						1.30784,
						1.30129,
						1.2972,
						1.30071,
						1.30055,
						1.31062,
						1.31523,
						1.30631,
						1.30918,
						1.31317,
						1.30786,
						1.29337,
						1.30754,
						1.31374,
						1.31814,
						1.32376,
						1.31272,
						1.31905,
						1.31604,
						1.3159,
						1.26567,
						1.31861,
						1.32221,
						1.31656,
						1.31801,
						1.32364,
						1.32579,
						1.3332,
						1.33313,
						1.34116,
						1.33861,
						1.34862,
						1.36366};
	Double_t ex3[m] = {0.};
    Double_t ey3[m] = {0.};
    					
	const Int_t k = 25; // simulation
 
	Double_t x2[k] = {	1.300,
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
    Double_t y2[k] = {0.};
    Double_t ex2[k] = {0.};
    Double_t ey2[k] = {0.};

    //****************** Create Histo ************************************//
    auto c1 = new TCanvas("c1","c1", 1000, 600);
	c1->SetTitle("Mean cluster width for 2S mini-module");
	gStyle->SetOptStat(0);
	gPad->SetGridx(1);
    gPad->SetGridy(1);

    gPad->SetTitle("Mean cluster width for 2S mini-module");

    // loop on all files
    for (int j = 0; j < k; ++j)
    {
        //CopyFile(j);

        string file_path = "/media/matthieu/ssd1/Geant4/Full-size/Data_figure17-19B/data_";
        file_path += std::to_string(j);
        file_path += ".root";
        char const *pchar = file_path.c_str();
        TFile f(pchar, "read");
        //cout << "File " << j << " opened:" << endl;

	    //************* variable ***************//
	    const int NBR_STRIP = 1016;
	    const int MAX_CLUSTER_WIDTH = 3;
	    const double CLUSTER_WINDOW = 3.5;
	    const double THRESHOLD = 5.1975; // (14 * 375 * 3.62) / 1000000
	    
	    Int_t nbrCAT, nbrCA, nbrCBT, nbrCB; // for A and B detectors
	    Double_t mCWAT, mCWBT, mCWA, mCWB;
	    Int_t CA1, CA2, CA3, CA4, CA5;
	    Int_t CB1, CB2, CB3, CB4, CB5;
	    Double_t momentum;

	    //************ Get the TTree ******************************************
	    auto data = f.Get<TTree>("data");

	    // Get the number of entries in TTree
	    const int ENTRIES = data->GetEntries();
	    cout << std::scientific << "Number of entries: " << ENTRIES << endl;

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
	    int index = 0;
	    bool stop = true;
	    double cluster_width = 0.;
	    double nbr_clusters = 0.;
	    std::vector<double> mClusWidth(100, 0);
	    for (int k = 0; k < ENTRIES; k++)
	    {
	        // fill variables with datas from entry i
	        data->GetEntry(k);

	        std::vector<double> res(2, 0);
	        std::vector<double> res_A(10, 0);
	        std::vector<double> res_B(10, 0);

	        bool stub = CBC3(strip_A, strip_B, res, res_A, res_B, MAX_CLUSTER_WIDTH, CLUSTER_WINDOW, THRESHOLD);

	        nbrCAT = (int)res_A.at(5);
	        mCWAT = (double)res_A.at(6);
	        nbrCA = (int)res_A.at(7);
	        mCWA = (double)res_A.at(8);

	        nbrCBT = (int)res_B.at(5);
	        mCWBT = (double)res_B.at(6);
	        nbrCB = (int)res_B.at(7);
	        mCWB = (double)res_B.at(8);

	        if(!isnan(mCWAT) && !isnan(nbrCAT))
	        {
	        	cluster_width += (double) mCWAT * nbrCAT;
	        	nbr_clusters += (double) nbrCAT;
	        }
	        /*if(!isnan(mCWBT) && !isnan(nbrCBT))
	        {
	        	cluster_width += (double) mCWBT * nbrCBT;
	        	nbr_clusters += (double) nbrCBT;
	        }*/

	        count_loop += 1;
	        if (count_loop == ENTRIES /100)
	        {
	            count_loop = 0;
	            if(!isnan(cluster_width))
	            {
	            	mClusWidth.at(index) = cluster_width / nbr_clusters;
	            }else
	            {
	            	cout << "error! at file " << j << " index " << index << endl;
	            }
	            
	            index += 1;
	            cluster_width = 0.;
	            nbr_clusters = 0.;
	        }
	    }
	    //********************* fig 17 computation and printing **********************************//
	    // for A
	    double variance, deviation, average;
	    average = std::accumulate(mClusWidth.begin(), mClusWidth.end(), 0.0) / 100;
	    for (int i = 0; i < 100; ++i){variance += pow((mClusWidth.at(i) - average), 2);}
	    variance /= 99;
	    deviation = sqrt(variance);
	    cout << std::scientific << "File " << j << " mean cluster width:\t" << average << "\t" << deviation << endl;

	    //*********************** 
	    y2[j] = average;
	    ey2[j] = deviation;

	    // Close file when finished
	    f.Close();
	}   
	// Fill graphs
	TGraphErrors *gr1 = new TGraphErrors(n,x1,y1,ex1,ey1); // non-irradiated
    gr1->SetName("gr1");
    gr1->SetMarkerColor(12);
    gr1->SetMarkerStyle(24);
    gr1->SetMarkerSize(1.);

    /*TGraphErrors *gr3 = new TGraphErrors(m,x3,y3,ex3,ey3); // non-irradiated
    gr3->SetName("gr3");
    gr3->SetMarkerColor(12);
    gr3->SetMarkerStyle(25);
    gr3->SetMarkerSize(1.);*/

    TGraphErrors *gr2 = new TGraphErrors(k,x2,y2,ex2,ey2);
    gr2->SetName("gr2");
    gr2->SetMarkerColor(kBlue+2);
    gr2->SetMarkerStyle(20);
    gr2->SetMarkerSize(1.);

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(gr1);
    //mg->Add(gr3);
    mg->Add(gr2);
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
    yaxis->SetTitle("Largeur moyenne des clusters [strip]");
    yaxis->SetTitleFont(22);
	yaxis->SetTitleSize(0.05);
	yaxis->SetTitleOffset(0.9);

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

    auto legend = new TLegend(0.15,0.6,0.45,0.85);
    legend->AddEntry("gr1","Mini module","ep");
    //legend->AddEntry("gr3","Adam et al. - irr.","p");
    legend->AddEntry("gr2","Grand module","ep");
    legend->Draw();

    gPad->Modified();

    c1->SaveAs("figure17-full.pdf");
    //c1->SaveAs("figure17.png");

    SaveData(k, x2, y2, ey2);
}