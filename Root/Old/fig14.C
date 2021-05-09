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
    res_A.at(5) = clus_pos_A.size();
    res_B.at(5) = clus_pos_B.size();
    res_A.at(6) = std::accumulate(clus_size_A.begin(), clus_size_A.end(), 0.0) / clus_size_A.size();
    res_B.at(6) = std::accumulate(clus_size_B.begin(), clus_size_B.end(), 0.0) / clus_size_B.size();

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

void fig14() {
    // open data file
    TFile f("/media/matthieu/ssd1/Geant4/Data/fig14.root", "read");

    //************* variable ***************//
    const int NBR_STRIP = 101;
    const int MAX_CLUSTER_WIDTH = 3;
    const int CLUSTER_WINDOW = 5;
    const double THRESHOLD = 0.0144; // MeV -> = 4 * (1000 * 3.6 keV)
    
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
    double total_clusters_A, total_clusters_B;
    double total_clusters_AB;
    std::vector<double> mNbrClusA;
    std::vector<double> mNbrClusB;
    std::vector<double> mNbrClusAB;
    for (int k = 0; k < ENTRIES; k++)
    {
        // fill variables with datas from entry i
        data->GetEntry(k);

        std::vector<double> res_A(9, 0);
        std::vector<double> res_B(9, 0);

        bool stub = CBC2(strip_A, strip_B, res_A, res_B, MAX_CLUSTER_WIDTH, CLUSTER_WINDOW, THRESHOLD);

        nbrCAT = (int)res_A.at(5);

        nbrCBT = (int)res_B.at(5);

        // fig 14 (mean number of cluster / event)
        total_clusters_A += nbrCAT;
        total_clusters_B += nbrCBT;
        count_loop += 1;
        if (count_loop == ENTRIES /10)
        {
            count_loop = 0;
            mNbrClusA.push_back((double) total_clusters_A / (ENTRIES / 10));
            mNbrClusB.push_back((double) total_clusters_B / (ENTRIES / 10));
            total_clusters_A = total_clusters_B = 0.;
        }
    }
    //********************* fig 14 computation and printing **********************************//
    // for A
    double varianceA, deviationA, averageA;
    averageA = std::accumulate(mNbrClusA.begin(), mNbrClusA.end(), 0.0) / 10;
    for (int i = 0; i < 10; ++i){varianceA += pow((mNbrClusA.at(i) - averageA), 2);}
    varianceA /= 9;
    deviationA = sqrt(varianceA);
    cout << std::scientific << "Mean number of cluster/event for sensor A:\t" << averageA << "\t" << deviationA << endl;
    // for B
    double varianceB, deviationB, averageB;
    averageB = std::accumulate(mNbrClusB.begin(), mNbrClusB.end(), 0.0) / 10;
    for (int i = 0; i < 10; ++i){varianceB += pow((mNbrClusB.at(i) - averageB), 2);}
    varianceB /= 9;
    deviationB = sqrt(varianceB);
    cout << std::scientific << "Mean number of cluster/event for sensor B:\t" << averageB << "\t" << deviationB << endl;

    //*********************** 

    // Close file when finished
    f.Close();
}