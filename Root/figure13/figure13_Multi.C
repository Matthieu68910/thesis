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
	vector<double> &res_A, // [0:4]-> 1-5 strip wide clusters,[5]-> number of clusters (tot), [6]-> mean cluster width (tot), [7:8]-> accepted, [9] -> nbr hits
	vector<double> &res_B,
	const int MAX_CLUSTER_WIDTH = 3,
	const int CLUSTER_WINDOW = 5,
	double THRESHOLD = 5.1975
    //double kill_value = 0.04331
	){

	const int NBR_STRIP = strip_A.size();

    double noise = 0.;
	// Clusters in sensor A
	std::vector<double> clus_pos_A;
	std::vector<double> clus_size_A;
	bool inside = false;
	int size = 0;
    int nbr_hits_A = 0;
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

void figure13_Multi() {
	// data for Adam2020
    const int NBR_BINS = 14;
    const int NBR_BINS_IRR = 16;

    // array creation
    Double_t x1[NBR_BINS] = {	0.69750,
                                1.46438,
                                2.23125,
                                2.95800,
                                5.21850,
                                5.98538,
                                7.47874,
                                9.69881,
                                11.23268,
                                16.48005,
                                19.46704,
                                20.97015,
                                22.45403,
                                23.94750}; // in keV !!!
    Double_t y1[NBR_BINS] = {   4.9381E+02,
                                1.4626E+02,
                                2.7247E+01,
                                4.2766E+00,
                                1.1251E+00,
                                1.1216E+00,
                                1.1135E+00,
                                1.0999E+00,
                                1.0909E+00,
                                1.0647E+00,
                                1.0335E+00,
                                9.7480E-01,
                                8.6380E-01,
                                7.1770E-01,
                                };
    Double_t ex1[NBR_BINS] = {  1.4583E+00,
                                1.1206E+00,
                                4.7950E-01,
                                1.8336E-01,
                                8.0309E-02,
                                7.9082E-02,
                                7.2938E-02,
                                6.0495E-02,
                                5.3221E-02,
                                4.0388E-02,
                                4.0060E-02,
                                4.5249E-02,
                                5.1223E-02,
                                6.2294E-02
                                };
    Double_t ey1[NBR_BINS] = {  1.4583E+00,
                                1.1206E+00,
                                4.7950E-01,
                                1.8336E-01,
                                8.0309E-02,
                                7.9082E-02,
                                7.2938E-02,
                                6.0495E-02,
                                5.3221E-02,
                                4.0388E-02,
                                4.0060E-02,
                                4.5249E-02,
                                5.1223E-02,
                                6.2294E-02
                                };

    Double_t x2[NBR_BINS] = {   0.73800,
                                1.50488,
                                2.23125,
                                2.99813,
                                5.25863,
                                5.98538,
                                7.47874,
                                9.73916,
                                11.23268,
                                16.49768,
                                19.50739,
                                21.00090,
                                22.49438,
                                23.98785}; // in keV !!!
    Double_t y2[NBR_BINS] = {   4.6856E+02,
                                1.3529E+02,
                                2.7232E+01,
                                3.9563E+00,
                                1.1620E+00,
                                1.1577E+00,
                                1.1464E+00,
                                1.1258E+00,
                                1.1121E+00,
                                1.0772E+00,
                                1.0446E+00,
                                9.8210E-01,
                                8.7070E-01,
                                7.1990E-01
                                };
    Double_t ex2[NBR_BINS] = {  1.3936E+00,
                                1.1558E+00,
                                4.9097E-01,
                                2.1856E-01,
                                1.1649E-01,
                                1.1456E-01,
                                1.0593E-01,
                                9.2629E-02,
                                8.1777E-02,
                                5.6373E-02,
                                5.1843E-02,
                                5.3904E-02,
                                6.2897E-02,
                                6.5912E-02
                                };
    Double_t ey2[NBR_BINS] = {0.};

    // adam D0
    Double_t x3[NBR_BINS] = {   0.69750,
                                1.46438,
                                2.23125,
                                2.95800,
                                5.21850,
                                5.98538,
                                7.47874,
                                9.69881,
                                11.23268,
                                16.48005,
                                19.46704,
                                20.97015,
                                22.45403,
                                23.94750}; // in keV !!!
    Double_t y3[NBR_BINS] = {   4.9378E+02,
                                1.4621E+02,
                                2.7183E+01,
                                4.2168E+00,
                                1.0647E+00,
                                1.0627E+00,
                                1.0594E+00,
                                1.0517E+00,
                                1.0460E+00,
                                1.0329E+00,
                                1.0051E+00,
                                9.4920E-01,
                                8.3810E-01,
                                6.9090E-01
                                };
    Double_t ex3[NBR_BINS] = {  1.4664E+00,
                                1.1246E+00,
                                4.7346E-01,
                                1.8076E-01,
                                4.0363E-02,
                                3.9128E-02,
                                3.7600E-02,
                                3.2037E-02,
                                2.8955E-02,
                                2.2487E-02,
                                2.3932E-02,
                                3.3505E-02,
                                4.1795E-02,
                                5.2399E-02
                                };
    Double_t ey3[NBR_BINS] = {  1.4664E+00,
                                1.1246E+00,
                                4.7346E-01,
                                1.8076E-01,
                                4.0363E-02,
                                3.9128E-02,
                                3.7600E-02,
                                3.2037E-02,
                                2.8955E-02,
                                2.2487E-02,
                                2.3932E-02,
                                3.3505E-02,
                                4.1795E-02,
                                5.2399E-02};

    // adam D1
    Double_t x4[NBR_BINS] = {   0.73800,
                                1.50488,
                                2.23125,
                                2.99813,
                                5.25863,
                                5.98538,
                                7.47874,
                                9.73916,
                                11.23268,
                                16.49768,
                                19.50739,
                                21.00090,
                                22.49438,
                                23.98785}; // in keV !!!
    Double_t y4[NBR_BINS] = {   4.6854E+02,
                                1.3523E+02,
                                2.7169E+01,
                                3.8931E+00,
                                1.1025E+00,
                                1.0984E+00,
                                1.0921E+00,
                                1.0803E+00,
                                1.0722E+00,
                                1.0492E+00,
                                1.0166E+00,
                                9.5410E-01,
                                8.4560E-01,
                                6.9930E-01
                                };
    Double_t ex4[NBR_BINS] = {  1.3972E+00,
                                1.1291E+00,
                                4.9971E-01,
                                1.8169E-01,
                                6.2448E-02,
                                5.9742E-02,
                                5.5438E-02,
                                4.7746E-02,
                                4.2866E-02,
                                3.2214E-02,
                                2.8292E-02,
                                3.7473E-02,
                                4.4593E-02,
                                5.5619E-02
                                };
    Double_t ey4[NBR_BINS] = {  1.3972E+00,
                                1.1291E+00,
                                4.9971E-01,
                                1.8169E-01,
                                6.2448E-02,
                                5.9742E-02,
                                5.5438E-02,
                                4.7746E-02,
                                4.2866E-02,
                                3.2214E-02,
                                2.8292E-02,
                                3.7473E-02,
                                4.4593E-02,
                                5.5619E-02};

    // adam D0_IRR
    Double_t x5[NBR_BINS] = {   0.69750,
                                1.46438,
                                2.23125,
                                2.95800,
                                5.21850,
                                5.98538,
                                7.47874,
                                9.69881,
                                11.23268,
                                16.48005,
                                19.46704,
                                20.97015,
                                22.45403,
                                23.94750}; // in keV !!!
    Double_t y5[NBR_BINS] = {   4.9377E+02,
                                1.4620E+02,
                                2.7187E+01,
                                4.2165E+00,
                                1.0665E+00,
                                1.0644E+00,
                                1.0604E+00,
                                1.0543E+00,
                                1.0496E+00,
                                1.0355E+00,
                                1.0076E+00,
                                9.4680E-01,
                                8.3690E-01,
                                6.7960E-01
                                };
    Double_t ex5[NBR_BINS] = {  1.4648E+00,
                                1.1293E+00,
                                4.7364E-01,
                                1.7664E-01,
                                4.6414E-02,
                                4.5711E-02,
                                4.3552E-02,
                                4.1615E-02,
                                3.8321E-02,
                                2.7464E-02,
                                2.9987E-02,
                                3.9283E-02,
                                4.6725E-02,
                                5.7541E-02
                                };
    Double_t ey5[NBR_BINS] = {  1.4648E+00,
                                1.1293E+00,
                                4.7364E-01,
                                1.7664E-01,
                                4.6414E-02,
                                4.5711E-02,
                                4.3552E-02,
                                4.1615E-02,
                                3.8321E-02,
                                2.7464E-02,
                                2.9987E-02,
                                3.9283E-02,
                                4.6725E-02,
                                5.7541E-02};

    // adam D1_IRR
    Double_t x6[NBR_BINS] = {   0.73800,
                                1.50488,
                                2.23125,
                                2.99813,
                                5.25863,
                                5.98538,
                                7.47874,
                                9.73916,
                                11.23268,
                                16.49768,
                                19.50739,
                                21.00090,
                                22.49438,
                                23.98785}; // in keV !!!
    Double_t y6[NBR_BINS] = {   4.6852E+02,
                                1.3522E+02,
                                2.7157E+01,
                                3.8824E+00,
                                1.0921E+00,
                                1.0894E+00,
                                1.0840E+00,
                                1.0727E+00,
                                1.0656E+00,
                                1.0464E+00,
                                1.0126E+00,
                                9.5120E-01,
                                8.4390E-01,
                                6.9700E-01
                                };
    Double_t ex6[NBR_BINS] = {  1.3963E+00,
                                1.1539E+00,
                                4.9567E-01,
                                1.8543E-01,
                                5.1704E-02,
                                5.1909E-02,
                                4.9093E-02,
                                4.2660E-02,
                                3.9653E-02,
                                2.8587E-02,
                                2.6232E-02,
                                3.3972E-02,
                                4.8386E-02,
                                5.8990E-02
                                };
    Double_t ey6[NBR_BINS] = {  1.3963E+00,
                                1.1539E+00,
                                4.9567E-01,
                                1.8543E-01,
                                5.1704E-02,
                                5.1909E-02,
                                4.9093E-02,
                                4.2660E-02,
                                3.9653E-02,
                                2.8587E-02,
                                2.6232E-02,
                                3.3972E-02,
                                4.8386E-02,
                                5.8990E-02};


    //****************** Create Histo ************************************//
    auto c1 = new TCanvas("c1","c1",1000,600);
	c1->SetTitle("Figure 1: threshold variation");
	gStyle->SetOptStat(0);
	gPad->SetGridx(1);
    gPad->SetGridy(1);
    gPad->SetLogy();

    gPad->SetTitle("Figure 1: threshold variation");



    for (int i = 0; i < NBR_BINS; ++i)
    {
        x1[i] *= 1000.;
        x2[i] *= 1000.;
        x3[i] *= 1000.;
        x4[i] *= 1000.;
        x5[i] *= 1000.;
        x6[i] *= 1000.;
    }

    TGraphErrors *gr1 = new TGraphErrors(NBR_BINS,x1,y1,ex1,ey1);
    gr1->SetName("gr1");
    gr1->SetMarkerColor(kBlue+2);
    gr1->SetMarkerStyle(20);
    gr1->SetMarkerSize(1.0);

    TGraphErrors *gr2 = new TGraphErrors(NBR_BINS,x2,y2,ex2,ey2);
    gr2->SetName("gr2");
    gr2->SetMarkerColor(kBlue+2);
    gr2->SetMarkerStyle(25);
    gr2->SetMarkerSize(1.0);

    TGraphErrors *gr3 = new TGraphErrors(NBR_BINS,x3,y3,ex3,ey3);
    gr3->SetName("gr3");
    gr3->SetMarkerColor(kRed+2);
    gr3->SetMarkerStyle(20);
    gr3->SetMarkerSize(1.0);

    TGraphErrors *gr4 = new TGraphErrors(NBR_BINS,x4,y4,ex4,ey4);
    gr4->SetName("gr4");
    gr4->SetMarkerColor(kRed+2);
    gr4->SetMarkerStyle(25);
    gr4->SetMarkerSize(1.0);

    TGraphErrors *gr5 = new TGraphErrors(NBR_BINS,x5,y5,ex5,ey5);
    gr5->SetName("gr5");
    gr5->SetMarkerColor(kGreen+2);
    gr5->SetMarkerStyle(20);
    gr5->SetMarkerSize(1.0);

    TGraphErrors *gr6 = new TGraphErrors(NBR_BINS,x6,y6,ex6,ey6);
    gr6->SetName("gr6");
    gr6->SetMarkerColor(kGreen+2);
    gr6->SetMarkerStyle(25);
    gr6->SetMarkerSize(1.0);

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(gr3);
    mg->Add(gr4);
    mg->Add(gr5);
    mg->Add(gr6);
    mg->Add(gr1);
    mg->Add(gr2);
    mg->SetTitle("");
    mg->Draw("AP");

    TAxis *xaxis = mg->GetXaxis();
    TAxis *yaxis = mg->GetYaxis();
    xaxis->SetMaxDigits(3);
    xaxis->SetLabelFont(42);
    xaxis->SetLabelSize(0.04);
    xaxis->SetTitle("Seuil [e]");
    xaxis->SetTitleFont(22);
    xaxis->SetTitleSize(0.05);
    xaxis->SetTitleOffset(0.95);
    xaxis->SetRangeUser(0., 25000);
    yaxis->SetLabelFont(42);
    yaxis->SetLabelSize(0.04);
    yaxis->SetTitle("Nombre moyen de hit/evenement");
    yaxis->SetTitleFont(22);
    yaxis->SetTitleSize(0.05);
    yaxis->SetTitleOffset(0.9);

    c1->RedrawAxis();

    TF1* f1 = new TF1("f1", "x", 0, 25000);
    TGaxis* A1 = new TGaxis(0., yaxis->GetXmax(), 25000, yaxis->GetXmax(), "f1", 510, "-");
    A1->SetLabelFont(42);
    A1->SetLabelSize(0.04);
    A1->SetTitle("Seuil [keV]");
    A1->SetTitleFont(22);
    A1->SetTitleSize(0.05);
    A1->SetTitleOffset(0.95);
    A1->ChangeLabel(1, -1, -1, -1, -1, -1, "0.0");
    A1->ChangeLabel(2, -1, -1, -1, -1, -1, "18.1");
    A1->ChangeLabel(3, -1, -1, -1, -1, -1, "36.2");
    A1->ChangeLabel(4, -1, -1, -1, -1, -1, "54.3");
    A1->ChangeLabel(5, -1, -1, -1, -1, -1, "72.4");
    A1->ChangeLabel(6, -1, -1, -1, -1, -1, "90.5");
    //A1->ChangeLabel(7, -1, -1, -1, -1, -1, "1.4");
    A1->Draw("SAME");

    auto legend = new TLegend(0.55,0.9,0.9,0.6);
    legend->AddEntry("gr1","B = 0 - capteur 1","ep");
    legend->AddEntry("gr2","B = 0 - capteur 2","ep");
    legend->AddEntry("gr3","By = 3T - capteur 1","ep");
    legend->AddEntry("gr4","By = 3T - capteur 2","ep");
    legend->AddEntry("gr5","Bz = 3T - capteur 1","ep");
    legend->AddEntry("gr6","Bz = 3T - capteur 2","ep");
    legend->Draw();

    gPad->Modified();
    //*********************** 

    c1->SaveAs("figure13.pdf");
    // Close file when finished
    //f.Close();
}