void Adam2020() {
    TFile f("/media/matthieu/ssd1/Geant4/Data/data-new.root", "update");
    // rootcp -r /media/matthieu/ssd1/Geant4/Data/data.root /media/matthieu/ssd1/Geant4/Data/data-new.root

    //************ Create variables ******************************************
    int nbr_strip = 21;
    // double B_field = 3.8; // telsa
    int max_cluster_width = 3;
    int cluster_window = 5;
    double threshold = 0.0144; // MeV -> = 4 * (1000 * 3.6 keV)


    // double c_v = 2.99792458e+8;
    double nbr_pass = 0.;

    // for A and B detectors
    Int_t nbrCA;
    Int_t nbrCB;
    Double_t mCWA, mCWB;
    // For A: cluster of 1, 2 or more than 2 strips
    Int_t CA1 = 0;
    Int_t CA2 = 0;
    Int_t CAP = 0;
    Int_t CB1 = 0;
    Int_t CB2 = 0;
    Int_t CBP = 0;

    // for general
    Bool_t pass, no_corr_stub;

    // for existing datas in TTree
    Double_t x0, y0, z0, theta_i;

    //************ Get the TTree ******************************************
    auto data = f.Get<TTree>("data");

    //************ Create new Branches ************************************
    auto newBranch_A1 = data->Branch("nbrCA", &nbrCA, "nbrCA/I");
    auto newBranch_A2 = data->Branch("CA1", &CA1, "CA1/I");
    auto newBranch_A3 = data->Branch("CA2", &CA2, "CA2/I");
    auto newBranch_A4 = data->Branch("CAP", &CAP, "CAP/I");
    auto newBranch_A5 = data->Branch("mCWA", &mCWA, "mCWA/D");

    auto newBranch_B1 = data->Branch("nbrCB", &nbrCB, "nbrCB/I");
    auto newBranch_B2 = data->Branch("CB1", &CB1, "CB1/I");
    auto newBranch_B3 = data->Branch("CB2", &CB2, "CB2/I");
    auto newBranch_B4 = data->Branch("CBP", &CBP, "CBP/I");
    auto newBranch_B5 = data->Branch("mCWB", &mCWB, "mCWB/D");

    auto newBranch_T8 = data->Branch("no_corr_stub", &no_corr_stub, "no_corr_stub/O");

    // Get the number of entries in TTree
    Int_t entries = data->GetEntries();

    //**************** Set BranchAddress for datas recovery ***************
    // for strips
    double strip_A [nbr_strip];
    double strip_B [nbr_strip];
    for (int i = 0; i < nbr_strip; ++i)
    {
        string strip_name = "s" + std::to_string(i);
        char const *pchar = strip_name.c_str();
        data->SetBranchAddress(pchar, &strip_B[i]);
    }
    for (int i = nbr_strip; i < (2*nbr_strip); ++i)
    {
        string strip_name = "s" + std::to_string(i);
        char const *pchar = strip_name.c_str();
        data->SetBranchAddress(pchar, &strip_A[i-nbr_strip]);
    }

    // for other datas
    data->SetBranchAddress("theta_i", &theta_i);
    data->SetBranchAddress("x0", &x0);
    data->SetBranchAddress("y0", &y0);
    data->SetBranchAddress("z0", &z0);

    //****************** Main loop over all entries **********************
    for (int k = 0; k < entries; k++)
    {
    	// fill variables with datas from entry i
        data->GetEntry(k);

        //*************** process here *****************
        // identify clusters
        std::vector<double> clus_pos_A;
        std::vector<double> clus_pos_B;
        std::vector<double> clus_size_A;
        std::vector<double> clus_size_B;
        bool inside = false;
        double size;
        // clusters for sensor A
        for (int i = 0; i < nbr_strip; ++i)
        {
            if (strip_A[i] < threshold && !inside)        
            {} else if (strip_A[i] < threshold && inside)
            {
                if (size <= max_cluster_width)
                {
                    clus_size_A.push_back(size);
                    if(size == 1){
                        CA1 += 1;
                    } else if(size == 2){
                        CA2 += 1;
                    } else{
                        CAP += 1;
                    }
                    double pos = floor((i - 1) - (size / 2) + 0.5);
                    clus_pos_A.push_back(pos);
                    size = 0;
                    inside = false;
                } else {
                    size = 0;
                    inside = false;
                }
            } else if (strip_A[i] >= threshold && !inside)
            {
                size += 1;
                inside = true;
                if (i == (nbr_strip - 1))
                {
                    clus_size_A.push_back(size);
                    if(size == 1){
                        CA1 += 1;
                    } else if(size == 2){
                        CA2 += 1;
                    } else{
                        CAP += 1;
                    }
                    clus_pos_A.push_back(i);
                    size = 0;
                    inside = false;
                }
            } else if (strip_A[i] >= threshold && inside)
            {
                size += 1;
                if (i == (nbr_strip - 1))
                {
                    if (size <= max_cluster_width)
                    {
                        clus_size_A.push_back(size);
                        if(size == 1){
                            CA1 += 1;
                        } else if(size == 2){
                            CA2 += 1;
                        } else{
                            CAP += 1;
                        }
                        clus_pos_A.push_back(floor(i - (size / 2) + 0.5));
                        size = 0;
                        inside = false;
                    } else {
                        size = 0;
                        inside = false;
                    } 
                }
            }
        }
        // clusters for sensor B
        for (int i = 0; i < nbr_strip; ++i)
        {
            if (strip_B[i] < threshold && !inside)        
            {} else if (strip_B[i] < threshold && inside)
            {
                if (size <= max_cluster_width)
                {
                    clus_size_B.push_back(size);
                    if(size == 1){
                        CB1 += 1;
                    } else if(size == 2){
                        CB2 += 1;
                    } else{
                        CBP += 1;
                    }
                    double pos = floor((i - 1) - (size / 2) + 0.5);
                    clus_pos_B.push_back(pos);
                    size = 0;
                    inside = false;
                } else {
                    size = 0;
                    inside = false;
                }         
            } else if (strip_B[i] >= threshold && !inside)
            {
                size += 1;
                inside = true;
                if (i == (nbr_strip - 1))
                {
                    clus_size_B.push_back(size);
                    if(size == 1){
                        CB1 += 1;
                    } else if(size == 2){
                        CB2 += 1;
                    } else{
                        CBP += 1;
                    }
                    clus_pos_B.push_back(i);
                    size = 0;
                    inside = false;
                }
            } else if (strip_B[i] >= threshold && inside)
            {
                size += 1;
                if (i == (nbr_strip - 1))
                {
                    if (size <= max_cluster_width)
                    {
                        clus_size_B.push_back(size);
                        if(size == 1){
                            CB1 += 1;
                        } else if(size == 2){
                            CB2 += 1;
                        } else{
                            CBP += 1;
                        }
                        clus_pos_B.push_back(floor(i - (size / 2) + 0.5));
                        size = 0;
                        inside = false;
                    } else {
                        size = 0;
                        inside = false;
                    }  
                }
            }
        }
        // find correlations between two layers
        
        std::vector<double> stub_pos_a;
        std::vector<double> stub_pos_b;
        std::vector<double> stub_angle;
        int index_of_stub = 0;
        no_corr_stub = true;
        if (clus_pos_A.size() != 0 && clus_pos_B.size() != 0)
        {
            // stub finding
            for (int i = 0; i < clus_pos_A.size(); ++i)
            {
                double seed = clus_pos_A.at(i);
                for (int j = 0; j < clus_pos_B.size(); ++j)
                {
                    double candidat = clus_pos_B.at(j);
                    if (abs(candidat - seed) <= cluster_window)
                    {
                        stub_pos_a.push_back(seed);
                        stub_pos_b.push_back(candidat);
                        stub_angle.push_back(abs(candidat - seed));
                    } else {}
                }
            }
            if (stub_pos_a.size() != 0)
            {
                no_corr_stub = false;
            }
        }
        if (!no_corr_stub) { // if at least one stub is found
            index_of_stub = std::min_element(stub_angle.begin(), stub_angle.end()) - stub_angle.begin(); // index of lowest angle stub
            
            nbrCA = clus_pos_A.size();
            nbrCB = clus_pos_B.size();
            mCWA = std::accumulate(clus_size_A.begin(), clus_size_A.end(), 0.0) / clus_size_A.size();
            mCWB = std::accumulate(clus_size_B.begin(), clus_size_B.end(), 0.0) / clus_size_B.size();
    
        } else {
            nbrCA = clus_pos_A.size();
            nbrCB = clus_pos_B.size();
            if (clus_size_A.size() != 0)
            {
                mCWA = std::accumulate(clus_size_A.begin(), clus_size_A.end(), 0.0) / clus_size_A.size();
            } else {
                mCWA = 0.;
            }
            if (clus_size_B.size() != 0)
            {
                mCWB = std::accumulate(clus_size_B.begin(), clus_size_B.end(), 0.0) / clus_size_B.size();
            } else {
                mCWB = 0.;
            }
        }

    	// Fill Branches
        newBranch_A1->Fill();
        newBranch_A2->Fill();
        newBranch_A3->Fill();
        newBranch_A4->Fill();
        newBranch_A5->Fill();

        newBranch_B1->Fill();
        newBranch_B2->Fill();
        newBranch_B3->Fill();
        newBranch_B4->Fill();
        newBranch_B5->Fill();
        
        newBranch_T8->Fill();

        if(!no_corr_stub) nbr_pass += 1;
        
        // reset variables (if needed only)
        CA1 = CA2 = CAP = CB1 = CB2 = CBP = 0;
    }
    cout << "Stub efficiency: " << nbr_pass / entries << endl;

    // Overwrite file with new data 
    data->Write("", TObject::kOverwrite);
}