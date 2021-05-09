void TestFunc2() {
    TFile f("/media/matthieu/ssd1/Geant4/Data/data-new.root", "update");
    // rootcp -r /media/matthieu/ssd1/Geant4/Data/data.root /media/matthieu/ssd1/Geant4/Data/data-new.root

    //************ Create variables ******************************************
    int nbr_strip = 21;
    double sensor_sep = 1.6; // mm
    double sensor_dist = 600.; // mm
    double B_field = 3.8; // telsa
    double PT_real = 2000.;
    double PT_limmit = 2000.;

    double c_v = 2.99792458e+8;

    // for A and B detectors
    Int_t nbrCA;
    Int_t nbrCB;
    Double_t mCWA, mCWB;
    Double_t EtotA, EtotB;
    Double_t xAmes, xBmes;
    Double_t err_xA, err_xB;

    // for center
    Double_t x0mes, err_x0, theta_mes, err_theta, PTmes, errPT;
    Bool_t pass, no_corr_stub;

    // for existing datas in TTree
    Double_t theta_i, phi_i, E_i, x_0, y_0, x_A, y_A, z_A, x_B, y_B, z_B;

    //************ Get the TTree ******************************************
    auto data = f.Get<TTree>("data");

    //************ Create new Branches ************************************
    auto newBranch_A1 = data->Branch("nbrCA", &nbrCA, "nbrCA/I");
    auto newBranch_A2 = data->Branch("mCWA", &mCWA, "mCWA/D");
    auto newBranch_A3 = data->Branch("EtotA", &EtotA, "EtotA/D");
    auto newBranch_A4 = data->Branch("xAmes", &xAmes, "xAmes/D");
    auto newBranch_A5 = data->Branch("err_xA", &err_xA, "err_xA/D");

    auto newBranch_B1 = data->Branch("nbrCB", &nbrCB, "nbrCB/I");
    auto newBranch_B2 = data->Branch("mCWB", &mCWB, "mCWB/D");
    auto newBranch_B3 = data->Branch("EtotB", &EtotB, "EtotB/D");
    auto newBranch_B4 = data->Branch("xBmes", &xBmes, "xBmes/D");
    auto newBranch_B5 = data->Branch("err_xB", &err_xB, "err_xB/D");

    auto newBranch_T1 = data->Branch("x0mes", &x0mes, "x0mes/D");
    auto newBranch_T2 = data->Branch("err_x0", &err_x0, "err_x0/D");
    auto newBranch_T3 = data->Branch("theta_mes", &theta_mes, "theta_mes/D");
    auto newBranch_T4 = data->Branch("err_theta", &err_theta, "err_theta/D");
    auto newBranch_T5 = data->Branch("PTmes", &PTmes, "PTmes/D");
    auto newBranch_T6 = data->Branch("errPT", &errPT, "errPT/D");
    auto newBranch_T7 = data->Branch("pass", &pass, "pass/O");
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
    data->SetBranchAddress("phi_i", &phi_i);
    data->SetBranchAddress("E_i", &E_i);
    data->SetBranchAddress("x_0", &x_0);
    data->SetBranchAddress("y_0", &y_0);
    data->SetBranchAddress("x_A", &x_A);
    data->SetBranchAddress("y_A", &y_A);
    data->SetBranchAddress("z_A", &z_A);
    data->SetBranchAddress("x_B", &x_B);
    data->SetBranchAddress("y_B", &y_B);
    data->SetBranchAddress("z_B", &z_B);

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
            if (strip_A[i] == 0 && !inside)        
            {} else if (strip_A[i] == 0 && inside)
            {
                if (size <= 4)
                {
                    clus_size_A.push_back(size);
                    double pos = (i - 1) - (size / 2) + 0.5;
                    clus_pos_A.push_back(pos);
                    size = 0;
                    inside = false;
                } else {
                    size = 0;
                    inside = false;
                }
            } else if (strip_A[i] != 0 && !inside)
            {
                size += 1;
                inside = true;
                if (i == (nbr_strip - 1))
                {
                    clus_size_A.push_back(size);
                    clus_pos_A.push_back(i);
                    size = 0;
                    inside = false;
                }
            } else if (strip_A[i] != 0 && inside)
            {
                size += 1;
                if (i == (nbr_strip - 1))
                {
                    if (size <= 4)
                    {
                        clus_size_A.push_back(size);
                        clus_pos_A.push_back(i - (size / 2) + 0.5);
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
            if (strip_B[i] == 0 && !inside)        
            {} else if (strip_B[i] == 0 && inside)
            {
                if (size <= 4)
                {
                    clus_size_B.push_back(size);
                    double pos = (i - 1) - (size / 2) + 0.5;
                    clus_pos_B.push_back(pos);
                    size = 0;
                    inside = false;
                } else {
                    size = 0;
                    inside = false;
                }         
            } else if (strip_B[i] != 0 && !inside)
            {
                size += 1;
                inside = true;
                if (i == (nbr_strip - 1))
                {
                    clus_size_B.push_back(size);
                    clus_pos_B.push_back(i);
                    size = 0;
                    inside = false;
                }
            } else if (strip_B[i] != 0 && inside)
            {
                size += 1;
                if (i == (nbr_strip - 1))
                {
                    if (size <= 4)
                    {
                        clus_size_B.push_back(size);
                        clus_pos_B.push_back(i - (size / 2) + 0.5);
                        size = 0;
                        inside = false;
                    } else {
                        size = 0;
                        inside = false;
                    }  
                }
            }
        }
        
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
                    if (abs(candidat - seed) <= 7)
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
        if (!no_corr_stub) {
            index_of_stub = std::min_element(stub_angle.begin(), stub_angle.end()) - stub_angle.begin();
            /*while (stub_angle.at(index_of_stub) == 0 && stub_angle.size() > 1)
            {
                stub_pos_a.erase(stub_pos_a.begin() + index_of_stub);
                stub_pos_b.erase(stub_pos_b.begin() + index_of_stub);
                stub_angle.erase(stub_angle.begin() + index_of_stub);
                index_of_stub = std::min_element(stub_angle.begin(), stub_angle.end()) - stub_angle.begin();
            }*/
            nbrCA = clus_pos_A.size();
            nbrCB = clus_pos_B.size();
            mCWA = std::accumulate(clus_size_A.begin(), clus_size_A.end(), 0.0) / clus_size_A.size();
            mCWB = std::accumulate(clus_size_B.begin(), clus_size_B.end(), 0.0) / clus_size_B.size();
            EtotA = accumulate(strip_A, strip_A+21, 0.0);
            EtotB = accumulate(strip_B, strip_B+21, 0.0);

            try {
                xAmes = (stub_pos_a.at(index_of_stub) * 0.09) - 0.9;
            } catch (const std::exception& e) {
                cout << k << "a" << " - Message : " << e.what() << endl;
                cout << "Index of stub = " << index_of_stub << endl;
                cout << "Vector stub_pos_a size = " << stub_pos_a.size() << endl;
            }
            try {
                xBmes = (stub_pos_b.at(index_of_stub) * 0.09) - 0.9;
            } catch (const std::exception& e) {
                cout << k << "b" << " - Message : " << e.what() << endl;
                cout << "Vector stub_pos_b size = " << stub_pos_b.size() << endl << endl;
            }
            err_xA = xAmes - x_A;
            err_xB = xBmes - x_B;

            // theta et err_theta en fonction de theta_0 !! (pas theta_i)
            theta_mes = atan((xBmes - xAmes) / sensor_sep);
            double theta_0 = atan((x_B - x_A) / sensor_sep);
            err_theta = theta_mes - theta_0;

            //************** x0mes computation ****************
            // triangle half-perimeter computation p = (a+b+c)/2
            double a = sqrt(pow((xBmes - xAmes), 2) + pow((z_B - z_A), 2));
            double b = sqrt(pow(xAmes, 2) + pow((z_A - sensor_dist), 2));
            double c = sqrt(pow(xBmes, 2) + pow((z_B - sensor_dist), 2));
            double p = (a + b + c) / 2;
            // radius computation with R = (abc)/(4*sqrt(p*(p-a)*(p-b)*(p-c)))
            double radius = (a * b * c) / (4 * sqrt(p * (p - a) * (p - b) * (p - c)));
            // deviation x 
            double l = sqrt(pow((xBmes - xAmes), 2) + pow(sensor_sep, 2)) / 2;
            double x = radius - radius * cos(0.5 * acos(1 - pow(l, 2) / (2 * pow(radius, 2))));
            // final displacement in x
            double dev_x = x / cos(theta_mes);
            // x0mes before correction
            double x0mes_before = (xAmes + xBmes) / 2;
            // direction of correction, and x0mes final determination
            if (xBmes >= xAmes)
            {
                x0mes = x0mes_before - dev_x;
            } else {
                x0mes = x0mes_before + dev_x;
            }
            // error computation
            err_x0 = x0mes - x_0;

            //****************** Momentum computation *********************
            PTmes = ((radius / 1000) * c_v * B_field) * 1e-6;
            /*if (PTmes > 5000.)
            {
                cout << "PTmes > 3000. at index " << k << endl;
            }*/
            errPT = PTmes - PT_real;
            if (PTmes >= PT_limmit)
            {
                pass = true;
            } else {
                pass = false;
            }
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
            EtotA = accumulate(strip_A, strip_A+21, 0.0);
            EtotB = accumulate(strip_B, strip_B+21, 0.0);
            xAmes = 0.;
            err_xA = 0.;
            xBmes = 0.;
            err_xB = 0.;
            x0mes = 0.;
            err_x0 = 0.;
            theta_mes = 0.;
            err_theta = 0.;
            PTmes = 0.;
            errPT = 0.;
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
        newBranch_T1->Fill();
        newBranch_T2->Fill();
        newBranch_T3->Fill();
        newBranch_T4->Fill();
        newBranch_T5->Fill();
        newBranch_T6->Fill();
        newBranch_T7->Fill();
        newBranch_T8->Fill();
        
        // reset variables (if needed only)
    }

    // Overwrite file with new data 
    data->Write("", TObject::kOverwrite);
}