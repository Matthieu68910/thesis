void SCurve(){
	TFile *f = TFile::Open("/media/matthieu/ssd1/Geant4/Data/S-curve/2.8/data-new.root");

	Bool_t pass, no_corr_stub;

	TTree *data = (TTree *)f->Get("data");

	Int_t entries = data->GetEntries();
	cout << "Nbr entries: " << entries << endl;

	data->SetBranchAddress("pass", &pass);
	data->SetBranchAddress("no_corr_stub", &no_corr_stub);

	std::vector<double> frac_pass;
	double total_average;
	double nbr_no_corr_stub = 0;

	for (int i = 0; i < 10; ++i)
	{
		double nbr_pass = 0;
		
		for (int j = 0; j < entries / 10; ++j)
		{
			data->GetEntry((i * 100000) + j);
			if (pass){nbr_pass += 1;}
			if (no_corr_stub){nbr_no_corr_stub += 1;}
		}
		frac_pass.push_back(nbr_pass / 100000);
		total_average = total_average + (nbr_pass / 100000);
		cout << "Loop nbr " << i << ", nbr_pass = " << nbr_pass << endl;
	}
	total_average = total_average / 10;

	cout << "total_average = " << total_average << endl;

	double variance, deviation;
	for (int i = 0; i < 10; ++i)
	{
		variance += pow((frac_pass.at(i) - total_average), 2);
	}
	variance /= 9;
	deviation = sqrt(variance);

	cout << "Average pass frac.: " << total_average << endl;
	cout << "Standard deviation: " << deviation << endl;
	cout << "Total incorect stub: " << (nbr_no_corr_stub / 1000000) << endl;
	cout << total_average << "\t" << deviation << "\t" << (nbr_no_corr_stub / 1000000) << endl;
}