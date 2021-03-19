void HistoTest1(){
	TFile *f = TFile::Open("/media/matthieu/ssd1/Geant4/Data/data-new.root");

	Double_t err;
	Bool_t no_corr_stub;

	TH1I *hErr;

	TTree *data = (TTree *)f->Get("data");

	Int_t entries = data->GetEntries();

	data->SetBranchAddress("errPT", &err);
	data->SetBranchAddress("no_corr_stub", &no_corr_stub);

	hErr = new TH1I("hErr", "Error in pT", 1000, -2000, 8000);

	hErr->GetXaxis()->SetTitle("MeV");

	for (int i = 0; i < entries; ++i)
	{
		data->GetEntry(i);
		if (!no_corr_stub)
		{
			hErr->Fill(err);
		}
		//hErr->Fill(err_xA * 1000);
	}
	//hErr->FillRandom("gaus",10000);
	hErr->Draw();
}
