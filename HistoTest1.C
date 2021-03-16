void HistoTest1(){
	TFile *f = TFile::Open("/media/matthieu/ssd1/Geant4/Data/data-new.root");

	Double_t err_xA;
	Bool_t no_corr_stub;

	TH1D *hErrXA;

	TTree *data = (TTree *)f->Get("data");

	Int_t entries = data->GetEntries();

	data->SetBranchAddress("err_xA", &err_xA);
	data->SetBranchAddress("no_corr_stub", &no_corr_stub);

	hErrXA = new TH1D("hErrXA", "Error in x_A", 1000, -50, 50);

	hErrXA->GetXaxis()->SetTitle("µm");

	for (int i = 0; i < entries; ++i)
	{
		data->GetEntry(i);
		if (!no_corr_stub)
		{
			hErrXA->Fill(err_xA * 1000);
		}
		//hErrXA->Fill(err_xA * 1000);
	}
	//hErrXA->FillRandom("gaus",10000);
	hErrXA->Draw();
}