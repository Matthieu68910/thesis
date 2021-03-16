/// \file eventgenerator/exgps/src/HistoManager.cc
/// \brief Implementation of the HistoManager class

#include "HistoManager.hh"
#include "G4UnitsTable.hh"
#include "B1DetectorConstruction.hh"
#include "G4RunManager.hh"
#include <string>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
  : fFileName("/media/matthieu/ssd1/Geant4/Data/data"),
    strip_nbr(0)
{
  Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  //
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  analysis->SetVerboseLevel(0);
  analysis->SetNtupleMerging(true);
  //analysis->SetNtupleDirectoryName("/media/matthieu/ssd1/Geant4/Data/data");
  
  analysis->SetFileName(fFileName);
  analysis->SetActivation(true);     //enable inactivation of histos, nTuples

  // nTuples
  //
  //analysis->SetNtupleDirectoryName("ntuple");
  analysis->CreateNtuple("data", "Primary Particle Tuple");
  analysis->CreateNtupleDColumn("x_i");         //column 0
  analysis->CreateNtupleDColumn("y_i");         //column 1
  analysis->CreateNtupleDColumn("z_i");         //column 2
  analysis->CreateNtupleDColumn("theta_i");     //column 3
  analysis->CreateNtupleDColumn("phi_i");       //column 4
  analysis->CreateNtupleDColumn("E_i");         //column 5
  analysis->CreateNtupleDColumn("x_0");         //column 6
  analysis->CreateNtupleDColumn("y_0");         //column 7
  analysis->CreateNtupleDColumn("x_A");         //column 8
  analysis->CreateNtupleDColumn("y_A");         //column 9
  analysis->CreateNtupleDColumn("z_A");         //column 10
  analysis->CreateNtupleDColumn("x_B");         //column 11
  analysis->CreateNtupleDColumn("y_B");         //column 12
  analysis->CreateNtupleDColumn("z_B");         //column 13

  if (strip_nbr == 0) {
    const B1DetectorConstruction* detectorConstruction
          = static_cast<const B1DetectorConstruction*>
            (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    strip_nbr = detectorConstruction->GetStripNumber();
  }

  for (int x = 0; x < (2*strip_nbr); x++){
    G4String name = "s" + std::to_string(x);
    analysis->CreateNtupleDColumn(name);        // strip x
  }
  analysis->FinishNtuple();
  
  analysis->SetNtupleActivation(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
