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
  analysis->CreateNtupleDColumn("x0");         //column 0
  analysis->CreateNtupleDColumn("y0");         //column 1
  analysis->CreateNtupleDColumn("z0");         //column 2
  analysis->CreateNtupleDColumn("theta_i");     //column 3
  analysis->CreateNtupleDColumn("momentum");     //column 4

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
