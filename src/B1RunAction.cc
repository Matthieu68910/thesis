/// \file B1RunAction.cc
/// \brief Implementation of the B1RunAction class

#include "B1RunAction.hh"
#include "B1PrimaryGeneratorAction.hh"
#include "B1DetectorConstruction.hh"
#include "B1EventAction.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <ctime>
#include <chrono>


B1RunAction::B1RunAction()
: G4UserRunAction(),
  fHistoManager(0)
{ 
  // histo manager
  //
  fHistoManager = new HistoManager(); 
  // add new units for dose
  // 
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;  
  const G4double picogray  = 1.e-12*gray;
   
  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray); 
}

B1RunAction::~B1RunAction()
{
  delete fHistoManager;
}

void B1RunAction::BeginOfRunAction(const G4Run*)
{ 
  //histograms
  //
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();

  analysis->OpenFile();

  if (IsMaster()) {
   start = std::chrono::system_clock::now();
  }
  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  // reset accumulables to their initial values
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();

}

void B1RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  // Merge accumulables 
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

  const B1PrimaryGeneratorAction* generatorAction
   = static_cast<const B1PrimaryGeneratorAction*>
     (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  if (generatorAction)
  {
    const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += " of ";
    G4double particleEnergy = particleGun->GetParticleEnergy();
    runCondition += G4BestUnit(particleEnergy,"Energy");
  }


  // Print
  //  
  if (IsMaster()) {
    G4cout
     << G4endl
     << "--------------------End of Global Run-----------------------";

    auto end = std::chrono::system_clock::now();
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    std::chrono::duration<double> elapsed_seconds = end-start;

    // Save elapsed Time

    std::ofstream MyFile("timer.txt", std::ios_base::app);
    MyFile << "Simulation started on " << std::ctime(&end_time);
    MyFile << "The run consists of " << nofEvents << " "<< runCondition << "\n";
    MyFile << "\t-----> elapsed time: " << elapsed_seconds.count() << " secondes\n\n";
    MyFile.close();
  }
  else {
    G4cout
     << G4endl
     << "--------------------End of Local Run------------------------";
  }
  
  G4cout
     << G4endl
     << " The run consists of " << nofEvents << " "<< runCondition
     << G4endl
     << "------------------------------------------------------------"
     << G4endl
     << G4endl;

  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  analysis->Write();
  analysis->CloseFile();
  G4cout << "Ntuples saved !" << G4endl;
}

void B1RunAction::AddEdep(G4int strip, G4double edep)
{
  StripVector.push_back(strip);
  EdepVector.push_back(edep);
}
