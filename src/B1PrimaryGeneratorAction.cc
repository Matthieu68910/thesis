/// \file B1PrimaryGeneratorAction.cc
/// \brief Implementation of the B1PrimaryGeneratorAction class

#include "B1PrimaryGeneratorAction.hh"

#include "B1PrimaryGeneratorMessenger.hh"
#include "B1DetectorConstruction.hh"
#include "HistoManager.hh"

#include "G4Event.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "G4Exception.hh"
#include <math.h>

B1PrimaryGeneratorAction::B1PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0),
  space(0),
  pTMomentum(0.),
  RdmPT(false),
  fGunMessenger(0),
  min_pT(0.),
  max_pT(0.)
{
  // Create particle Gun
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName="pi+"); // pi+ 120, proton 120, e+ 5
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
  fParticleGun->SetParticleEnergy(120.*GeV);

  //create a messenger for this class
  fGunMessenger = new B1PrimaryGeneratorMessenger(this);
}



B1PrimaryGeneratorAction::~B1PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fGunMessenger;
}

//this function is called at the begining of ecah event
void B1PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // pTMomentum and angle
  G4double theta_i;
  G4double distance = 0.6; // [m]

  if(RdmPT)
  {
      pTMomentum = min_pT + (G4UniformRand() * (max_pT - min_pT));
      theta_i = asin((0.57 * distance) / pTMomentum);  // radian
      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(tan(theta_i),0,1.));
  } else if (pTMomentum > 0.)
  {
      theta_i = asin((0.57 * distance) / pTMomentum);  // radian
      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(tan(theta_i),0,1.));
  } else
  {
      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
      //G4cout << "!!! WARNING !!! No pT input. Set to défaut" << G4endl;
  }

  // 0 point selection
  if (space == 0) {
    const B1DetectorConstruction* detectorConstruction
          = static_cast<const B1DetectorConstruction*>
            (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    space = detectorConstruction->GetSpace();
  }


  G4double x0 = 900*um * (G4UniformRand()-0.5) - (tan(theta_i) * space);
  G4double y0 = 900*um * (G4UniformRand()-0.5);

  fParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, -(space + 1.*um)));

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  analysisManager->FillNtupleDColumn(0, x0);
  analysisManager->FillNtupleDColumn(1, y0);
  analysisManager->FillNtupleDColumn(2, -(space + 1.*um));
  analysisManager->FillNtupleDColumn(3, theta_i);
  analysisManager->FillNtupleDColumn(4, pTMomentum);

  fParticleGun->GeneratePrimaryVertex(anEvent);
}
