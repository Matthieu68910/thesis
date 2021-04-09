/// \file B1PrimaryGeneratorAction.cc
/// \brief Implementation of the B1PrimaryGeneratorAction class

#include "B1PrimaryGeneratorAction.hh"
#include "B1DetectorConstruction.hh"
#include "HistoManager.hh"

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
  fEnvelopeBox(0),
  space(0),
  theta_i(0)
{
  // variables
  G4double pTMomentum = 1.87568; // simulated transverse momentum [GeV]
  G4double distance = 0.6; // [m]

  // compute theta_i
  theta_i = asin((0.57 * distance) / pTMomentum);  // radian

  // Create particle Gun
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName="e+");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(tan(theta_i),0,1.));
  fParticleGun->SetParticleEnergy(5.*GeV);
}



B1PrimaryGeneratorAction::~B1PrimaryGeneratorAction()
{
  delete fParticleGun;
}


//this function is called at the begining of ecah event
void B1PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // 0 point selection
  G4double x0 = 360*um * (G4UniformRand()-1); //360*um * (G4UniformRand()-0.5)
  G4double y0 = 360*um * (G4UniformRand()-0.5);
  if (space == 0) {
    const B1DetectorConstruction* detectorConstruction
          = static_cast<const B1DetectorConstruction*>
            (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    space = detectorConstruction->GetSpace();
  }

  fParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, -space));

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  analysisManager->FillNtupleDColumn(0, x0);
  analysisManager->FillNtupleDColumn(1, y0);
  analysisManager->FillNtupleDColumn(2, -space);
  analysisManager->FillNtupleDColumn(3, theta_i);

  /*G4cout << "************* computations *************" << G4endl;
  G4cout << "x_0 " << x0 << "  y_0 " << y0 << G4endl;
  G4cout << "x_A " << x_A << "  y_A " << y_A << "  z_A " << z_A << G4endl;
  G4cout << "x_B " << x_B << "  y_B " << y_B << "  z_B " << z_B << G4endl;
  G4cout << "****************************************" << G4endl;*/


  fParticleGun->GeneratePrimaryVertex(anEvent);
}


