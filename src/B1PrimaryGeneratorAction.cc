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
  posAB(0),
  distance(0),
  rayon(0),
  pTMomentum(0),
  masse(0),
  c_speed(2.99792458e+8)
{
  G4int n_particle = 1;
  G4String particleName = "e-";
  pTMomentum = 2.e+9; // simulated transverse momentum
  G4double B = 3.8;
  distance = 0.6*m;
  // compute radius
  rayon = (pTMomentum / (c_speed * B))*m;

  // Create particle Gun
  fParticleGun  = new G4ParticleGun(n_particle);
  // Particle type selection
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName="e-");
  fParticleGun->SetParticleDefinition(particle);

  if(particle->GetParticleName() == "e-"){
      masse = 0.510998910*MeV;
  }else if (particle->GetParticleName() == "proton"){
      masse = 938.272013*MeV;
  }else {
      G4String excep("Waring! No momentum computaiton. Particle mass unknown! See B1PrimiaryGeneratorAction.cc");
      G4Exception("B1PrimaryGeneratorAction::B1PrimaryGeneratorAction()", "InvalidSetup", FatalException, excep);
  }
}



B1PrimaryGeneratorAction::~B1PrimaryGeneratorAction()
{
  delete fParticleGun;
}


//this function is called at the begining of ecah event
void B1PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // O point selection
  G4double x0 = 90*um; //360*um * (G4UniformRand()-0.5)
  G4double y0 = 90*um * (G4UniformRand()-0.5);
  // theta computation
  G4double D = sqrt(pow(x0, 2) + pow(distance, 2));
  G4double beta = asin(D / (2 * rayon));
  // phi computation
  G4double phi = atan(y0 / (2 * beta * rayon));
  // total momentum computation
  G4double momemtum = pTMomentum / cos(phi);
  // kinetic Energy computation [eV]
  G4double E_tot = (sqrt(pow(momemtum, 2) + pow(masse, 2)))*1e-6;
  G4double E_kin = (E_tot - masse);
  fParticleGun->SetParticleEnergy(E_kin);
  // circle equation computation (x;z) plane
  G4double A = (pow(x0, 2) / pow(distance, 2)) + 1;
  G4double B = -((pow(x0, 3) / pow(distance, 2)) + x0);
  G4double C = (pow(x0, 4) / (4*pow(distance, 2))) + 0.5*pow(x0, 2) + 0.25*pow(distance, 2) - pow(rayon, 2);
  G4double Delta = pow(B, 2) - 4*A*C;
  G4double a = 0;
  if(fParticleGun->GetParticleDefinition()->GetParticleName() == "e-"){
      a = (-B + sqrt(Delta)) / (2*A);
  } else if(fParticleGun->GetParticleDefinition()->GetParticleName() == "proton"){
      a = (-B - sqrt(Delta)) / (2*A);
  } else {
      G4String excep("Waring! No momentum computaiton. Particle type unknown! See B1PrimiaryGeneratorAction.cc");
      G4Exception("B1PrimaryGeneratorAction::GeneratePrimaries()", "InvalidSetup", FatalException, excep);
  }
  G4double b = sqrt(pow(rayon, 2) - pow(a, 2)) - distance;
  if (space == 0) {
    const B1DetectorConstruction* detectorConstruction
          = static_cast<const B1DetectorConstruction*>
            (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    space = detectorConstruction->GetSpace();
    posAB = detectorConstruction->GetPosAB();
  }
  // computation of x_I, x_A and x_B from z_I, z_A and z_B
  G4double z_I = -space;
  G4double z_A = -posAB;
  G4double z_B = posAB;
  G4double x_I = -sqrt(pow(rayon, 2) - pow((z_I - b), 2)) + a;
  G4double x_A = -sqrt(pow(rayon, 2) - pow((z_A - b), 2)) + a;
  G4double x_B = -sqrt(pow(rayon, 2) - pow((z_B - b), 2)) + a;
  // computation of y_I, y_A and y_B
  // y_I
  G4double D_I = sqrt(pow(x_I, 2) + pow((distance + z_I), 2));
  G4double beta_I = asin(D_I / (2 * rayon));
  G4double alpha_I = 2 * beta_I;
  G4double y_I = alpha_I * rayon * tan(phi);
  // y_A
  G4double D_A = sqrt(pow(x_A, 2) + pow((distance + z_A), 2));
  G4double beta_A = asin(D_A / (2 * rayon));
  G4double alpha_A = 2 * beta_A;
  G4double y_A = alpha_A * rayon * tan(phi);
  // y_B
  G4double D_B = sqrt(pow(x_B, 2) + pow((distance + z_B), 2));
  G4double beta_B = asin(D_B / (2 * rayon));
  G4double alpha_B = 2 * beta_B;
  G4double y_B = alpha_B * rayon * tan(phi);
  // theta_I (initial) computation
  G4double m_dir = (z_I - b) / (x_I - a);
  G4double theta_I = atan(abs(m_dir));

  fParticleGun->SetParticlePosition(G4ThreeVector(x_I, y_I, z_I));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(theta_I,phi,1.));

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  analysisManager->FillNtupleDColumn(0, x_I);
  analysisManager->FillNtupleDColumn(1, y_I);
  analysisManager->FillNtupleDColumn(2, z_I);
  analysisManager->FillNtupleDColumn(3, theta_I);
  analysisManager->FillNtupleDColumn(4, phi);
  analysisManager->FillNtupleDColumn(5, E_kin);
  analysisManager->FillNtupleDColumn(6, x0);
  analysisManager->FillNtupleDColumn(7, y0);
  analysisManager->FillNtupleDColumn(8, x_A);
  analysisManager->FillNtupleDColumn(9, y_A);
  analysisManager->FillNtupleDColumn(10, z_A);
  analysisManager->FillNtupleDColumn(11, x_B);
  analysisManager->FillNtupleDColumn(12, y_B);
  analysisManager->FillNtupleDColumn(13, z_B);

  fParticleGun->GeneratePrimaryVertex(anEvent);
}


