//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B1SteppingAction.cc
/// \brief Implementation of the B1SteppingAction class

#include "B1SteppingAction.hh"
#include "B1EventAction.hh"
#include "B1DetectorConstruction.hh"
#include "HistoManager.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::B1SteppingAction(B1EventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  posAB(-1.)
  //fScoringVolume(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::~B1SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1SteppingAction::UserSteppingAction(const G4Step* step)
{
  /*if (!fScoringVolume) {
    const B1DetectorConstruction* detectorConstruction
      = static_cast<const B1DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detectorConstruction->GetScoringVolume();   
  }

  // get volume of the current step
  G4LogicalVolume* volume 
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();
      
  // check if we are in scoring volume
  if (volume != fScoringVolume) return;

  // collect energy deposited in this step
  G4double edepStep = step->GetTotalEnergyDeposit();
  //fEventAction->AddEdep(edepStep);*/
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  if(posAB == -1.){
      const B1DetectorConstruction* detectorConstruction
            = static_cast<const B1DetectorConstruction*>
              (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
      posAB = detectorConstruction->GetPosAB();
  }
  if(step->GetTrack()->GetTrackID() == 1){
      G4ThreeVector pre = step->GetPreStepPoint()->GetPosition();
      G4ThreeVector post = step->GetPostStepPoint()->GetPosition();
      G4double pre_z = pre.getZ();
      G4double post_z = post.getZ();

      if(pre_z < -posAB && post_z > -posAB){
          //G4cout << "PreA = " << pre_z << " and PostA = " << post_z << G4endl;
          G4double pre_x = pre.getX();
          G4double post_x = post.getX();
          G4double pre_y = pre.getY();
          G4double post_y = post.getY();
          G4double coef = (-posAB - pre_z) / (post_z - pre_z);
          G4double x = pre_x + coef * (post_x - pre_x);
          G4double y = pre_y + coef * (post_y - pre_y);

          analysisManager->FillNtupleDColumn(8, x);
          analysisManager->FillNtupleDColumn(9, y);
          analysisManager->FillNtupleDColumn(10, -posAB);
          //G4cout << "x_A " << x << "  y_A " << y << "  z_A " << -posAB << G4endl;
      }
      if(pre_z < 0. && post_z > 0.){
          //G4cout << "Pre0 = " << pre_z << " and Post0 = " << post_z << G4endl;
          G4double pre_x = pre.getX();
          G4double post_x = post.getX();
          G4double pre_y = pre.getY();
          G4double post_y = post.getY();
          G4double coef = - pre_z / (post_z - pre_z);
          G4double x = pre_x + coef * (post_x - pre_x);
          G4double y = pre_y + coef * (post_y - pre_y);

          analysisManager->FillNtupleDColumn(6, x);
          analysisManager->FillNtupleDColumn(7, y);
          //G4cout << "x_0 " << x << "  y_0 " << y << G4endl;
      }
      if(pre_z < posAB && post_z > posAB){
          //G4cout << "PreB = " << pre_z << " and PostB = " << post_z << G4endl;
          G4double pre_x = pre.getX();
          G4double post_x = post.getX();
          G4double pre_y = pre.getY();
          G4double post_y = post.getY();
          G4double coef = (posAB - pre_z) / (post_z - pre_z);
          G4double x = pre_x + coef * (post_x - pre_x);
          G4double y = pre_y + coef * (post_y - pre_y);
          //G4cout << "x_B " << x << "  y_B " << y << "  z_B " << posAB << G4endl;

          analysisManager->FillNtupleDColumn(11, x);
          analysisManager->FillNtupleDColumn(12, y);
          analysisManager->FillNtupleDColumn(13, posAB);
      }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

