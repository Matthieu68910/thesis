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
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::B1SteppingAction(B1EventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::~B1SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1SteppingAction::UserSteppingAction(const G4Step* step)
{
      // get volume of the current step
      G4LogicalVolume* volume
        = step->GetPreStepPoint()->GetTouchableHandle()
          ->GetVolume()->GetLogicalVolume();

      if(volume->GetName() == "logicStrip"){
          // get strip number
          G4int stripNbr = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetCopyNo();
          // get energy deposited in this step
          G4double edepStep = step->GetTotalEnergyDeposit();
          // get position
          G4double pos_x = step->GetPreStepPoint()->GetPosition().getX();
          //G4cout << "pos_x = " << pos_x << G4endl;
          G4double pos_strip;
          if(stripNbr <= 253){
              pos_strip = (stripNbr - 127 + 0.5) * 0.09*mm;
          }else{
              pos_strip = (stripNbr - 381 + 0.5) * 0.09*mm;
          }
          //G4cout << "pos_strip = " << pos_strip << G4endl;
          G4double difference = pos_x - pos_strip; // < +/- 45*um
          //G4cout << "difference = " << difference << G4endl;
          if(difference > 0){
              G4double factor = pow((difference / (45.*um)), 20) * 0.5;
              //G4cout << "Factor+ = " << factor << G4endl;
              fEventAction->AddEnergy(stripNbr, (edepStep * (1-factor)));
              if(stripNbr < 507){fEventAction->AddEnergy(stripNbr+1, (edepStep * factor));}
          }else if(difference < 0){
              G4double factor = pow((-difference / (45.*um)), 20) * 0.5;
              //G4cout << "Factor- = " << factor << G4endl;
              fEventAction->AddEnergy(stripNbr, (edepStep * (1-factor)));
              if(stripNbr > 0){fEventAction->AddEnergy(stripNbr-1, (edepStep * factor));}
          }
          //G4cout << G4endl;


          /*
          G4cout << "Volume name: "
                 << volume->GetName()
                 << " number "
                 << stripNbr
                 << ", at: "
                 << G4BestUnit(pos_x, "Length")
                 << ", energy: "
                 << G4BestUnit(edepStep, "Energy")
                 << G4endl;
         */
      }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

