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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "B1PrimaryGeneratorMessenger.hh"

#include "B1PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithABool.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorMessenger::B1PrimaryGeneratorMessenger(
                                             B1PrimaryGeneratorAction* Gun)
:G4UImessenger(),fAction(Gun),
 fGunDir(0),
 fpTCmd(0),
 fRndmCmd(0),
 fSet_min(0),
 fSet_max(0)
{
  fGunDir = new G4UIdirectory("/generator/momentum/");
  fGunDir->SetGuidance("Primary Generator Control");

  fRndmCmd = new G4UIcmdWithABool("/generator/momentum/SetRandom",this);
  fRndmCmd->SetGuidance("Random pT generation");
  fRndmCmd->SetGuidance("true or fasle");
  fRndmCmd->SetParameterName("Logic",false);
  fRndmCmd->SetDefaultValue(true);
  fRndmCmd->AvailableForStates(G4State_Idle);
  
  fpTCmd = new G4UIcmdWithADouble("/generator/momentum/SetMomentum",this);
  fpTCmd->SetGuidance("Select simulated pT");
  fpTCmd->SetGuidance("unit is GeV (max 1000)");
  fpTCmd->SetParameterName("pTMomentum",false);
  fpTCmd->SetRange("pTMomentum>=0.4&&pTMomentum<=1000.");
  fpTCmd->AvailableForStates(G4State_Idle);

  fSet_min = new G4UIcmdWithADouble("/generator/momentum/SetRandomMinValue",this);
  fSet_min->SetGuidance("Set min randum value");
  fSet_min->SetGuidance("unit is GeV (min 0.4)");
  fSet_min->SetParameterName("min_value",false);
  fSet_min->SetRange("min_value>=0.4&&min_value<=1000.");
  fSet_min->AvailableForStates(G4State_Idle);

  fSet_max = new G4UIcmdWithADouble("/generator/momentum/SetRandomMaxValue",this);
  fSet_max->SetGuidance("Set max randum value");
  fSet_max->SetGuidance("unit is GeV (max 1000)");
  fSet_max->SetParameterName("max_value",false);
  fSet_max->SetRange("max_value>=0.5&&max_value<=1000.");
  fSet_max->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorMessenger::~B1PrimaryGeneratorMessenger()
{
  delete fpTCmd;
  delete fGunDir;
  delete fRndmCmd;
  delete fSet_min;
  delete fSet_max;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command,
                                               G4String newValue)
{
  if( command == fpTCmd )
  {
      fAction->SetTransverseMomentum(fpTCmd->GetNewDoubleValue(newValue));
  }
  if( command == fRndmCmd )
  {
      fAction->SetTMRandom(fRndmCmd->GetNewBoolValue(newValue));
  }
  if( command == fSet_min )
  {
      fAction->SetMinValue(fSet_min->GetNewDoubleValue(newValue));
  }
  if( command == fSet_max )
  {
      fAction->SetMaxValue(fSet_max->GetNewDoubleValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

