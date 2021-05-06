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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "B1Field.hh"
#include "B1DetectorConstruction.hh"

#include "G4SystemOfUnits.hh"
#include "G4GeometryManager.hh"
#include "G4Track.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1Field::B1Field() : G4ElectroMagneticField(),
    backplane(0)
{
    if(backplane == 0.){
        const B1DetectorConstruction* detectorConstruction
              = static_cast<const B1DetectorConstruction*>
                (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        backplane = detectorConstruction->GetBackPlane();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1Field::~B1Field()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1Field::GetFieldValue( const G4double Point[4], G4double* Bfield ) const
{
  // Point[0],Point[1],Point[2] are x-, y-, z-cordinates, Point[3] is time



  if(Point[2] >= 0)     // Ez1(-9.0E+5*volt/m)
  {
      Bfield[0]=0.;
      Bfield[1]=0.;
      Bfield[2]=0.;
      Bfield[3]=0.;
      Bfield[4]=0.;
      Bfield[5]= - (4.6394E+9 * ((Point[2] - backplane + 0.0585) / 1000)) * volt/m;
      //G4cout << "    calc = " << G4BestUnit((Point[2] - backplane + 0.0585), "Length") << G4endl;
  } else                // Ez1(9.0E+5*volt/m)
  {
      Bfield[0]=0.;
      Bfield[1]=0.;
      Bfield[2]=0.;
      Bfield[3]=0.;
      Bfield[4]=0.;
      Bfield[5]= (4.6394E+9 * ((0.0585 - (Point[2] + backplane)) / 1000)) * volt/m;
      //G4cout << "    calc = " << G4BestUnit((0.0585 - (Point[2] + backplane)), "Length") << G4endl;
  }

  //G4cout << "z = " << G4BestUnit(Point[2], "Length") << " Field = " << G4BestUnit(Bfield[5], "Electric field") << G4endl;
  //G4cout << "    bkpl = " << G4BestUnit(backplane, "Length") << G4endl;

  return;
}
