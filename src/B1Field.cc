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

#include "G4SystemOfUnits.hh"
#include "G4GeometryManager.hh"
#include "G4Track.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1Field::B1Field() : G4ElectroMagneticField(),
    Bx1(0), By1(0), Bz1(0), Ex1(0), Ey1(0), Ez1(-9.0E+5*volt/m),
    Bx2(0), By2(0), Bz2(0), Ex2(0), Ey2(0), Ez2(9.0E+5*volt/m)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1Field::~B1Field()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1Field::GetFieldValue( const G4double Point[4], G4double* Bfield ) const
{
  // Point[0],Point[1],Point[2] are x-, y-, z-cordinates, Point[3] is time

  if(Point[2] >= 0)
  {
      Bfield[0]=Bx1;
      Bfield[1]=By1;
      Bfield[2]=Bz1;
      Bfield[3]=Ex1;
      Bfield[4]=Ey1;
      Bfield[5]=Ez1;
  } else
  {
      Bfield[0]=Bx2;
      Bfield[1]=By2;
      Bfield[2]=Bz2;
      Bfield[3]=Ex2;
      Bfield[4]=Ey2;
      Bfield[5]=Ez2;
  }


  return;
}
