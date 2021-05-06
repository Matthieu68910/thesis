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
/// \file B1DetectorConstruction.hh
/// \brief Definition of the B1DetectorConstruction class

#ifndef B1DetectorConstruction_h
#define B1DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "tls.hh"
#include "G4FieldManager.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4UniformMagField;
class G4ChordFinder;
class G4GenericMessenger;
class B1Field;
class B1WorldField;

/// Detector construction class to define materials and geometry.

class B1DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    B1DetectorConstruction();
    virtual ~B1DetectorConstruction();

  public:
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();
    
    G4double    GetSpace() const { return space; }
    G4double    GetPosAB() const { return posAB; }
    G4double    GetBackPlane() const { return backplane; }
    G4int       GetStripNumber() const {return strip_nbr; }

  protected:
    G4double    space;
    G4double    posAB;
    G4int       strip_nbr;
    G4double    backplane;

  private:
    G4LogicalVolume* logicSiUp;
    G4LogicalVolume* logicSiDo;
    G4LogicalVolume* logicWorld;

    static G4ThreadLocal B1Field* fFieldUp;
    static G4ThreadLocal B1Field* fFieldDo;
    static G4ThreadLocal B1WorldField* fFieldWorld;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

