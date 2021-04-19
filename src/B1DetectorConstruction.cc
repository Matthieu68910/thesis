/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"
#include "B1Field.hh"
#include "B1WorldField.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4Color.hh"
#include "G4VisAttributes.hh"
#include "G4UnionSolid.hh"
#include "G4PVDivision.hh"
#include "G4UnitsTable.hh"

#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ClassicalRK4.hh"
#include "G4DormandPrince745.hh"
#include "G4BogackiShampine45.hh"
#include "G4PropagatorInField.hh"
#include "G4UniformMagField.hh"
#include "G4ChordFinder.hh"
#include "G4ElectroMagneticField.hh"
#include "G4EqMagElectricField.hh"

#include "G4MultiFunctionalDetector.hh"
#include "G4SDManager.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  space(0),
  posAB(0),
  strip_nbr(254),
  logicSiUp(nullptr),
  logicSiDo(nullptr),
  logicWorld(nullptr)
{ }


B1DetectorConstruction::~B1DetectorConstruction()
{ }


G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
   
  // Option to switch on/off checking of volumes overlaps
  G4bool checkOverlaps = true;

  // Elements and Materials

  /*G4double A, Z, nel, natoms, density;
  G4String name;

  G4Element* elSi  = new G4Element("Silicium",  "Si",  Z=14., A= 28.085*g/mole);
  G4Element* elO  = new G4Element("Oxygene",  "O",  Z=8., A= 15.9994*g/mole);
  G4Element* elAl = new G4Element("Aluminium",  "Al", Z=13., A= 26.981539*g/mole);
  G4Element* elC = new G4Element("Carbon",  "C", Z=6., A= 12.0106*g/mole);
  G4Element* elH  = new G4Element("Hydrogene",  "H",  Z=1., A= 1.0000*g/mole);*/

  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR"); // G4_Galactic
  G4Material* Silicium = nist->FindOrBuildMaterial("G4_Si");
  G4Material* Aluminium = nist->FindOrBuildMaterial("G4_Al");

  // Colors

  /*G4Color red(1.0, 0.0, 0.0);
  G4VisAttributes* attrired = new G4VisAttributes(red);
  //attrired -> SetForceWireframe(true);*/
  
  G4Color green(0.0, 1.0, 0.0);
  G4VisAttributes* attrigreen = new G4VisAttributes(green);

  G4Color blue(0.0, 0.0, 1.0);
  G4VisAttributes* attriblue = new G4VisAttributes(blue);

  /*G4Color purple(0.6, 0.3, 0.9);
  G4VisAttributes* attripurple = new G4VisAttributes(purple);
  //attripurple -> SetForceWireframe(true);*/

  G4Color orange(0.95, 0.6, 0.05);
  G4VisAttributes* attriorange = new G4VisAttributes(orange);

  //************************************************************//
  //********************** VARIABLES ***************************//
  //************************************************************//

  // strips
  G4double strip_width = 90.*um;
  G4double strip_thickness = 270.*um;
  strip_nbr = 254; // 1016 for full sensor !!! change in line 44 !!!
  G4double strip_length = 5*cm;// strip_nbr * strip_width;

  // silicon backplane
  G4double Si_bp_thickness = 30.*um;

  // aluminum backplane
  G4double Al_bp_thickness = 1.*um;

  // sensors separation and tilt
  posAB = 1.375*mm; // from mid-plane to center 1.375 * 2 = 2.75 mm (375 -> 370 -> 365)
  G4double sensor_sep = 2*posAB; //2*posAB - Si_bp_thickness - Al_bp_thickness

  G4RotationMatrix tilt1  = G4RotationMatrix();
  tilt1.rotateX(0*deg); // tilt first (viewed from particle)
  tilt1.rotateY(0*deg);
  tilt1.rotateZ(0*deg);

  G4RotationMatrix tilt2  = G4RotationMatrix();
  tilt2.rotateX(0*deg); // tilt second (viewed from particle)
  tilt2.rotateY(0*deg);
  tilt2.rotateZ(0*deg);

  G4double sensor_thickness = strip_thickness + Si_bp_thickness + Al_bp_thickness;

  space = (sensor_sep + sensor_thickness) / 2;

  // world safety margin
  G4double WorldSafeX = 2*sensor_sep;
  G4double WorldSafeY = 2*sensor_sep;
  G4double WorldSafeZ = 2*sensor_sep;

  // Champ magnétique : see bellow (bottom of this file)

  //************************************************************//
  //*********************** GEOMETRY ***************************//
  //************************************************************//

  //************************ World *****************************//

  G4double world_sizeX = ((strip_nbr * strip_width) + WorldSafeX) / 2;
  G4double world_sizeY = (strip_length + WorldSafeY) / 2;
  G4double world_sizeZ = (sensor_thickness + sensor_sep + WorldSafeZ) / 2;

  G4cout
          << G4endl
          << "World dimensions are : X = "
          << 2*world_sizeX
          << " mm, Y = "
          << 2*world_sizeY
          << " mm and Z = "
          << 2*world_sizeZ
          << " mm."
          << G4endl;
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       world_sizeX, world_sizeY, world_sizeZ);     //its size
      
  logicWorld =
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking

  //************************ Sensor Up *****************************//

  // Al Up

  G4double Al_sizeX = (strip_nbr * strip_width) / 2;
  G4double Al_sizeY = strip_length / 2;
  G4double Al_sizeZ = Al_bp_thickness / 2;

  G4Box* solidAlUp =
    new G4Box("AlUp",                    //its name
        Al_sizeX, Al_sizeY, Al_sizeZ); //its size

  G4LogicalVolume* logicAlUp =
    new G4LogicalVolume(solidAlUp,            //its solid
                        Aluminium,             //its material
                        "AlUp");         //its name

  logicAlUp -> SetVisAttributes(attriblue);

  G4double AlUp_posZ = (sensor_sep / 2) - (sensor_thickness / 2) + (Al_bp_thickness / 2);
  G4cout << "Aluminium Up pos Z = " << G4BestUnit(AlUp_posZ, "Length") << G4endl;

  G4ThreeVector posAlUp = G4ThreeVector(0, 0, AlUp_posZ);
  G4Transform3D transformAlUp = G4Transform3D(tilt2,posAlUp);
  new G4PVPlacement(transformAlUp,         //at (0,0,0)
                    logicAlUp,                //its logical volume
                    "AlUp",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  // SiBp Up

  G4double SiBp_sizeX = (strip_nbr * strip_width) / 2;
  G4double SiBp_sizeY = strip_length / 2;
  G4double SiBp_sizeZ = Si_bp_thickness / 2;

  G4Box* solidSiBpUp =
    new G4Box("SiBpUp",                    //its name
        SiBp_sizeX, SiBp_sizeY, SiBp_sizeZ); //its size

  G4LogicalVolume* logicSiBpUp =
    new G4LogicalVolume(solidSiBpUp,            //its solid
                        Silicium,             //its material
                        "SiBpUp");         //its name

  logicSiBpUp -> SetVisAttributes(attriorange);

  G4double SiBpUp_posZ = AlUp_posZ + (Al_bp_thickness / 2) + (Si_bp_thickness / 2);
  G4cout << "Silicium backplane Up pos Z = " << G4BestUnit(SiBpUp_posZ, "Length") << G4endl;

  G4ThreeVector posSiBpUp = G4ThreeVector(0, 0, SiBpUp_posZ);
  G4Transform3D transformSiBpUp = G4Transform3D(tilt2,posSiBpUp);
  new G4PVPlacement(transformSiBpUp,         //at (0,0,0)
                    logicSiBpUp,                //its logical volume
                    "SiBpUp",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  // Si Up

  G4double Si_sizeX = (strip_nbr * strip_width) / 2;
  G4double Si_sizeY = strip_length / 2;
  G4double Si_sizeZ = strip_thickness / 2;

  G4Box* solidSiUp =
    new G4Box("SiUp",                    //its name
        Si_sizeX, Si_sizeY, Si_sizeZ); //its size

  logicSiUp =
    new G4LogicalVolume(solidSiUp,            //its solid
                        Silicium,             //its material
                        "SiUp");         //its name

  logicSiUp -> SetVisAttributes(attrigreen);

  G4double SiUp_posZ = SiBpUp_posZ + (Si_bp_thickness / 2) + (strip_thickness / 2);
  G4cout << "Silicium Up pos Z = " << G4BestUnit(SiUp_posZ, "Length") << G4endl;

  G4Region* siliconUpRegion = new G4Region("SiliconUp");
  logicSiUp->SetRegion(siliconUpRegion);
  siliconUpRegion->AddRootLogicalVolume(logicSiUp);

  // define strip

  G4Box* strip =
    new G4Box("Strip",                    //its name
        (strip_width/2), (strip_length/2), (strip_thickness/2)); //its size
  G4LogicalVolume* logicStrip =
    new G4LogicalVolume(strip,          //its solid
                        Silicium,           //its material
                        "logicStrip");        //its name

  // place strips

  for (G4int i = 0; i < strip_nbr ; i++) {
      G4double X_pos = (i * strip_width) - Si_sizeX + (strip_width / 2);
      G4double Y_pos = 0;
      G4double Z_pos = 0;
      G4ThreeVector position = G4ThreeVector(X_pos, Y_pos, Z_pos);
      G4Transform3D transform = G4Transform3D(G4RotationMatrix(), position);

      new G4PVPlacement(transform,             //rotation,position
                        logicStrip,            //its logical volume
                        "strip",             //its name
                        logicSiUp,             //its mother  volume
                        false,                 //no boolean operation
                        i,                 //copy number
                        checkOverlaps);       // checking overlaps
    }

  // place SiUp with strips

  G4ThreeVector posSiUp = G4ThreeVector(0, 0, SiUp_posZ);
  G4Transform3D transformSiUp = G4Transform3D(tilt2,posSiUp);
  new G4PVPlacement(transformSiUp,         //at (0,0,0)
                      logicSiUp,           //its logical volume
                      "SiUp",              //its name
                      logicWorld,              //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);         // checking overlaps
                

  //************************ Sensor Do *****************************//

  // Al Do

  G4LogicalVolume* logicAlDo =
    new G4LogicalVolume(solidAlUp,            //its solid
                        Aluminium,             //its material
                        "AlDo");         //its name

  logicAlDo -> SetVisAttributes(attriblue);

  G4double AlDo_posZ = (sensor_thickness / 2) - (sensor_sep / 2) - (Al_bp_thickness / 2);

  G4ThreeVector posAlDo = G4ThreeVector(0, 0, AlDo_posZ);
  G4Transform3D transformAlDo = G4Transform3D(tilt1,posAlDo);
  new G4PVPlacement(transformAlDo,         //at (0,0,0)
                    logicAlDo,                //its logical volume
                    "AlDo",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  // SiBp Do

  G4LogicalVolume* logicSiBpDo =
    new G4LogicalVolume(solidSiBpUp,            //its solid
                        Silicium,             //its material
                        "SiBpDo");         //its name

  logicSiBpDo -> SetVisAttributes(attriorange);

  G4double SiBpDo_posZ = AlDo_posZ - (Al_bp_thickness / 2) - (Si_bp_thickness / 2);

  G4ThreeVector posSiBpDo = G4ThreeVector(0, 0, SiBpDo_posZ);
  G4Transform3D transformSiBpDo = G4Transform3D(tilt1,posSiBpDo);
  new G4PVPlacement(transformSiBpDo,         //at (0,0,0)
                    logicSiBpDo,                //its logical volume
                    "SiBpDo",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  // Si Do

  logicSiDo =
    new G4LogicalVolume(solidSiUp,            //its solid
                        Silicium,             //its material
                        "SiDo");         //its name

  logicSiDo -> SetVisAttributes(attrigreen);

  G4double SiDo_posZ = SiBpDo_posZ - (Si_bp_thickness / 2) - (strip_thickness / 2);

  G4Region* siliconDoRegion = new G4Region("SiliconDo");
  logicSiDo->SetRegion(siliconDoRegion);
  siliconDoRegion->AddRootLogicalVolume(logicSiDo);

  // place strips

  for (G4int i = 0; i < strip_nbr ; i++) {
      G4double X_pos = (i * strip_width) - Si_sizeX + (strip_width / 2);
      G4double Y_pos = 0;
      G4double Z_pos = 0;
      G4ThreeVector position = G4ThreeVector(X_pos, Y_pos, Z_pos);
      G4Transform3D transform = G4Transform3D(G4RotationMatrix(), position);

      new G4PVPlacement(transform,             //rotation,position
                        logicStrip,            //its logical volume
                        "strip",             //its name
                        logicSiDo,             //its mother  volume
                        false,                 //no boolean operation
                        (strip_nbr+i),                 //copy number
                        checkOverlaps);       // checking overlaps
    }

  // place SiDo with strips

  G4ThreeVector posSiDo = G4ThreeVector(0, 0, SiDo_posZ);
  G4Transform3D transformSiDo = G4Transform3D(tilt1,posSiDo);
  new G4PVPlacement(transformSiDo,         //at (0,0,0)
                      logicSiDo,           //its logical volume
                      "SiDo",              //its name
                      logicWorld,              //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);         // checking overlaps

  return physWorld;
}

G4ThreadLocal B1Field* B1DetectorConstruction::fFieldUp = 0;
G4ThreadLocal B1Field* B1DetectorConstruction::fFieldDo = 0;
G4ThreadLocal B1WorldField* B1DetectorConstruction::fFieldWorld = 0;

void B1DetectorConstruction::ConstructSDandField()
{
  //******************* MF Detector **********************//
  // declare strips as a MultiFunctionalDetector scorer

  G4MultiFunctionalDetector* det = new G4MultiFunctionalDetector("stripDetector");
  G4SDManager::GetSDMpointer()->AddNewDetector(det);
  G4VPrimitiveScorer* primitiv;
  primitiv = new G4PSEnergyDeposit("edep");
  det->RegisterPrimitive(primitiv);
  SetSensitiveDetector("logicStrip",det);

  //****************** MAGNETIC FIELD *******************//
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  if(!fFieldWorld){

      fFieldWorld = new B1WorldField();

      G4EqMagElectricField* equation0 = new G4EqMagElectricField(fFieldWorld);

      G4FieldManager* WorldFieldManager = new G4FieldManager;
      WorldFieldManager->SetDetectorField(fFieldWorld);

      G4MagIntegratorStepper* stepper = new G4DormandPrince745(equation0,8);

      G4double minStep           = 100.*nm;

      G4ChordFinder* chordFinder =
                  new G4ChordFinder((G4MagneticField*)fFieldWorld,minStep,stepper);

      // Set accuracy parameters
      G4double deltaChord        = 100.*nm;
      chordFinder->SetDeltaChord( deltaChord );

      G4double deltaOneStep      = 10.*nm;
      WorldFieldManager->SetAccuraciesWithDeltaOneStep(deltaOneStep);

      G4double deltaIntersection = 10.*nm;
      WorldFieldManager->SetDeltaIntersection(deltaIntersection);

      G4TransportationManager* transportManager =
                         G4TransportationManager::GetTransportationManager();

      G4PropagatorInField* fieldPropagator =
                                    transportManager->GetPropagatorInField();

      G4double epsMin            = 1.0e-6;
      G4double epsMax            = 1.0e-5;

      fieldPropagator->SetMinimumEpsilonStep(epsMin);
      fieldPropagator->SetMaximumEpsilonStep(epsMax);

      WorldFieldManager->SetChordFinder(chordFinder);

      G4bool allLocal = true;
      logicWorld->SetFieldManager(WorldFieldManager, allLocal);
  }


  if(!fFieldUp){

      fFieldUp = new B1Field();

      G4EqMagElectricField* equation1 = new G4EqMagElectricField(fFieldUp);

      G4FieldManager* SiUpFieldManager = new G4FieldManager;
      SiUpFieldManager->SetDetectorField(fFieldUp);

      G4MagIntegratorStepper* stepper = new G4DormandPrince745(equation1,8);

      G4double minStep           = 100*nm;

      G4ChordFinder* chordFinder =
                  new G4ChordFinder((G4MagneticField*)fFieldUp,minStep,stepper);

      // Set accuracy parameters
      G4double deltaChord        = 100.*nm;
      chordFinder->SetDeltaChord( deltaChord );

      G4double deltaOneStep      = 10.*nm;
      SiUpFieldManager->SetAccuraciesWithDeltaOneStep(deltaOneStep);

      G4double deltaIntersection = 10.*nm;
      SiUpFieldManager->SetDeltaIntersection(deltaIntersection);

      G4TransportationManager* transportManager =
                         G4TransportationManager::GetTransportationManager();

      G4PropagatorInField* fieldPropagator =
                                    transportManager->GetPropagatorInField();

      G4double epsMin            = 1.0e-6;
      G4double epsMax            = 1.0e-5;

      fieldPropagator->SetMinimumEpsilonStep(epsMin);
      fieldPropagator->SetMaximumEpsilonStep(epsMax);

      SiUpFieldManager->SetChordFinder(chordFinder);

      G4bool allLocal = true;
      logicSiUp->SetFieldManager(SiUpFieldManager, allLocal);
  }

  if(!fFieldDo){

      fFieldDo = new B1Field();

      G4EqMagElectricField* equation2 = new G4EqMagElectricField(fFieldDo);

      G4FieldManager* SiDoFieldManager = new G4FieldManager;
      SiDoFieldManager->SetDetectorField(fFieldDo);

      G4MagIntegratorStepper* stepper = new G4DormandPrince745(equation2,8);

      G4double minStep           = 100.*nm;

      G4ChordFinder* chordFinder =
                  new G4ChordFinder((G4MagneticField*)fFieldDo,minStep,stepper);

      // Set accuracy parameters
      G4double deltaChord        = 100*nm;
      chordFinder->SetDeltaChord( deltaChord );

      G4double deltaOneStep      = 10.*nm;
      SiDoFieldManager->SetAccuraciesWithDeltaOneStep(deltaOneStep);

      G4double deltaIntersection = 10.*nm;
      SiDoFieldManager->SetDeltaIntersection(deltaIntersection);

      G4TransportationManager* transportManager =
                         G4TransportationManager::GetTransportationManager();

      G4PropagatorInField* fieldPropagator =
                                    transportManager->GetPropagatorInField();

      G4double epsMin            = 1.0e-6;
      G4double epsMax            = 1.0e-5;

      fieldPropagator->SetMinimumEpsilonStep(epsMin);
      fieldPropagator->SetMaximumEpsilonStep(epsMax);

      SiDoFieldManager->SetChordFinder(chordFinder);

      G4bool allLocal = true;
      logicSiDo->SetFieldManager(SiDoFieldManager, allLocal);
  }
  //**********************************************************//
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
