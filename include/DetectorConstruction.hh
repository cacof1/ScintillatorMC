#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "OrganicMaterial.hh"
#include "pCTconfig.hh"
#include "SensitiveDetectorTracker.hh"
#include "SensitiveDetectorScintillator.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class OrganicMaterial;
class G4ProductionCuts;
class G4PhantomParameterisation;
class G4PVParameterised;
class TH3S;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction();//G4String,G4double,G4double,G4String);
  ~DetectorConstruction();
  virtual G4VPhysicalVolume* Construct();
  static inline DetectorConstruction* GetInstance() { return theDetector; }
  G4String thePhantom;
  G4String theCTFileName;
  G4String theScintillator;
  G4double theThickness;
  OrganicMaterial* theMaterial;

  G4int theAngle;
  G4double PhantomHalfX,PhantomHalfY,PhantomHalfZ;
  G4double VoxelHalfX,VoxelHalfY,VoxelHalfZ;
  G4int NbinsX,NbinsY,NbinsZ;
  G4double Xmin, Ymin, Zmin;
  G4double Xmax, Ymax, Zmax;
  G4double midX,midY,midZ ;
  G4ThreeVector shift;
  TH3S* hu;
  G4double ScintHalfX,ScintHalfY,ScintHalfZ;
  G4double ScintPosX,ScintPosY,ScintPosZ;
  G4Material* water;
  G4ProductionCuts*  fCuts;

  pCTconfig* theConfig;
  //XCAT
  G4PhantomParameterisation* param;
  G4PVParameterised* phantomPhys;
  vector<G4Material*> theMaterialList;
  map<G4int,G4int> hu2id;
  vector<G4int> huList;
  map<G4int,G4double> hu2density;


private:
  static DetectorConstruction* theDetector;

};


#endif
