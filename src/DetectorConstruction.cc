#include "G4NistManager.hh"
#include "G4EmCalculator.hh"
#include "DetectorConstruction.hh"
#include "G4RunManager.hh"
#include "TString.h"
#include "TH3F.h"
#include "TFile.h"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4SubtractionSolid.hh"
#include "G4PhantomParameterisation.hh"
#include "G4PVParameterised.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4EllipticalTube.hh"
#include "G4Ellipsoid.hh"
#include "G4LogicalVolume.hh"
#include "G4Trap.hh"
#include "G4UnionSolid.hh"
#include "SensitiveDetectorTracker.hh"
#include "SensitiveDetectorScintillator.hh"
#include "G4PVPlacement.hh"
#include "OrganicMaterial.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UserLimits.hh"
#include "G4Region.hh"
#include "G4UserLimits.hh"

using namespace CLHEP;
DetectorConstruction* DetectorConstruction::theDetector=NULL;
DetectorConstruction::~DetectorConstruction()
{
  delete theDetector;
}

DetectorConstruction::DetectorConstruction()//G4String theModel,G4double angle,G4double thick, G4String CTFileName)
: G4VUserDetectorConstruction()
{
  theConfig = pCTconfig::GetInstance();
  theCTFileName = theConfig->item_str["CTPath"];
  theThickness  = theConfig->item_float["Thickness"];
  theAngle      = theConfig->item_float["Angle"];
  theDetector   = this;
  thePhantom    = theConfig->item_str["Model"];
  theScintillator    = theConfig->item_str["Scintillator"];
  theMaterial = OrganicMaterial::GetInstance();

  if(thePhantom=="XCAT"){
    TFile *f  = new TFile(theCTFileName,"update");
    hu        = (TH3S*)f->Get("hu");
    NbinsX =  hu->GetNbinsX();
    NbinsY =  hu->GetNbinsY();
    NbinsZ =  hu->GetNbinsZ();
    Xmax   =  hu->GetXaxis()->GetXmax();
    Ymax   =  hu->GetYaxis()->GetXmax();
    Zmax   =  hu->GetZaxis()->GetXmax();
    Xmin   =  hu->GetXaxis()->GetXmin();
    Ymin   =  hu->GetYaxis()->GetXmin();
    Zmin   =  hu->GetZaxis()->GetXmin();
    midX   =  (Xmax+Xmin)/2.;
    midY   =  (Ymax+Ymin)/2.;
    midZ   =  (Zmax+Zmin)/2.;
    VoxelHalfX    =  hu->GetXaxis()->GetBinWidth(1)/2.*cm;
    VoxelHalfY    =  hu->GetYaxis()->GetBinWidth(1)/2.*cm;
    VoxelHalfZ    =  hu->GetZaxis()->GetBinWidth(1)/2.*cm;
    PhantomHalfX  =  NbinsX*VoxelHalfX;
    PhantomHalfY  =  NbinsY*VoxelHalfY;
    PhantomHalfZ  =  NbinsZ*VoxelHalfZ;

    shift  =  G4ThreeVector(midX,midY,midZ);
    G4cout<<"Centroid shift from (0,0,0) ="<<shift<<G4endl;
    cout<<"NbinsX - Xmax - Xmin - BinSizeX/2 :"<<NbinsX<<" "<<Xmax<<" "<<Xmin<<" "<<PhantomHalfX<<endl;
    cout<<"NbinsY - Ymax - Ymin - BinSizeY/2 :"<<NbinsY<<" "<<Ymax<<" "<<Ymin<<" "<<PhantomHalfY<<endl;
    cout<<"NbinsZ - Zmax - Zmin - BinSizeZ/2 :"<<NbinsZ<<" "<<Zmax<<" "<<Zmin<<" "<<PhantomHalfZ<<endl;
  }

  else{
    midX = midY = midZ = 0;
    PhantomHalfX     = theThickness/2.*cm;
    PhantomHalfY     = 30./2.*cm;
    PhantomHalfZ     = 30./2.*cm;
  }

  ScintHalfX      = 30./2*cm;
  ScintHalfY      = 30./2*cm;
  ScintHalfZ      = 30./2*cm;
  ScintPosX       = midX + PhantomHalfX + ScintHalfX + 2*mm; // in the middle
  ScintPosY       = midY; // in the middle
  ScintPosZ       = midZ; // in the middle

  G4double cut    = 1.*mm;
  fCuts   = new G4ProductionCuts();
  fCuts->SetProductionCut(cut);
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  water = theMaterial->ConstructMaterial("Water",1.0);//water;
  G4Material* air   = theMaterial->ConstructMaterial("Air",0.0001025);
  G4Material* theWorldMaterial = air;
  G4double world_size       = 5*m;
  G4double PhantomPositionX = 0*cm;
  G4double PhantomPositionY = 0*cm;
  G4double PhantomPositionZ = 0*cm;

  G4double InRadius         = 5.5 *cm;
  G4double OutRadius        = 10.5*cm;
  G4double InsertRadius     = 1.45*cm;

  //Cubic world
  G4Box* boxWorld = new G4Box("World",world_size,world_size,world_size);
  G4LogicalVolume* logicWorld = new G4LogicalVolume(boxWorld, theWorldMaterial,"World");
  G4VPhysicalVolume* physWorld = new G4PVPlacement(0,G4ThreeVector(),logicWorld,"physWorld",0,false,0);
  G4VisAttributes* world_att = new G4VisAttributes(G4Colour(1,0,0));
  world_att->SetVisibility(true);
  logicWorld->SetVisAttributes(world_att);

  // Tracking Detectors
  SensitiveDetectorTracker* sd1        = new SensitiveDetectorTracker("FrontTracker"); // Front Tracker
  G4Box* rad_vol1                      = new G4Box("rad_vol1",1.0*mm,PhantomHalfY,PhantomHalfZ);
  G4LogicalVolume * rad_log1           = new G4LogicalVolume(rad_vol1, theWorldMaterial,"rad_log1",0,0,0);
  G4VisAttributes* sd_att = new G4VisAttributes(G4Colour(0,1,1));
  sd_att->SetVisibility(true);
  rad_log1->SetVisAttributes(sd_att);
  rad_log1->SetSensitiveDetector(sd1);
  new G4PVPlacement(0,G4ThreeVector(-1*PhantomHalfX - 1*mm+midX ,midY,midZ),"rad_phys1",rad_log1,physWorld,false,0);// 2.0 mm thick so the edge fit with the box edge

  SensitiveDetectorTracker* sd2        = new SensitiveDetectorTracker("RearTracker"); // Rear Tracker
  G4Box* rad_vol2                      = new G4Box("rad_vol2",1*mm,PhantomHalfY,PhantomHalfZ);
  G4LogicalVolume * rad_log2           = new G4LogicalVolume(rad_vol2, theWorldMaterial,"rad_log2",0,0,0);
  rad_log2->SetVisAttributes(sd_att);
  rad_log2->SetSensitiveDetector(sd2);
  new G4PVPlacement(0,G4ThreeVector(PhantomHalfX + 1*mm + midX ,midY,midZ),"rad_phys2",rad_log2,physWorld,false,0);// 2.0 mm thick so the edge fit with the box edge

  // Container box which is the mother of all the phantom
  G4Box* cont_vol = new G4Box("cont_vol",PhantomHalfX,PhantomHalfY,PhantomHalfZ);
  G4LogicalVolume* cont_log  = new G4LogicalVolume(cont_vol, theWorldMaterial,"cont_log");
  G4PVPlacement*   cont_phys = new G4PVPlacement(0,G4ThreeVector(midX,midY,midZ),"cont_phys",cont_log,physWorld,false,0);
  G4VisAttributes* cont_att  = new G4VisAttributes(G4Colour(0,1,0));
  cont_att->SetVisibility(true);
  cont_log->SetVisAttributes(cont_att);

  //Scintillator
  SensitiveDetectorScintillator* scintillatorDetector = new SensitiveDetectorScintillator("Scintillator");
  G4Material* theScintillatorMaterial     = theMaterial->ConstructMaterial(theScintillator,1.0);
  //theScintillatorMaterial->GetIonisation()->SetMeanExcitationEnergy(64.7*eV);

  G4Box* boxScintillator = new G4Box("boxScintillator",ScintHalfX,ScintHalfY,ScintHalfZ);
  G4LogicalVolume* logicScintillator = new G4LogicalVolume(boxScintillator, theScintillatorMaterial,"logScintillator");
  new G4PVPlacement(0,G4ThreeVector(ScintPosX, ScintPosY, ScintPosZ),"physScintillator", logicScintillator, physWorld,false,0);
  G4VisAttributes* scint_att = new G4VisAttributes(G4Colour(1,0,0));
  scint_att->SetVisibility(true);
  logicScintillator->SetVisAttributes(scint_att);
  logicScintillator->SetSensitiveDetector(scintillatorDetector);

  // Set the region with no mcs
  //G4Region* noMCS_region = new G4Region("noMCS_region");
  //noMCS_region->SetProductionCuts(fCuts);
  //noMCS_region->AddRootLogicalVolume(logicScintillator);
  //noMCS_region->AddRootLogicalVolume(box_log);

  // Set the region with small steps
  G4Region* smallSteps = new G4Region("smallSteps");
  smallSteps->SetUserLimits(new G4UserLimits(1*mm));
  smallSteps->AddRootLogicalVolume(logicScintillator);

  G4RotationMatrix* rotExt = new G4RotationMatrix();
  rotExt->rotateZ(theAngle*pi/180.);

  //----------------------------------------------------------------------------------------------------------------
  // Water Cylinder phantom
  //----------------------------------------------------------------------------------------------------------------
  if(thePhantom == "WaterCylinder"){
    G4VisAttributes* linepair_att  = new G4VisAttributes(G4Colour(0,0,1));
    linepair_att->SetVisibility(true);
    linepair_att->SetForceSolid(true);
    G4Tubs* WaterCylinder  = new G4Tubs("WaterCylinder",0, 20./2*cm, 30./2*cm, 0, 2*pi);
    G4LogicalVolume* WaterCyl_log = new G4LogicalVolume(WaterCylinder, water, "WaterCyl_log");
    new G4PVPlacement(rotExt,G4ThreeVector(0,0,0),"WaterCyl_phys",WaterCyl_log,cont_phys,false,0);

  }

  //----------------------------------------------------------------------------------------------------------------
  // Water Cylinder phantom
  //----------------------------------------------------------------------------------------------------------------
  else if(thePhantom == "WaterCylinder_Air"){
    G4VisAttributes* linepair_att  = new G4VisAttributes(G4Colour(0,0,1));
    linepair_att->SetVisibility(true);
    linepair_att->SetForceSolid(true);
    G4Tubs* WaterCylinder  = new G4Tubs("WaterCylinder",0, 30./2*cm, 30./2*cm, 0, 2*pi);
    G4LogicalVolume* WaterCyl_log = new G4LogicalVolume(WaterCylinder, water, "WaterCyl_log");
    G4PVPlacement* WaterCyl_phys = new G4PVPlacement(0,G4ThreeVector(0,0,0),"WaterCyl_phys",WaterCyl_log,cont_phys,false,0);

    G4Tubs* AirCylinder  = new G4Tubs("WaterCylinder",0, 2./2*cm, 30./2*cm, 0, 2*pi);
    G4LogicalVolume* AirCyl_log = new G4LogicalVolume(AirCylinder, air, "AirCyl_log");
    new G4PVPlacement(rotExt,G4ThreeVector(0,0,0),"AirCyl_phys",AirCyl_log,WaterCyl_phys,false,0);
  }

  //----------------------------------------------------------------------------------------------------------------
  // Bone box phantom
  //----------------------------------------------------------------------------------------------------------------
  else if(thePhantom == "BoneBox"){
    G4VisAttributes* linepair_att  = new G4VisAttributes(G4Colour(0,0,1));
    linepair_att->SetVisibility(true);
    linepair_att->SetForceSolid(true);
    G4Material* cb230_bone = theMaterial->ConstructMaterial("CB230",1.34);
    G4Box* BoneBox  = new G4Box("BoneBox", 2.5/2*cm, 30/2*cm, 0.5/2*cm);
    G4LogicalVolume* BoneBox_log = new G4LogicalVolume(BoneBox, cb230_bone, "BoneBox_log");
    new G4PVPlacement(rotExt,G4ThreeVector(0,0,0),"BoneBox_phys",BoneBox_log,cont_phys,false,0);
  }

  //----------------------------------------------------------------------------------------------------------------
  // Water box phantom
  //----------------------------------------------------------------------------------------------------------------
  else if(thePhantom == "WaterBox"){
    G4VisAttributes* linepair_att  = new G4VisAttributes(G4Colour(0,0,1));
    linepair_att->SetVisibility(true);
    linepair_att->SetForceSolid(true);
    G4Box* WaterBox  = new G4Box("WaterBox", 15/2*cm, 30/2*cm, 30/2*cm);
    G4LogicalVolume* WaterBox_log = new G4LogicalVolume(WaterBox, water, "WaterBox_log");
    new G4PVPlacement(rotExt,G4ThreeVector(0,0,0),"WaterBox_phys",WaterBox_log,cont_phys,false,0);
  }

  //----------------------------------------------------------------------------------------------------------------
  // Gammex phantom
  //----------------------------------------------------------------------------------------------------------------

  else if(thePhantom == "Gammex"){

    //Plastic Material
    std::string mat_outer[8]  = {"CB250","MineralBone","CorticalBone","CTsolidwater","LN300","CTsolidwater","Liver","CTsolidwater"};
    double dens_outer[8] = {1.56   , 1.145       , 1.84         , 1.015        , 0.3   , 1.015        , 1.08  , 1.015};
    std::string mat_inner[8]  = {"CB230","Water","Breast","LN450","Brain","AP6Adipose","InnerBone","CTsolidwater"};
    double dens_inner[8] = {1.34   , 1.0   , 0.99   , 0.45  , 1.045 , 0.92       , 1.12      , 1.015};

    //Loma Linda Material
    /*std::string mat_outer[8]  = {"SoftTissue_LL","CorticalBone_LL","TrabecularBone_LL","SpinalDisc_LL","BrainTissue_LL","ToothEnamel_LL","ToothDentin_LL","SinusCavities_LL"};
    double dens_outer[8] = {1.03, 1.80, 1.18, 1.10, 1.04, 2.89, 2.14, 0.00120};
    std::string mat_inner[8]  = {"SoftTissue_LL","CorticalBone_LL","TrabecularBone_LL","SpinalDisc_LL","BrainTissue_LL","ToothEnamel_LL","ToothDentin_LL","SinusCavities_LL"};
    double dens_inner[8] = {1.03, 1.80, 1.18, 1.10, 1.04, 2.89, 2.14, 0.00120};*/

    //Tissue Material
    //std::string mat_inner[8] = {"Cartilage","HumerusWholeSpecimen","Ribs2and6","FemurWholeSpecimen","Ribs10","Cranium","FemurCylindricalShaft","ConnectiveTissue"};
    //std::string mat_inner[8] = {"Adipose3","Adipose2","Adipose1","SoftTissue","Muscle3","Liver3","Skin2","LungInflated"};

    //External Phantom
    G4VSolid* ExtPhantom = new G4Tubs("ExtPhantom",0,(330./2)*mm,(50./2)*mm,0,2*pi);
    G4LogicalVolume *PhantomLog = new G4LogicalVolume(ExtPhantom,theMaterial->ConstructMaterial("CTsolidwater",1.015), "PhantomLog",0,0,0);

    // Inserts
    G4VSolid* insert = new G4Tubs("insert",0,InsertRadius,(50./2)*mm,0,2*pi);
    G4VisAttributes* ins_att  = new G4VisAttributes(G4Colour(0,1,0));
    ins_att->SetVisibility(true);
    ins_att->SetForceSolid(true);

    G4VisAttributes* ins_att2  = new G4VisAttributes(G4Colour(1,1,0));
    ins_att2->SetVisibility(true);
    ins_att2->SetForceSolid(true);

    //Inner Circle
    for (int i = 0; i<8; i++){
      G4LogicalVolume* insert_log = new G4LogicalVolume(insert,theMaterial->ConstructMaterial(mat_inner[i],dens_inner[i]),mat_inner[i],0,0,0);
      if(i==4){
	insert_log->SetVisAttributes(ins_att2);
      }
      else insert_log->SetVisAttributes(ins_att);
      G4ThreeVector placement(0,InRadius,0);
      placement.rotateZ(pi/4-pi*i/4);
      new G4PVPlacement(0,placement,insert_log,mat_inner[i],PhantomLog,false,0);
    }
    //Outer Circle
    for (int i = 0; i<8; i++){
      G4LogicalVolume* insert_log = new G4LogicalVolume(insert,theMaterial->ConstructMaterial(mat_outer[i],dens_outer[i]),mat_outer[i],0,0,0);
      insert_log->SetVisAttributes(ins_att);
      G4ThreeVector placement(0,OutRadius,0);
      placement.rotateZ(3*pi/8-pi*i/4);
      new G4PVPlacement(0,placement,insert_log,mat_outer[i],PhantomLog,false,0);
    }


    new G4PVPlacement(rotExt,G4ThreeVector(PhantomPositionX,PhantomPositionY,PhantomPositionZ),"Phantom",PhantomLog,cont_phys,false,0);
  }

  //----------------------------------------------------------------------------------------------------------------
  // Wedge phantom
  //----------------------------------------------------------------------------------------------------------------

  else if (thePhantom =="Wedge"){
    // Wedge 1 & 2

    //1
    G4Trap* wedge_1 =  new G4Trap("wedge",33*cm,33*cm,16.5*cm,0.0001*cm);
    G4LogicalVolume* wedg1_log  = new G4LogicalVolume(wedge_1,theMaterial->ConstructMaterial("SoftTissue",1.03),"wedg1_log");
    G4RotationMatrix *rot1 = new G4RotationMatrix(0,180*degree,0*degree);
    new G4PVPlacement(rot1,G4ThreeVector(4.125*cm,0,0),"wedg1_phys",wedg1_log,cont_phys,false,0);
    G4VisAttributes* wedg1_att  = new G4VisAttributes(G4Colour(0,0,1));
    wedg1_att->SetVisibility(true);
    wedg1_att->SetForceSolid(true);
    wedg1_log->SetVisAttributes(wedg1_att);

    //2
    G4Trap* wedge_2 =  new G4Trap("wedge",33*cm,33*cm,16.5*cm,0.0001*cm);
    G4LogicalVolume* wedg2_log  = new G4LogicalVolume(wedge_2,theMaterial->ConstructMaterial("LungInflated",0.26),"wedg2_log");
    G4RotationMatrix *rot2 = new G4RotationMatrix(0,180*degree,180*degree);
    new G4PVPlacement(rot2,G4ThreeVector(12.375*cm,0,0),"wedg2_phys",wedg2_log,cont_phys,false,0);
    G4VisAttributes* wedg2_att  = new G4VisAttributes(G4Colour(1,0,0));
    wedg2_att->SetVisibility(true);
    wedg2_att->SetForceSolid(true);
    wedg2_log->SetVisAttributes(wedg2_att);

    // Wedge 3 & 4

    //3
    G4Trap* _wedge3 =  new G4Trap("wedge",33*cm,16.5*cm,16.5*cm,0.0001*cm);
    G4LogicalVolume* wedg3_log  = new G4LogicalVolume(_wedge3,theMaterial->ConstructMaterial("HumerusWholeSpecimen",1.39),"wedg3_log");
    G4RotationMatrix *rot3 = new G4RotationMatrix(0,0,180*degree);
    new G4PVPlacement(rot3,G4ThreeVector(-4.125*cm,8.25*cm,0),"wedg3_phys",wedg3_log,cont_phys,false,0);
    G4VisAttributes* wedg3_att  = new G4VisAttributes(G4Colour(1,1,1));
    wedg3_att->SetVisibility(true);
    wedg3_att->SetForceSolid(true);
    wedg3_log->SetVisAttributes(wedg3_att);

    //4
    G4Trap* _wedge4 =  new G4Trap("wedge",33*cm,16.5*cm,16.5*cm,0.0001*cm);
    G4LogicalVolume* wedg4_log  = new G4LogicalVolume(_wedge4,theMaterial->ConstructMaterial("HumerusWholeSpecimen",1.39),"wedg4_log");
    G4RotationMatrix *rot4 = new G4RotationMatrix(0,0,90*degree);
    new G4PVPlacement(rot4,G4ThreeVector(-8.25*cm,-12.375*cm,0),"wedg4_phys",wedg4_log,cont_phys,false,0);
    G4VisAttributes* wedg4_att  = new G4VisAttributes(G4Colour(1,1,1));
    wedg4_att->SetVisibility(true);
    wedg4_att->SetForceSolid(true);
    wedg4_log->SetVisAttributes(wedg4_att);

    // Wedge 5 & 6

    //5
    G4Trap* _wedge5 =  new G4Trap("wedge",33*cm,16.5*cm,16.5*cm,0.0001*cm);
    G4LogicalVolume* wedg5_log  = new G4LogicalVolume(_wedge5,theMaterial->ConstructMaterial("Adipose3",0.93),"wedg5_log");
    G4RotationMatrix *rot5 = new G4RotationMatrix(0,0,270*degree);
    new G4PVPlacement(rot5,G4ThreeVector(-8.25*cm,-4.125*cm,0),"wedg5_phys",wedg5_log,cont_phys,false,0);
    G4VisAttributes* wedg5_att  = new G4VisAttributes(G4Colour(1,0,1));
    wedg5_att->SetVisibility(true);
    wedg5_att->SetForceSolid(true);
    wedg5_log->SetVisAttributes(wedg5_att);

    //6
    G4Trap* _wedge6 =  new G4Trap("wedge",33*cm,16.5*cm,16.5*cm,0.0001*cm);
    G4LogicalVolume* wedg6_log  = new G4LogicalVolume(_wedge6,theMaterial->ConstructMaterial("Adipose3",0.93),"wedg6_log");
    new G4PVPlacement(0,G4ThreeVector(-12.375*cm,8.25*cm,0),"wedg6_phys",wedg6_log,cont_phys,false,0);
    G4VisAttributes* wedg6_att  = new G4VisAttributes(G4Colour(1,0,1));
    wedg6_att->SetVisibility(true);
    wedg6_att->SetForceSolid(true);
    wedg6_log->SetVisAttributes(wedg6_att);

  }

  //----------------------------------------------------------------------------------------------------------------
  // MTF phantom (Line pair for radiography)
  //----------------------------------------------------------------------------------------------------------------
  else if (thePhantom =="MTF"){
    G4Box* WaterBox  = new G4Box("WaterBox", 15./2*cm, PhantomHalfY, PhantomHalfZ);
    G4LogicalVolume* WaterBox_log = new G4LogicalVolume(WaterBox, water, "WaterBox_log");
    G4PVPlacement* watercont_phys = new G4PVPlacement(0,G4ThreeVector(0,0,0),"Watercont_phys",WaterBox_log,cont_phys,false,0);

    G4VisAttributes* linepair_att  = new G4VisAttributes(G4Colour(0,0,1));
    linepair_att->SetVisibility(true);
    linepair_att->SetForceSolid(true);

    G4VisAttributes* linepair_att2  = new G4VisAttributes(G4Colour(0,1,0));
    linepair_att2->SetVisibility(true);
    linepair_att2->SetForceSolid(true);
    G4Material *Aluminium = new G4Material("Aluminum", 13, 26.98*g/mole, 2.7*g/cm3);

    // 1lp\cm
    G4double halfX1 = (5./1.)/2.;
    G4Box* linePairBox1 = new G4Box("LinePair",2*cm,halfX1*mm,20*cm);
    G4LogicalVolume* linePairLog1 = new G4LogicalVolume(linePairBox1,Aluminium ,"linePairLog");
    linePairLog1->SetVisAttributes(linepair_att);

    new G4PVPlacement(0,G4ThreeVector(0*cm,140*mm,0*cm),"linePairPhys",linePairLog1,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(140-4*halfX1)*mm,0*cm),"linePairPhys",linePairLog1,watercont_phys,false,0);

    // 2lp\cm
    G4double halfX2 = (5./2)/2.;
    G4Box* linePairBox2 = new G4Box("LinePair",2*cm,halfX2*mm,20*cm);
    G4LogicalVolume* linePairLog2 = new G4LogicalVolume(linePairBox2,Aluminium ,"linePairLog");
    linePairLog2->SetVisAttributes(linepair_att);

    new G4PVPlacement(0,G4ThreeVector(0*cm,120*mm,0*cm),"linePairPhys",linePairLog2,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(120-4*halfX2)*mm,0*cm),"linePairPhys",linePairLog2,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(120-8*halfX2)*mm,0*cm),"linePairPhys",linePairLog2,watercont_phys,false,0);


    // 3lp\cm
    G4double halfX3 = (5./3.)/2.;
    G4Box* linePairBox3 = new G4Box("LinePair",2*cm,halfX3*mm,20*cm);
    G4LogicalVolume* linePairLog3 = new G4LogicalVolume(linePairBox3,Aluminium ,"linePairLog");
    linePairLog3->SetVisAttributes(linepair_att);

    new G4PVPlacement(0,G4ThreeVector(0*cm,100*mm,0*cm),"linePairPhys",linePairLog3,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(100-4*halfX3)*mm,0*cm),"linePairPhys",linePairLog3,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(100-8*halfX3)*mm,0*cm),"linePairPhys",linePairLog3,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(100-12*halfX3)*mm,0*cm),"linePairPhys",linePairLog3,watercont_phys,false,0);

    // 4lp\cm
    G4double halfX4 = (5./4.)/2.;
    G4Box* linePairBox4 = new G4Box("LinePair",2*cm,halfX4*mm,20*cm);
    G4LogicalVolume* linePairLog4 = new G4LogicalVolume(linePairBox4,Aluminium ,"linePairLog");
    linePairLog4->SetVisAttributes(linepair_att);

    new G4PVPlacement(0,G4ThreeVector(0*cm,80*mm,0*cm),"linePairPhys",linePairLog4,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(80-4*halfX4)*mm,0*cm),"linePairPhys",linePairLog4,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(80-8*halfX4)*mm,0*cm),"linePairPhys",linePairLog4,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(80-12*halfX4)*mm,0*cm),"linePairPhys",linePairLog4,watercont_phys,false,0);

    //5lp\cm
    G4double halfX5 = (5./5.)/2.;
    G4Box* linePairBox5 = new G4Box("LinePair",2*cm,halfX5*mm,20*cm);
    G4LogicalVolume* linePairLog5 = new G4LogicalVolume(linePairBox5,Aluminium ,"linePairLog");
    linePairLog5->SetVisAttributes(linepair_att);

    new G4PVPlacement(0,G4ThreeVector(0*cm,60*mm,0*cm),"linePairPhys",linePairLog5,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(60-4*halfX5)*mm,0*cm),"linePairPhys",linePairLog5,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(60-8*halfX5)*mm,0*cm),"linePairPhys",linePairLog5,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(60-12*halfX5)*mm,0*cm),"linePairPhys",linePairLog5,watercont_phys,false,0);

    //6lp\cm
    G4double halfX6 = (5./6.)/2.;
    G4Box* linePairBox6 = new G4Box("LinePair",2*cm,halfX6*mm,20*cm);
    G4LogicalVolume* linePairLog6 = new G4LogicalVolume(linePairBox6,Aluminium ,"linePairLog");
    linePairLog6->SetVisAttributes(linepair_att);

    new G4PVPlacement(0,G4ThreeVector(0*cm,40*mm,0*cm),"linePairPhys",linePairLog6,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(40-4*halfX6)*mm,0*cm),"linePairPhys",linePairLog6,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(40-8*halfX6)*mm,0*cm),"linePairPhys",linePairLog6,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(40-12*halfX6)*mm,0*cm),"linePairPhys",linePairLog6,watercont_phys,false,0);

    //7lp\cm
    G4double halfX7 = (5./7.)/2.;
    G4Box* linePairBox7 = new G4Box("LinePair",2*cm,halfX7*mm,20*cm);
    G4LogicalVolume* linePairLog7 = new G4LogicalVolume(linePairBox7,Aluminium ,"linePairLog");
    linePairLog7->SetVisAttributes(linepair_att);

    new G4PVPlacement(0,G4ThreeVector(0*cm,20*mm,0*cm),"linePairPhys",linePairLog7,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(20-4*halfX7)*mm,0*cm),"linePairPhys",linePairLog7,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(20-8*halfX7)*mm,0*cm),"linePairPhys",linePairLog7,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(20-12*halfX7)*mm,0*cm),"linePairPhys",linePairLog7,watercont_phys,false,0);

    //8lp\cm
    G4double halfX8 = (5./8.)/2.;
    G4Box* linePairBox8 = new G4Box("LinePair",2*cm,halfX8*mm,20*cm);
    G4LogicalVolume* linePairLog8 = new G4LogicalVolume(linePairBox8,Aluminium ,"linePairLog");
    linePairLog8->SetVisAttributes(linepair_att);

    new G4PVPlacement(0,G4ThreeVector(0*cm,0*mm,0*cm),"linePairPhys",linePairLog8,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(0-4*halfX8)*mm,0*cm),"linePairPhys",linePairLog8,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(0-8*halfX8)*mm,0*cm),"linePairPhys",linePairLog8,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(0-12*halfX8)*mm,0*cm),"linePairPhys",linePairLog8,watercont_phys,false,0);

    //9lp\cm
    G4double halfX9 = (5./9.)/2.;
    G4Box* linePairBox9 = new G4Box("LinePair",2*cm,halfX9*mm,20*cm);
    G4LogicalVolume* linePairLog9 = new G4LogicalVolume(linePairBox9,Aluminium ,"linePairLog");
    linePairLog9->SetVisAttributes(linepair_att);

    new G4PVPlacement(0,G4ThreeVector(0*cm,-20*mm,0*cm),"linePairPhys",linePairLog9,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-20-4*halfX9)*mm,0*cm),"linePairPhys",linePairLog9,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-20-8*halfX9)*mm,0*cm),"linePairPhys",linePairLog9,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-20-12*halfX9)*mm,0*cm),"linePairPhys",linePairLog9,watercont_phys,false,0);

    //10lp\cm
    G4double halfX10 = (5./10.)/2.;
    G4Box* linePairBox10 = new G4Box("LinePair",2*cm,halfX10*mm,20*cm);
    G4LogicalVolume* linePairLog10 = new G4LogicalVolume(linePairBox10,Aluminium ,"linePairLog");
    linePairLog10->SetVisAttributes(linepair_att2);

    new G4PVPlacement(0,G4ThreeVector(0*cm,(-40 -0*halfX10)*mm,0*cm),"linePairPhys",linePairLog10,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-40 -4*halfX10)*mm,0*cm),"linePairPhys",linePairLog10,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-40 -8*halfX10)*mm,0*cm),"linePairPhys",linePairLog10,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-40-12*halfX10)*mm,0*cm),"linePairPhys",linePairLog10,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-40-16*halfX10)*mm,0*cm),"linePairPhys",linePairLog10,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-40-20*halfX10)*mm,0*cm),"linePairPhys",linePairLog10,watercont_phys,false,0);

    //11lp\cm
    G4double halfX11 = (5./11.)/2.;
    G4Box* linePairBox11 = new G4Box("LinePair",2*cm,halfX11*mm,20*cm);
    G4LogicalVolume* linePairLog11 = new G4LogicalVolume(linePairBox11,Aluminium ,"linePairLog");
    linePairLog11->SetVisAttributes(linepair_att);

    new G4PVPlacement(0,G4ThreeVector(0*cm,-60*mm,0*cm),"linePairPhys",linePairLog11,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-60-4*halfX11)*mm,0*cm),"linePairPhys",linePairLog11,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-60-8*halfX11)*mm,0*cm),"linePairPhys",linePairLog11,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-60-12*halfX11)*mm,0*cm),"linePairPhys",linePairLog11,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-60-16*halfX11)*mm,0*cm),"linePairPhys",linePairLog11,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-60-20*halfX11)*mm,0*cm),"linePairPhys",linePairLog11,watercont_phys,false,0);

    //12lp\cm
    G4double halfX12 = (5./12.)/2.;
    G4Box* linePairBox12 = new G4Box("LinePair",2*cm,halfX12*mm,20*cm);
    G4LogicalVolume* linePairLog12 = new G4LogicalVolume(linePairBox12,Aluminium ,"linePairLog");
    linePairLog12->SetVisAttributes(linepair_att);

    new G4PVPlacement(0,G4ThreeVector(0*cm,-80*mm,0*cm),"linePairPhys",linePairLog12,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-80-4*halfX12)*mm,0*cm),"linePairPhys",linePairLog12,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-80-8*halfX12)*mm,0*cm),"linePairPhys",linePairLog12,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-80-12*halfX12)*mm,0*cm),"linePairPhys",linePairLog12,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-80-16*halfX12)*mm,0*cm),"linePairPhys",linePairLog12,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-80-20*halfX12)*mm,0*cm),"linePairPhys",linePairLog12,watercont_phys,false,0);

    //13lp\cm
    G4double halfX13 = (5./13.)/2.;
    G4Box* linePairBox13 = new G4Box("LinePair",2*cm,halfX13*mm,20*cm);
    G4LogicalVolume* linePairLog13 = new G4LogicalVolume(linePairBox13,Aluminium ,"linePairLog");
    linePairLog13->SetVisAttributes(linepair_att);

    new G4PVPlacement(0,G4ThreeVector(0*cm,-100*mm,0*cm),"linePairPhys",linePairLog13,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-100-4*halfX13)*mm,0*cm),"linePairPhys",linePairLog13,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-100-8*halfX13)*mm,0*cm),"linePairPhys",linePairLog13,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-100-12*halfX13)*mm,0*cm),"linePairPhys",linePairLog13,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-100-16*halfX13)*mm,0*cm),"linePairPhys",linePairLog13,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-100-20*halfX13)*mm,0*cm),"linePairPhys",linePairLog13,watercont_phys,false,0);

    //14lp\cm
    G4double halfX14 = (5./14.)/2.;
    G4Box* linePairBox14 = new G4Box("LinePair",2*cm,halfX14*mm,20*cm);
    G4LogicalVolume* linePairLog14 = new G4LogicalVolume(linePairBox14,Aluminium ,"linePairLog");
    linePairLog14->SetVisAttributes(linepair_att);

    new G4PVPlacement(0,G4ThreeVector(0*cm,-120*mm,0*cm),"linePairPhys",linePairLog14,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-120-4*halfX14)*mm,0*cm),"linePairPhys",linePairLog14,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-120-8*halfX14)*mm,0*cm),"linePairPhys",linePairLog14,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-120-12*halfX14)*mm,0*cm),"linePairPhys",linePairLog14,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-120-16*halfX14)*mm,0*cm),"linePairPhys",linePairLog14,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-120-20*halfX14)*mm,0*cm),"linePairPhys",linePairLog14,watercont_phys,false,0);

    //15lp\cm
    G4double halfX15 = (5./15.)/2.;
    G4Box* linePairBox15 = new G4Box("LinePair",2*cm,halfX15*mm,20*cm);
    G4LogicalVolume* linePairLog15 = new G4LogicalVolume(linePairBox15,Aluminium ,"linePairLog");
    linePairLog15->SetVisAttributes(linepair_att);

    new G4PVPlacement(0,G4ThreeVector(0*cm,-140*mm,0*cm),"linePairPhys",linePairLog15,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-140-4*halfX15)*mm,0*cm),"linePairPhys",linePairLog15,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-140-8*halfX15)*mm,0*cm),"linePairPhys",linePairLog15,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-140-12*halfX15)*mm,0*cm),"linePairPhys",linePairLog15,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-140-16*halfX15)*mm,0*cm),"linePairPhys",linePairLog15,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-140-20*halfX15)*mm,0*cm),"linePairPhys",linePairLog15,watercont_phys,false,0);
    /*
    //16lp\cm
    G4double halfX16 = (5./16.)/2.;
    G4Box* linePairBox16 = new G4Box("LinePair",2*cm,halfX16*mm,4*cm);
    G4LogicalVolume* linePairLog16 = new G4LogicalVolume(linePairBox16,Aluminium ,"linePairLog");
    linePairLog16->SetVisAttributes(linepair_att);

    new G4PVPlacement(0,G4ThreeVector(0*cm,140*mm,7.5*cm),"linePairPhys",linePairLog16,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(140-4*halfX16)*mm,7.5*cm),"linePairPhys",linePairLog16,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(140-8*halfX16)*mm,7.5*cm),"linePairPhys",linePairLog16,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(140-12*halfX16)*mm,7.5*cm),"linePairPhys",linePairLog16,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(140-16*halfX16)*mm,7.5*cm),"linePairPhys",linePairLog16,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(140-20*halfX16)*mm,7.5*cm),"linePairPhys",linePairLog16,watercont_phys,false,0);

    //17lp\cm
    G4double halfX17 = (5./17.)/2.;
    G4Box* linePairBox17 = new G4Box("LinePair",2*cm,halfX17*mm,4*cm);
    G4LogicalVolume* linePairLog17 = new G4LogicalVolume(linePairBox17,Aluminium ,"linePairLog");
    linePairLog17->SetVisAttributes(linepair_att);

    new G4PVPlacement(0,G4ThreeVector(0*cm,120*mm,7.5*cm),"linePairPhys",linePairLog17,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(120-4*halfX17)*mm,7.5*cm),"linePairPhys",linePairLog17,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(120-8*halfX17)*mm,7.5*cm),"linePairPhys",linePairLog17,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(120-12*halfX17)*mm,7.5*cm),"linePairPhys",linePairLog17,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(120-16*halfX17)*mm,7.5*cm),"linePairPhys",linePairLog17,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(120-20*halfX17)*mm,7.5*cm),"linePairPhys",linePairLog17,watercont_phys,false,0);

    //18lp\cm
    G4double halfX18 = (5./18.)/2.;
    G4Box* linePairBox18 = new G4Box("LinePair",2*cm,halfX18*mm,4*cm);
    G4LogicalVolume* linePairLog18 = new G4LogicalVolume(linePairBox18,Aluminium ,"linePairLog");
    linePairLog18->SetVisAttributes(linepair_att);

    new G4PVPlacement(0,G4ThreeVector(0*cm,100*mm,7.5*cm),"linePairPhys",linePairLog18,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(100-4*halfX18)*mm,7.5*cm),"linePairPhys",linePairLog18,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(100-8*halfX18)*mm,7.5*cm),"linePairPhys",linePairLog18,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(100-12*halfX18)*mm,7.5*cm),"linePairPhys",linePairLog18,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(100-16*halfX18)*mm,7.5*cm),"linePairPhys",linePairLog18,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(100-20*halfX18)*mm,7.5*cm),"linePairPhys",linePairLog18,watercont_phys,false,0);

    //19lp\cm
    G4double halfX19 = (5./19.)/2.;
    G4Box* linePairBox19 = new G4Box("LinePair",2*cm,halfX19*mm,4*cm);
    G4LogicalVolume* linePairLog19 = new G4LogicalVolume(linePairBox19,Aluminium ,"linePairLog");
    linePairLog19->SetVisAttributes(linepair_att);

    new G4PVPlacement(0,G4ThreeVector(0*cm,80*mm,7.5*cm),"linePairPhys",linePairLog19,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(80-4*halfX19)*mm,7.5*cm),"linePairPhys",linePairLog19,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(80-8*halfX19)*mm,7.5*cm),"linePairPhys",linePairLog19,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(80-12*halfX19)*mm,7.5*cm),"linePairPhys",linePairLog19,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(80-16*halfX19)*mm,7.5*cm),"linePairPhys",linePairLog19,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(80-20*halfX19)*mm,7.5*cm),"linePairPhys",linePairLog19,watercont_phys,false,0);

    //20lp\cm
    G4double halfX20 = (5./20.)/2.;
    G4Box* linePairBox20 = new G4Box("LinePair",2*cm,halfX20*mm,4*cm);
    G4LogicalVolume* linePairLog20 = new G4LogicalVolume(linePairBox20,Aluminium ,"linePairLog");
    linePairLog20->SetVisAttributes(linepair_att);

    new G4PVPlacement(0,G4ThreeVector(0*cm,60*mm,7.5*cm),"linePairPhys",linePairLog20,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(60-4*halfX20)*mm,7.5*cm),"linePairPhys",linePairLog20,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(60-8*halfX20)*mm,7.5*cm),"linePairPhys",linePairLog20,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(60-12*halfX20)*mm,7.5*cm),"linePairPhys",linePairLog20,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(60-16*halfX20)*mm,7.5*cm),"linePairPhys",linePairLog20,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(60-20*halfX20)*mm,7.5*cm),"linePairPhys",linePairLog20,watercont_phys,false,0);

    //21lp\cm
    G4double halfX21 = (5./21.)/2.;
    G4Box* linePairBox21 = new G4Box("LinePair",2*cm,halfX21*mm,4*cm);
    G4LogicalVolume* linePairLog21 = new G4LogicalVolume(linePairBox21,Aluminium ,"linePairLog");
    linePairLog21->SetVisAttributes(linepair_att);

    new G4PVPlacement(0,G4ThreeVector(0*cm,40*mm,7.5*cm),"linePairPhys",linePairLog21,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(40-4*halfX21)*mm,7.5*cm),"linePairPhys",linePairLog21,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(40-8*halfX21)*mm,7.5*cm),"linePairPhys",linePairLog21,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(40-12*halfX21)*mm,7.5*cm),"linePairPhys",linePairLog21,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(40-16*halfX21)*mm,7.5*cm),"linePairPhys",linePairLog21,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(40-20*halfX21)*mm,7.5*cm),"linePairPhys",linePairLog21,watercont_phys,false,0);

    //22lp\cm
    G4double halfX22 = (5./22.)/2.;
    G4Box* linePairBox22 = new G4Box("LinePair",2*cm,halfX22*mm,4*cm);
    G4LogicalVolume* linePairLog22 = new G4LogicalVolume(linePairBox22,Aluminium ,"linePairLog");
    linePairLog22->SetVisAttributes(linepair_att);

    new G4PVPlacement(0,G4ThreeVector(0*cm,20*mm,7.5*cm),"linePairPhys",linePairLog22,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(20-4*halfX22)*mm,7.5*cm),"linePairPhys",linePairLog22,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(20-8*halfX22)*mm,7.5*cm),"linePairPhys",linePairLog22,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(20-12*halfX22)*mm,7.5*cm),"linePairPhys",linePairLog22,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(20-16*halfX22)*mm,7.5*cm),"linePairPhys",linePairLog22,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(20-20*halfX22)*mm,7.5*cm),"linePairPhys",linePairLog22,watercont_phys,false,0);

    //23lp\cm
    G4double halfX23 = (5./23.)/2.;
    G4Box* linePairBox23 = new G4Box("LinePair",2*cm,halfX23*mm,4*cm);
    G4LogicalVolume* linePairLog23 = new G4LogicalVolume(linePairBox23,Aluminium ,"linePairLog");
    linePairLog23->SetVisAttributes(linepair_att);

    new G4PVPlacement(0,G4ThreeVector(0*cm,0*mm,7.5*cm),"linePairPhys",linePairLog23,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(0-4*halfX23)*mm,7.5*cm),"linePairPhys",linePairLog23,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(0-8*halfX23)*mm,7.5*cm),"linePairPhys",linePairLog23,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(0-12*halfX23)*mm,7.5*cm),"linePairPhys",linePairLog23,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(0-16*halfX23)*mm,7.5*cm),"linePairPhys",linePairLog23,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(0-20*halfX23)*mm,7.5*cm),"linePairPhys",linePairLog23,watercont_phys,false,0);

    //24lp\cm
    G4double halfX24 = (5./24.)/2.;
    G4Box* linePairBox24 = new G4Box("LinePair",2*cm,halfX24*mm,4*cm);
    G4LogicalVolume* linePairLog24 = new G4LogicalVolume(linePairBox24,Aluminium ,"linePairLog");
    linePairLog24->SetVisAttributes(linepair_att);

    new G4PVPlacement(0,G4ThreeVector(0*cm,-20*mm,7.5*cm),"linePairPhys",linePairLog24,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-20-4*halfX24)*mm,7.5*cm),"linePairPhys",linePairLog24,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-20-8*halfX24)*mm,7.5*cm),"linePairPhys",linePairLog24,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-20-12*halfX24)*mm,7.5*cm),"linePairPhys",linePairLog24,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-20-16*halfX24)*mm,7.5*cm),"linePairPhys",linePairLog24,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-20-20*halfX24)*mm,7.5*cm),"linePairPhys",linePairLog24,watercont_phys,false,0);

    //25lp\cm
    G4double halfX25 = (5./25.)/2.;
    G4Box* linePairBox25 = new G4Box("LinePair",2*cm,halfX25*mm,4*cm);
    G4LogicalVolume* linePairLog25 = new G4LogicalVolume(linePairBox25,Aluminium ,"linePairLog");
    linePairLog25->SetVisAttributes(linepair_att);

    new G4PVPlacement(0,G4ThreeVector(0*cm,-40*mm,7.5*cm),"linePairPhys",linePairLog25,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-40-4*halfX25)*mm,7.5*cm),"linePairPhys",linePairLog25,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-40-8*halfX25)*mm,7.5*cm),"linePairPhys",linePairLog25,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-40-12*halfX25)*mm,7.5*cm),"linePairPhys",linePairLog25,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-40-16*halfX25)*mm,7.5*cm),"linePairPhys",linePairLog25,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-40-20*halfX25)*mm,7.5*cm),"linePairPhys",linePairLog25,watercont_phys,false,0);

    //26lp\cm
    G4double halfX26 = (5./26.)/2.;
    G4Box* linePairBox26 = new G4Box("LinePair",2*cm,halfX26*mm,4*cm);
    G4LogicalVolume* linePairLog26 = new G4LogicalVolume(linePairBox26,Aluminium ,"linePairLog");
    linePairLog26->SetVisAttributes(linepair_att);

    new G4PVPlacement(0,G4ThreeVector(0*cm,-60*mm,7.5*cm),"linePairPhys",linePairLog26,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-60-4*halfX26)*mm,7.5*cm),"linePairPhys",linePairLog26,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-60-8*halfX26)*mm,7.5*cm),"linePairPhys",linePairLog26,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-60-12*halfX26)*mm,7.5*cm),"linePairPhys",linePairLog26,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-60-16*halfX26)*mm,7.5*cm),"linePairPhys",linePairLog26,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-60-20*halfX26)*mm,7.5*cm),"linePairPhys",linePairLog26,watercont_phys,false,0);

    //27lp\cm
    G4double halfX27 = (5./27.)/2.;
    G4Box* linePairBox27 = new G4Box("LinePair",2*cm,halfX27*mm,4*cm);
    G4LogicalVolume* linePairLog27 = new G4LogicalVolume(linePairBox27,Aluminium ,"linePairLog");
    linePairLog27->SetVisAttributes(linepair_att);

    new G4PVPlacement(0,G4ThreeVector(0*cm,-80*mm,7.5*cm),"linePairPhys",linePairLog27,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-80-4*halfX27)*mm,7.5*cm),"linePairPhys",linePairLog27,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-80-8*halfX27)*mm,7.5*cm),"linePairPhys",linePairLog27,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-80-12*halfX27)*mm,7.5*cm),"linePairPhys",linePairLog27,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-80-16*halfX27)*mm,7.5*cm),"linePairPhys",linePairLog27,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-80-20*halfX27)*mm,7.5*cm),"linePairPhys",linePairLog27,watercont_phys,false,0);

    //28lp\cm
    G4double halfX28 = (5./28.)/2.;
    G4Box* linePairBox28 = new G4Box("LinePair",2*cm,halfX28*mm,4*cm);
    G4LogicalVolume* linePairLog28 = new G4LogicalVolume(linePairBox28,Aluminium ,"linePairLog");
    linePairLog28->SetVisAttributes(linepair_att);

    new G4PVPlacement(0,G4ThreeVector(0*cm,-100*mm,7.5*cm),"linePairPhys",linePairLog28,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-100-4*halfX28)*mm,7.5*cm),"linePairPhys",linePairLog28,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-100-8*halfX28)*mm,7.5*cm),"linePairPhys",linePairLog28,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-100-12*halfX28)*mm,7.5*cm),"linePairPhys",linePairLog28,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-100-16*halfX28)*mm,7.5*cm),"linePairPhys",linePairLog28,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-100-20*halfX28)*mm,7.5*cm),"linePairPhys",linePairLog28,watercont_phys,false,0);

    //29lp\cm
    G4double halfX29 = (5./29.)/2.;
    G4Box* linePairBox29 = new G4Box("LinePair",2*cm,halfX29*mm,4*cm);
    G4LogicalVolume* linePairLog29 = new G4LogicalVolume(linePairBox29,Aluminium ,"linePairLog");
    linePairLog29->SetVisAttributes(linepair_att);

    new G4PVPlacement(0,G4ThreeVector(0*cm,-120*mm,7.5*cm),"linePairPhys",linePairLog29,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-120-4*halfX29)*mm,7.5*cm),"linePairPhys",linePairLog29,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-120-8*halfX29)*mm,7.5*cm),"linePairPhys",linePairLog29,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-120-12*halfX29)*mm,7.5*cm),"linePairPhys",linePairLog29,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-120-16*halfX29)*mm,7.5*cm),"linePairPhys",linePairLog29,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-120-20*halfX29)*mm,7.5*cm),"linePairPhys",linePairLog29,watercont_phys,false,0);

    //30lp\cm
    G4double halfX30 = (5./30.)/2.;
    G4Box* linePairBox30 = new G4Box("LinePair",2*cm,halfX30*mm,4*cm);
    G4LogicalVolume* linePairLog30 = new G4LogicalVolume(linePairBox30,Aluminium ,"linePairLog");
    linePairLog30->SetVisAttributes(linepair_att);

    new G4PVPlacement(0,G4ThreeVector(0*cm,-140*mm,7.5*cm),"linePairPhys",linePairLog30,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-140-4*halfX30)*mm,7.5*cm),"linePairPhys",linePairLog30,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-140-8*halfX30)*mm,7.5*cm),"linePairPhys",linePairLog30,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-140-12*halfX30)*mm,7.5*cm),"linePairPhys",linePairLog30,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-140-16*halfX30)*mm,7.5*cm),"linePairPhys",linePairLog30,watercont_phys,false,0);
    new G4PVPlacement(0,G4ThreeVector(0*cm,(-140-20*halfX30)*mm,7.5*cm),"linePairPhys",linePairLog30,watercont_phys,false,0);
    */
  }

  //----------------------------------------------------------------------------------------------------------------
  // Las Vegas Constrast Phantom
  //----------------------------------------------------------------------------------------------------------------
  else if (thePhantom =="LasVegas"){

    G4VisAttributes* linepair_att  = new G4VisAttributes(G4Colour(0,0,1));
    linepair_att->SetVisibility(true);
    //linepair_att->SetForceSolid(true);
    G4Material *Aluminium = new G4Material("Aluminum", 13, 26.98*g/mole, 2.7*g/cm3);
    G4Box* LasVegasBox = new G4Box("LasVegasBox",(2.18/2)*cm,10*cm,10*cm);
    G4LogicalVolume* LasVegasLog = new G4LogicalVolume(LasVegasBox,Aluminium ,"LasVegasLog");
    LasVegasLog->SetVisAttributes(linepair_att);
    G4VPhysicalVolume* LasVegasPhys = new G4PVPlacement(rotExt,G4ThreeVector(),"LasVegasPhys",LasVegasLog,cont_phys,false,0);

    double PosX[6]      ={-75,-45,-15,15,45,75}; //mm
    double Diameter[6]  ={0.5/2, 2./2, 4./2, 7./2,   10./2, 15./2}; //mm
    double PosY[5]      ={-60,-30,0,30,60}; //mm
    double Depth[5]     ={0.5/2, 1./2, 2./2, 3.25/2, 4.5/2}; // mm

    G4RotationMatrix* rot1 = new G4RotationMatrix();
    rot1->rotateY(pi/2.);
    for(int ix = 0; ix<6; ix++){
      for(int iy = 0; iy<5; iy++){
	G4ThreeVector Pos = {-1.09*cm + Depth[iy], PosY[iy], PosX[ix]};
	G4Tubs* insert = new G4Tubs("insert",0, Diameter[ix] , Depth[iy] , 0. , twopi);
	G4LogicalVolume* insert_Log = new G4LogicalVolume(insert,air ,Form("insert_log_%d_%d",ix,iy));
	new G4PVPlacement(rot1, Pos,Form("insert_phys_%d_%d",ix,iy),insert_Log, LasVegasPhys,false,0);
      }
    }
  }

  //----------------------------------------------------------------------------------------------------------------
  // Slanted edge phantom
  //----------------------------------------------------------------------------------------------------------------

  else if(thePhantom =="SlantedEdge"){

    G4Box* WaterBox  = new G4Box("WaterBox", 10./2*cm, PhantomHalfY, PhantomHalfZ);
    G4LogicalVolume* WaterBox_log = new G4LogicalVolume(WaterBox, water, "WaterBox_log");
    G4PVPlacement* watercont_phys = new G4PVPlacement(0,G4ThreeVector(0,0,0),"Watercont_phys",WaterBox_log,cont_phys,false,0);


    G4VisAttributes* linepair_att  = new G4VisAttributes(G4Colour(0,0,1));
    linepair_att->SetVisibility(true);
    linepair_att->SetForceSolid(true);
    G4Material *Aluminium = new G4Material("Aluminum", 13, 26.98*g/mole, 2.7*g/cm3);
    //G4Material *InnerBone = theMaterial->ConstructMaterial("InnerBone",1.089);
    //InnerBone->GetIonisation()->SetMeanExcitationEnergy(67.138689*eV);
    G4Box* SlantedBox = new G4Box("SlantedBox",4./2*cm,10./2*cm,10./2*cm);
    G4LogicalVolume* SlantedLog = new G4LogicalVolume(SlantedBox,Aluminium ,"SlantedLog");
    SlantedLog->SetVisAttributes(linepair_att);
    rotExt->rotateX(2.5*pi/180.);
    new G4PVPlacement(rotExt,G4ThreeVector(),"SlantedPhys",SlantedLog,watercont_phys,false,0);
}

  //----------------------------------------------------------------------------------------------------------------
  // Edge phantom
  //----------------------------------------------------------------------------------------------------------------

  else if(thePhantom =="Edge"){
    G4VisAttributes* linepair_att  = new G4VisAttributes(G4Colour(0,0,1));
    linepair_att->SetVisibility(true);
    linepair_att->SetForceSolid(true);
    G4Material *Aluminium = new G4Material("Aluminum", 13, 26.98*g/mole, 2.7*g/cm3);
    //G4Material *InnerBone = theMaterial->ConstructMaterial("InnerBone",1.089);
    //InnerBone->GetIonisation()->SetMeanExcitationEnergy(67.138689*eV);
    G4Box* Box = new G4Box("Box",4./2*cm,10./2*cm,10./2*cm);
    G4LogicalVolume* Log = new G4LogicalVolume(Box,Aluminium ,"Log");
    Log->SetVisAttributes(linepair_att);

    new G4PVPlacement(rotExt,G4ThreeVector(),"Phys",Log,cont_phys,false,0);
}

  //----------------------------------------------------------------------------------------------------------------
  // Many Edge phantom
  //----------------------------------------------------------------------------------------------------------------

  else if(thePhantom =="ManyEdge"){
    G4Material *Aluminium = new G4Material("Aluminum", 13, 26.98*g/mole, 2.7*g/cm3);
    G4Box* Box = new G4Box("Box",4./2*cm,4./2*cm,4./2*cm);
    G4LogicalVolume* Log = new G4LogicalVolume(Box,Aluminium ,"Log");

    new G4PVPlacement(rotExt,G4ThreeVector(-7*cm,-10*cm,-10*cm),"Phys_0",Log,cont_phys,true,0);
    new G4PVPlacement(rotExt,G4ThreeVector(-3.5*cm,-5*cm,-5*cm),"Phys_1",Log,cont_phys,true,1);
    new G4PVPlacement(rotExt,G4ThreeVector(0*cm,0*cm,0*cm),"Phys_2",Log,cont_phys,true,2);
    new G4PVPlacement(rotExt,G4ThreeVector(3.5*cm,5*cm,5*cm),"Phys_3",Log,cont_phys,true,3);
    new G4PVPlacement(rotExt,G4ThreeVector(7*cm,10*cm,10*cm),"Phys_4",Log,cont_phys,true,4);


}

  //----------------------------------------------------------------------------------------------------------------
  // Catphan phantom
  //----------------------------------------------------------------------------------------------------------------

  else if(thePhantom =="Catphan"){

    //Specific list of material for the catphan phantom
    G4double z, a, density;
    G4String name, symbol;
    G4int ncomponents, natoms;

    G4Element* elH  = new G4Element(name="Hydrogen"  ,symbol="H" , z= 1., a=1.01*g/mole);
    G4Element* elC  = new G4Element(name="Carbon"    ,symbol="C" , z= 6., a=12.01*g/mole);
    G4Element* elO  = new G4Element(name="Oxygen"    ,symbol="O" , z= 8., a=16.00*g/mole);
    G4Element* elF  = new G4Element(name="Fluorine"  ,symbol="F" , z= 9., a=19.00*g/mole);
    G4Element* elAl = new G4Element(name="Aluminium" ,symbol="Al", z= 13.,a=26.98*g/mole);
    G4Element* elW  = new G4Element(name="Tungsten"  ,symbol="W",  z= 74.,a=183.84*g/mole);

    G4Material *Teflon = new G4Material("Teflon",density=2.16*g/cm3,ncomponents=2); //C2F4
    Teflon->AddElement(elC,natoms=2);
    Teflon->AddElement(elF,natoms=4);

    G4Material *PMP    = new G4Material("PMP",density=0.83*g/cm3,ncomponents=2); //H14C7
    PMP->AddElement(elH,natoms=14);
    PMP->AddElement(elC,natoms=7);

    G4Material *LDPE    = new G4Material("LDPE",density=0.92*g/cm3,ncomponents=2); //H4C2
    LDPE->AddElement(elH,natoms=4);
    LDPE->AddElement(elC,natoms=2);

    G4Material *Polystyrene    = new G4Material("Polystyrene",density=1.05*g/cm3,ncomponents=2); //H8C8
    Polystyrene->AddElement(elH,natoms=8);
    Polystyrene->AddElement(elC,natoms=8);

    G4Material *Delrin    = new G4Material("Delrin",density=1.41*g/cm3,ncomponents=3); //H2C1O1
    Delrin->AddElement(elH,natoms=2);
    Delrin->AddElement(elC,natoms=1);
    Delrin->AddElement(elO,natoms=1);

    G4Material *Acrylic    = new G4Material("Acrylic",density=1.18*g/cm3,ncomponents=3); //H8C5O2
    Acrylic->AddElement(elH,natoms=8);
    Acrylic->AddElement(elC,natoms=5);
    Acrylic->AddElement(elO,natoms=2);

    G4Material* Aluminium = new G4Material("Aluminium",density=2.70*g/cm3,ncomponents=1);
    Aluminium->AddElement(elAl,natoms=1);

    G4Material* TungstenCarbide = new G4Material("TungstenCarbide", density=15.63*g/cm3,ncomponents=2);
    TungstenCarbide->AddElement(elC,natoms=1);
    TungstenCarbide->AddElement(elW,natoms=1);

    G4VisAttributes* att3  = new G4VisAttributes(G4Colour(0,1,1));
    att3->SetVisibility(true);

    G4Tubs* Casing= new G4Tubs("Casing",0,200*mm/2,160*mm/2,0,2*pi);
    G4LogicalVolume* Casing_Log = new G4LogicalVolume(Casing,water,"Casing_Log");
    Casing_Log->SetVisAttributes(att3);

    G4VPhysicalVolume* CasingPhys = new G4PVPlacement(rotExt,G4ThreeVector(),"CasingWorld",Casing_Log,cont_phys,false,0);
    // CTP404 Slice Geometry and sensitometry module
    G4VisAttributes* att  = new G4VisAttributes(G4Colour(0,0,1));
    att->SetVisibility(true);

    G4VisAttributes* sph_att  = new G4VisAttributes(G4Colour(0,1,0));
    sph_att->SetVisibility(true);
    sph_att->SetForceSolid(true);

    G4Tubs* CTP404 = new G4Tubs("CTP404",0,150*mm/2,20/2*mm,0,2*pi);
    G4LogicalVolume* CTP404Log = new G4LogicalVolume(CTP404,water ,"CTP404Log"); // Water

    G4VPhysicalVolume* CTP404World = new G4PVPlacement(0,G4ThreeVector(0,0,70*mm),"CTP404World",CTP404Log,CasingPhys,false,0);

    // 4 Air and Teflons rods
    G4Tubs* small_rod  = new G4Tubs("small_rod", 0, 3*mm/2, 20/2*mm, 0,2*pi);
    G4LogicalVolume* small_rod_1 = new G4LogicalVolume(small_rod,air ,"small_rod_1"); // Air
    G4LogicalVolume* small_rod_2 = new G4LogicalVolume(small_rod,air ,"small_rod_2"); // Air
    G4LogicalVolume* small_rod_3 = new G4LogicalVolume(small_rod,air ,"small_rod_3"); // Air
    G4LogicalVolume* small_rod_4 = new G4LogicalVolume(small_rod,Teflon ,"small_rod_4"); // Teflon

    new G4PVPlacement(0,G4ThreeVector(-2.5*cm,-2.5*cm,0),"small_rod_1World",small_rod_1,CTP404World,false,0);
    new G4PVPlacement(0,G4ThreeVector(-2.5*cm, 2.5*cm,0),"small_rod_2World",small_rod_2,CTP404World,false,0);
    new G4PVPlacement(0,G4ThreeVector( 2.5*cm,-2.5*cm,0),"small_rod_3World",small_rod_3,CTP404World,false,0);
    new G4PVPlacement(0,G4ThreeVector( 2.5*cm, 2.5*cm,0),"small_rod_4World",small_rod_4,CTP404World,false,0);

    /*
    // 5 Acrylic Sphere
    G4Sphere *acrilSphere1 = new G4Sphere ("acrilSphere1", 0, 10/2*mm, 0,2*pi,0,2*pi);
    G4Sphere *acrilSphere2 = new G4Sphere ("acrilSphere2", 0, 8/2*mm, 0,2*pi,0,2*pi);
    G4Sphere *acrilSphere3 = new G4Sphere ("acrilSphere3", 0, 6/2*mm, 0,2*pi,0,2*pi);
    G4Sphere *acrilSphere4 = new G4Sphere ("acrilSphere4", 0, 4/2*mm, 0,2*pi,0,2*pi);
    G4Sphere *acrilSphere5 = new G4Sphere ("acrilSphere5", 0, 2/2*mm, 0,2*pi,0,2*pi);

    G4LogicalVolume *acrilLog1 = new G4LogicalVolume(acrilSphere1,Acrylic,"acrilLog1"); //Acrylic
    G4LogicalVolume *acrilLog2 = new G4LogicalVolume(acrilSphere2,Acrylic,"acrilLog2"); //Acrylic
    G4LogicalVolume *acrilLog3 = new G4LogicalVolume(acrilSphere3,Acrylic,"acrilLog3"); //Acrylic
    G4LogicalVolume *acrilLog4 = new G4LogicalVolume(acrilSphere4,Acrylic,"acrilLog4"); //Acrylic
    G4LogicalVolume *acrilLog5 = new G4LogicalVolume(acrilSphere5,Acrylic,"acrilLog5"); //Acrylic

    acrilLog1->SetVisAttributes(sph_att);
    acrilLog2->SetVisAttributes(sph_att);
    acrilLog3->SetVisAttributes(sph_att);
    acrilLog4->SetVisAttributes(sph_att);
    acrilLog5->SetVisAttributes(sph_att);

    G4ThreeVector placement(0,15*mm,0);
    placement.rotateZ(2*pi/5);
    new G4PVPlacement(0,placement,"acril_1_World",acrilLog1,CTP404World,false,0);
    placement.rotateZ(2*pi/5);
    new G4PVPlacement(0,placement,"acril_2_World",acrilLog2,CTP404World,false,0);
    placement.rotateZ(2*pi/5);
    new G4PVPlacement(0,placement,"acril_3_World",acrilLog3,CTP404World,false,0);
    placement.rotateZ(2*pi/5);
    new G4PVPlacement(0,placement,"acril_4_World",acrilLog4,CTP404World,false,0);
    placement.rotateZ(2*pi/5);
    new G4PVPlacement(0,placement,"acril_5_World",acrilLog5,CTP404World,false,0);

    */
    // 8 Sensitometry 12 mm Large Rod
    G4Tubs* large_rod  = new G4Tubs("large_rod", 0, 12*mm/2, 20/2*mm, 0,2*pi);
    G4LogicalVolume* large_rod_1 = new G4LogicalVolume(large_rod,Acrylic ,"large_rod_1");  //Acrylic
    G4LogicalVolume* large_rod_2 = new G4LogicalVolume(large_rod,Polystyrene ,"large_rod_2");  //Polystyrene
    G4LogicalVolume* large_rod_3 = new G4LogicalVolume(large_rod,LDPE ,"large_rod_3");  //LDPE
    G4LogicalVolume* large_rod_4 = new G4LogicalVolume(large_rod,Delrin ,"large_rod_4"); //Delrin
    G4LogicalVolume* large_rod_5 = new G4LogicalVolume(large_rod,Teflon ,"large_rod_5"); //Teflon
    G4LogicalVolume* large_rod_6 = new G4LogicalVolume(large_rod,air ,"large_rod_6"); //Air
    G4LogicalVolume* large_rod_7 = new G4LogicalVolume(large_rod,PMP ,"large_rod_7"); //PMP
    G4LogicalVolume* large_rod_8 = new G4LogicalVolume(large_rod,water ,"large_rod_8"); //Water

    G4ThreeVector placement2(0,6*cm,0);
    new G4PVPlacement(0,placement2,"large_rod_1World",large_rod_1,CTP404World,false,0);
    placement2.rotateZ(pi/6);
    new G4PVPlacement(0,placement2,"large_rod_4World",large_rod_2,CTP404World,false,0);
    placement2.rotateZ(2*pi/6);
    new G4PVPlacement(0,placement2,"large_rod_4World",large_rod_3,CTP404World,false,0);
    placement2.rotateZ(2*pi/6);
    new G4PVPlacement(0,placement2,"large_rod_4World",large_rod_4,CTP404World,false,0);
    placement2.rotateZ(pi/6);
    new G4PVPlacement(0,placement2,"large_rod_4World",large_rod_5,CTP404World,false,0);
    placement2.rotateZ(pi/6);
    new G4PVPlacement(0,placement2,"large_rod_4World",large_rod_6,CTP404World,false,0);
    placement2.rotateZ(2*pi/6);
    new G4PVPlacement(0,placement2,"large_rod_4World",large_rod_7,CTP404World,false,0);
    placement2.rotateZ(2*pi/6);
    new G4PVPlacement(0,placement2,"large_rod_4World",large_rod_8,CTP404World,false,0);


    G4Tubs* wire_ramp  = new G4Tubs("wire_ramp", 0, 1*mm, 5.12/2*cm, 0,2*pi); //5.12 come from the 23 degree wire ramp and the thickness Z of the catphan of 2.0cm (Hyp = 2.0/sin(23))
    G4LogicalVolume* wire_log = new G4LogicalVolume(wire_ramp,Aluminium ,"wire_log"); // Aluminum

    G4RotationMatrix* rotWire = new G4RotationMatrix();
    rotWire->rotateY((90-23)*pi/180.);
    new G4PVPlacement(rotWire,G4ThreeVector(0,4*cm,0),"wire_1_World",wire_log,CTP404World,false,0);
    new G4PVPlacement(rotWire,G4ThreeVector(0,-4*cm,0),"wire_1_World",wire_log,CTP404World,false,0);

    G4RotationMatrix* rotWire2 = new G4RotationMatrix();
    rotWire2->rotateX((90-23)*pi/180.);
    new G4PVPlacement(rotWire2,G4ThreeVector(4*cm,0,0),"wire_1_World",wire_log,CTP404World,false,0);
    new G4PVPlacement(rotWire2,G4ThreeVector(-4*cm,0,0),"wire_1_World",wire_log,CTP404World,false,0);

    //CTP528 High Resolution Module
    G4VisAttributes* att2  = new G4VisAttributes(G4Colour(1,0,1));
    att2->SetVisibility(true);

    G4Tubs* CTP528 = new G4Tubs("CTP528",0,150*mm/2,40/2*mm,0,2*pi);
    G4LogicalVolume* CTP528Log = new G4LogicalVolume(CTP528,water ,"CTP528Log");
    CTP528Log->SetVisAttributes(att2);
    G4VPhysicalVolume* CTP528World = new G4PVPlacement(0,G4ThreeVector(0,0,40*mm),"CTP528World",CTP528Log,CasingPhys,false,0);

    G4Sphere *bead = new G4Sphere ("acrilSphere1", 0, 0.28/2*mm, 0,2*pi,0,2*pi);
    G4LogicalVolume* beadLog = new G4LogicalVolume(bead,TungstenCarbide,"bead_log");  // Tungsten Carbide
    beadLog->SetVisAttributes(sph_att);
    new G4PVPlacement(0,G4ThreeVector(0,20*mm,0),"bead_1_world",beadLog,CTP528World,false,0);
    new G4PVPlacement(0,G4ThreeVector(0,-20*mm,0),"bead_1_world",beadLog,CTP528World,false,0);
    std::vector<float> ctp528_anglesCoefficients = { -3.8654e-6f, 1.9470e-4f, -2.3863e-3f, -3.6253e-2f, 1.3968f, -3.0373e1f, 1.9902e2f};// Lp/Cm to angle [deg]
    for(int i =1;i<21;i++){

      double GapSize = 0.5/float(i);
      G4ThreeVector placement4(0,47.5*mm,0);
      int nbLP = 5;
      G4RotationMatrix* rotLP = new G4RotationMatrix();
      float angle = 0.f;
      float x=1.;
      const int Ns = static_cast<int>(ctp528_anglesCoefficients.size());
      for(int j=0; j<Ns;++j) {
        angle += x*ctp528_anglesCoefficients[Ns-1-j];
        x *=i;
      }
      placement4.rotateZ((270+angle)*degree);
      rotLP->rotateZ( (270-angle)*degree);

      if(i==1) nbLP =2;
      else if(i==2) nbLP =3;
      else if(i<6)  nbLP = 4;
      G4Box* container_box = new G4Box("container_box", (nbLP-0.5)*GapSize*cm,4*mm,2.5*mm);
      G4LogicalVolume* container_log = new G4LogicalVolume(container_box,water,"container_log"); //Water
      G4VPhysicalVolume* container_world = new G4PVPlacement(rotLP,placement4,"container_World",container_log,CTP528World,false,0); // Container for the LP
      for(int j=0;j<nbLP;j++){

	G4Box* LinePair = new G4Box("LinePair",GapSize/2*cm,4*mm,4*mm);
	G4LogicalVolume* LinePair_Log = new G4LogicalVolume(LinePair,Aluminium,"LinePair_Log"); // Aluminium
	LinePair_Log->SetVisAttributes(sph_att);
	new G4PVPlacement(0,G4ThreeVector(-(nbLP-0.5)*GapSize*cm + GapSize/2*cm +  2*j*GapSize*cm,0,0),"LinePair_world",LinePair_Log,container_world,false,0);
      }
    }
    // CTP515 low contrast module
    G4Tubs* CTP515 = new G4Tubs("CTP515",0,150*mm/2,40/2*mm,0,2*pi);
    G4LogicalVolume* CTP515Log = new G4LogicalVolume(CTP515,water ,"CTP515Log"); // Water
    G4VPhysicalVolume* CTP515World = new G4PVPlacement(0,G4ThreeVector(0,0,0*mm),"CTP515World",CTP515Log,CasingPhys,false,0);

    // Supra-slices
    //--------------
    float ss_radius = 50*mm;
    vector<double> ss_insdia = { { 15*mm, 9*mm, 8*mm, 7*mm, 6*mm, 5*mm, 4*mm, 3*mm, 2*mm  } };
    vector<double> ss_insang = { { 0, 20, 35.5f, 48.5f, 62, 74, 82, 90, 95  } };

    // Supra-slice 1%
    float angleStart = -90;
    for(int i=0; i<9; ++i) {
      float angle =  angleStart - ss_insang[i];
      G4ThreeVector placement(ss_radius,0,0);
      placement.rotateZ(angle*pi/180);
      G4Material* temp_water = theMaterial->ConstructMaterial("Water",1.01);//water + 1.0%
      G4Tubs* insert   = new G4Tubs(Form("CT515_SupraSlice_1.0_%i_vol",i),0,ss_insdia[i]*mm/2,40/2*mm,0,2*pi);
      G4LogicalVolume* insertLog = new G4LogicalVolume(insert,temp_water ,Form("CT515_SupraSlice_1.0_%i_log",i));
      new G4PVPlacement(0, placement,Form("CT515_SupraSlice_1.0_%i_phys",i) ,insertLog,CTP515World,false,0);
    }

    // Supra-slice 0.3%
    angleStart = 150.f;
    for(int i=0; i<9; ++i) {
      float angle =  angleStart - ss_insang[i];
      G4ThreeVector placement(ss_radius,0,0);
      placement.rotateZ(angle*pi/180);
      G4Material* temp_water = theMaterial->ConstructMaterial("Water",1.003);// water + 0.3%
      G4Tubs* insert   = new G4Tubs(Form("CT515_SupraSlice_0.3_%i_vol",i),0,ss_insdia[i]*mm/2,40/2*mm,0,2*pi);
      G4LogicalVolume* insertLog = new G4LogicalVolume(insert,temp_water ,Form("CT515_SupraSlice_0.3_%i_log",i));
      new G4PVPlacement(0, placement,Form("CT515_SupraSlice_0.3_%i_phys",i) ,insertLog,CTP515World,false,0);
    }

    // Supra-slice 0.5%
    angleStart = 30.f;
    for(int i=0; i<9; ++i) {
      float angle =  angleStart - ss_insang[i];
      G4ThreeVector placement(ss_radius,0,0);
      placement.rotateZ(angle*pi/180);
      G4Material* temp_water = theMaterial->ConstructMaterial("Water",1.005);// water + 0.5%
      G4Tubs* insert   = new G4Tubs(Form("CT515_SupraSlice_0.5_%i_vol",i),0,ss_insdia[i]*mm/2,40/2*mm,0,2*pi);
      G4LogicalVolume* insertLog = new G4LogicalVolume(insert,temp_water ,Form("CT515_SupraSlice_0.5_%i_log",i));
      new G4PVPlacement(0, placement,Form("CT515_SupraSlice_0.5_%i_phys",i) ,insertLog,CTP515World,false,0);
    }

    // Subslices
    //--------------
    float su_radius = 25.*mm;
    vector<double> su_insdia = { { 9.*mm, 7.*mm, 5.*mm, 3.*mm  } };
    vector<double> su_insang = { { 21, 50, 74, 90  } };

    // Subslices 7mm
    angleStart = -90;
    for(int i=0; i<4; ++i) {
      float angle =  angleStart - su_insang[i];
      G4ThreeVector placement(su_radius,0,0);
      placement.rotateZ(angle*pi/180);
      G4Material* temp_water = theMaterial->ConstructMaterial("Water",1.01);//water + 1.0%
      G4Tubs* insert   = new G4Tubs(Form("CT515_SubSlice_7mm_%i_vol",i),0,su_insdia[i]*mm/2, 7./2*mm,0,2*pi);
      G4LogicalVolume* insertLog = new G4LogicalVolume(insert,temp_water ,Form("CT515_SubSlice_7mm_%i_log",i));
      new G4PVPlacement(0, placement,Form("CT515_SubSlice_7mm_%i_phys",i) ,insertLog,CTP515World,false,0);
    }

    // Subslices 3mm
    angleStart = 150;
    for(int i=0; i<4; ++i) {
      float angle =  angleStart - su_insang[i];
      G4ThreeVector placement(su_radius,0,0);
      placement.rotateZ(angle*pi/180);
      G4Material* temp_water = theMaterial->ConstructMaterial("Water",1.01);//water + 1.0%
      G4Tubs* insert   = new G4Tubs(Form("CT515_SubSlice_3mm_%i_vol",i),0,su_insdia[i]*mm/2,  3/2*mm,0,2*pi);
      G4LogicalVolume* insertLog = new G4LogicalVolume(insert,temp_water ,Form("CT515_SubSlice_3mm_%i_log",i));
      new G4PVPlacement(0, placement,Form("CT515_SubSlice_3mm_%i_phys",i) ,insertLog,CTP515World,false,0);
    }

    // Subslices 5mm
    angleStart = 30;
    for(int i=0; i<4; ++i) {
      float angle =  angleStart - su_insang[i];
      G4ThreeVector placement(su_radius,0,0);
      placement.rotateZ(angle*pi/180);
      G4Material* temp_water = theMaterial->ConstructMaterial("Water",1.01);//water + 1.0%
      G4Tubs* insert   = new G4Tubs(Form("CT515_SubSlice_5mm_%i_vol",i),0,su_insdia[i]*mm/2,  5/2*mm,0,2*pi);
      G4LogicalVolume* insertLog = new G4LogicalVolume(insert,temp_water ,Form("CT515_SubSlice_5mm_%i_log",i));
      new G4PVPlacement(0, placement,Form("CT515_SubSlice_5mm_%i_phys",i) ,insertLog,CTP515World,false,0);
    }

    // CTP486 Uniformity module
    G4Tubs* CTP486 = new G4Tubs("CTP486",0,150*mm/2,60/2*mm,0,2*pi);
    G4LogicalVolume* CTP486Log = new G4LogicalVolume(CTP486,water ,"CTP486Log"); // Water
    new G4PVPlacement(0,G4ThreeVector(0,0,-50*mm),"CTP486World",CTP486Log,CasingPhys,false,0);



  }
  //----------------------------------------------------------------------------------------------------------------
  // Lung phantom (Charles-Antoine Collins-Fekete)
  //----------------------------------------------------------------------------------------------------------------

  else if(thePhantom =="Lung"){
    G4Material *AverageMaleSoftTissue          = theMaterial->ConstructMaterial("Average_male_soft_tissue_Woodard_1986",1.03);
    G4Material *LungInflated                   = theMaterial->ConstructMaterial("Lung_Woodard_1986",0.26);
    G4Material *RedMarrow                      = theMaterial->ConstructMaterial("red_marrow_Woodard_1986",1.03);
    G4Material *Rib2ndICRUreport46             = theMaterial->ConstructMaterial("Rib_2nd_ICRU_report_46",1.41);
    G4Material *brainwhiteWoodard1986          = theMaterial->ConstructMaterial("brain_white_Woodard_1986",1.04);
    G4Material *vertebralcolumnC4ICRUreport46  = theMaterial->ConstructMaterial("vertebral_column_C4_ICRU_report_46",1.42);
    G4VisAttributes* att     = new G4VisAttributes(true,G4Colour(0,1,1));
    G4VisAttributes* att2    = new G4VisAttributes(true,G4Colour(1,0,1));
    G4VisAttributes* att3    = new G4VisAttributes(true,G4Colour(1,1,1));
    G4VisAttributes* att4    = new G4VisAttributes(true,G4Colour(0,0,1));

    //rotExt->rotateY(pi/2.);
    rotExt->rotateZ(pi/2.);

    //Torso
    G4EllipticalTube* Torso      = new G4EllipticalTube("Torso",135*mm,75*mm,150*mm);
    G4EllipticalTube* Arm        = new G4EllipticalTube("LeftArm",37.5*mm,37.5*mm,150*mm);

    //Left Arm
    G4UnionSolid* TorsoLeftArm   = new G4UnionSolid("Torso+Left", Torso, Arm, 0,G4ThreeVector(-150*mm,0,0*mm));

    //Right Arm
    G4UnionSolid* TorsoBothArm   = new G4UnionSolid("Torso+Both", TorsoLeftArm,Arm, 0, G4ThreeVector(150*mm,0,0*mm));
    G4LogicalVolume* Torso_Log   = new G4LogicalVolume(TorsoBothArm,AverageMaleSoftTissue,"Torso_Log");
    Torso_Log->SetVisAttributes(att);
    G4VPhysicalVolume* TorsoPhys = new G4PVPlacement(rotExt,G4ThreeVector(0,0,0),"TorsoWorld",Torso_Log,cont_phys,false,0);

    //Oesophagus
    G4EllipticalTube* Oesophagus = new G4EllipticalTube("Oesophagus",5*mm,5*mm,150*mm);
    G4LogicalVolume* Oesophagus_Log = new G4LogicalVolume(Oesophagus,air,"Oesophagus_Log");
    Oesophagus_Log->SetVisAttributes(att3);
    new G4PVPlacement(0,G4ThreeVector(10*mm, -12.5*mm,0*mm),"OesophagusWorld",Oesophagus_Log,TorsoPhys,false,0);

    //SpineBlock
    G4EllipticalTube* SpineBlock= new G4EllipticalTube("SpineBlock",7.5*mm,7.5*mm,150*mm);
    G4LogicalVolume* SpineBlock_Log = new G4LogicalVolume(SpineBlock, vertebralcolumnC4ICRUreport46,"SpineBlock_Log");
    SpineBlock_Log->SetVisAttributes(att3);
    new G4PVPlacement(0,G4ThreeVector(0, 66.0*mm,0),"SpineBlockWorld",SpineBlock_Log,TorsoPhys,false,0);

    //SpineVertebra
    G4EllipticalTube* SpineVertebra= new G4EllipticalTube("SpineVertebra",12*mm,12*mm,150*mm);
    G4LogicalVolume* SpineVertebra_Log = new G4LogicalVolume(SpineVertebra, vertebralcolumnC4ICRUreport46,"SpineVertebra_Log");
    SpineVertebra_Log->SetVisAttributes(att2);
    G4VPhysicalVolume* VertPhys = new G4PVPlacement(0,G4ThreeVector(0, 50*mm,0),"SpineVertebraWorld",SpineVertebra_Log,TorsoPhys,false,0);

    //SpineNerve
    G4EllipticalTube* SpineNerve= new G4EllipticalTube("SpineNerve",7.5*mm,7.5*mm,150*mm);
    G4LogicalVolume* SpineNerve_Log = new G4LogicalVolume(SpineNerve, brainwhiteWoodard1986,"SpineNerve_Log");
    SpineNerve_Log->SetVisAttributes(att3);
    new G4PVPlacement(0,G4ThreeVector(0, 0*mm,0),"SpineNerveWorld",SpineNerve_Log,VertPhys,false,0);

    //SpineLeft
    G4RotationMatrix* rotLeft = new G4RotationMatrix(pi/4,0,0);
    G4EllipticalTube* SpineLeft= new G4EllipticalTube("SpineLeft",3*mm,6*mm,150*mm);
    G4LogicalVolume* SpineLeft_Log = new G4LogicalVolume(SpineLeft, vertebralcolumnC4ICRUreport46,"SpineLeft_Log");
    SpineLeft_Log->SetVisAttributes(att3);
    new G4PVPlacement(rotLeft,G4ThreeVector(-12*mm, 58*mm,0),"SpineLeftWorld",SpineLeft_Log,TorsoPhys,false,0);

    //SpineRight
    G4RotationMatrix* rotRight = new G4RotationMatrix(-pi/4,0,0);
    G4EllipticalTube* SpineRight= new G4EllipticalTube("SpineRight",3*mm,6*mm,150*mm);
    G4LogicalVolume* SpineRight_Log = new G4LogicalVolume(SpineRight, vertebralcolumnC4ICRUreport46,"SpineRight_Log");
    SpineRight_Log->SetVisAttributes(att3);
    new G4PVPlacement(rotRight,G4ThreeVector(12*mm, 58*mm,0),"SpineRightWorld",SpineRight_Log,TorsoPhys,false,0);

    //Central Tissue
    G4EllipticalTube* CentralTissue= new G4EllipticalTube("CentralTissue",37.5*mm,48.0*mm,150*mm);
    G4LogicalVolume* CentralTissue_Log = new G4LogicalVolume(CentralTissue,AverageMaleSoftTissue,"CentralTissue_Log");
    CentralTissue_Log->SetVisAttributes(att2);
    new G4PVPlacement(0,G4ThreeVector(0,-2.5*mm,0),"CentralTissueWorld",CentralTissue_Log,TorsoPhys,false,0);

    //Ribs Volume
    G4EllipticalTube* InsideRibs     = new G4EllipticalTube("InsideRibs",10.5*mm,1.5*mm,150*mm);
    G4LogicalVolume*  InsideRibs_Log = new G4LogicalVolume(InsideRibs,RedMarrow,"InsideRibs_Log");
    G4EllipticalTube* OutsideRibs    = new G4EllipticalTube("InsideRibs",13.5*mm,4.5*mm,150*mm);
    G4LogicalVolume*  OutsideRibs_Log = new G4LogicalVolume(OutsideRibs,Rib2ndICRUreport46,"OutsideRibs_Log");
    InsideRibs_Log->SetVisAttributes(att4);
    OutsideRibs_Log->SetVisAttributes(att3);

    //Placement Ribs Left Lung
    //1
    new G4PVPlacement(0,G4ThreeVector(-60*mm,-49*mm,0*mm),"InsideRibsWorld_Left",InsideRibs_Log,TorsoPhys,false,0);
    new G4PVPlacement(0,G4ThreeVector(-60*mm,-49*mm,0*mm),"OutsideRibsWorld_Left",OutsideRibs_Log,TorsoPhys,false,0);
    //2
    G4RotationMatrix* rot2_l = new G4RotationMatrix(-75*pi/180,0,0);
    new G4PVPlacement(rot2_l,G4ThreeVector(-95*mm,-14*mm,0*mm),"InsideRibsWorld_Left",InsideRibs_Log,TorsoPhys,false,0);
    new G4PVPlacement(rot2_l,G4ThreeVector(-95*mm,-14*mm,0*mm),"OutsideRibsWorld_Left",OutsideRibs_Log,TorsoPhys,false,0);
    //3
    G4RotationMatrix* rot3_l = new G4RotationMatrix(-120*pi/180,0,0);
    new G4PVPlacement(rot3_l,G4ThreeVector(-87*mm,31*mm,0*mm),"InsideRibsWorld_Left",InsideRibs_Log,TorsoPhys,false,0);
    new G4PVPlacement(rot3_l,G4ThreeVector(-87*mm,31*mm,0*mm),"OutsideRibsWorld_Left",OutsideRibs_Log,TorsoPhys,false,0);
    //4
    G4RotationMatrix* rot4_l = new G4RotationMatrix(130*pi/180,0,0);
    new G4PVPlacement(rot4_l,G4ThreeVector(-37*mm,38*mm,0*mm),"InsideRibsWorld_Left",InsideRibs_Log,TorsoPhys,false,0);
    new G4PVPlacement(rot4_l,G4ThreeVector(-37*mm,38*mm,0*mm),"OutsideRibsWorld_Left",OutsideRibs_Log,TorsoPhys,false,0);

    //Placement Ribs Right Lung
    //1
    new G4PVPlacement(0,G4ThreeVector(60*mm,-49*mm,0*mm),"InsideRibsWorld_Right",InsideRibs_Log,TorsoPhys,false,0);
    new G4PVPlacement(0,G4ThreeVector(60*mm,-49*mm,0*mm),"OutsideRibsWorld_Right",OutsideRibs_Log,TorsoPhys,false,0);
    //2
    G4RotationMatrix* rot2_r = new G4RotationMatrix(75*pi/180,0,0);
    new G4PVPlacement(rot2_r,G4ThreeVector(95*mm,-14*mm,0*mm),"InsideRibsWorld_Right",InsideRibs_Log,TorsoPhys,false,0);
    new G4PVPlacement(rot2_r,G4ThreeVector(95*mm,-14*mm,0*mm),"OutsideRibsWorld_Right",OutsideRibs_Log,TorsoPhys,false,0);
    //3
    G4RotationMatrix* rot3_r = new G4RotationMatrix(120*pi/180,0,0);
    new G4PVPlacement(rot3_r,G4ThreeVector(87*mm,31*mm,0*mm),"InsideRibsWorld_Right",InsideRibs_Log,TorsoPhys,false,0);
    new G4PVPlacement(rot3_r,G4ThreeVector(87*mm,31*mm,0*mm),"OutsideRibsWorld_Right",OutsideRibs_Log,TorsoPhys,false,0);
    //4
    G4RotationMatrix* rot4_r = new G4RotationMatrix(-130*pi/180,0,0);
    new G4PVPlacement(rot4_r,G4ThreeVector(37*mm,38*mm,0*mm),"InsideRibsWorld_Right",InsideRibs_Log,TorsoPhys,false,0);
    new G4PVPlacement(rot4_r,G4ThreeVector(37*mm,38*mm,0*mm),"OutsideRibsWorld_Right",OutsideRibs_Log,TorsoPhys,false,0);

    //Right Lung Base
    G4EllipticalTube* RightLungBase       = new G4EllipticalTube("RightLungBase",37.5*mm,52.5*mm,150*mm);
    G4LogicalVolume* RightLungBase_Log = new G4LogicalVolume(RightLungBase,LungInflated,"RightLungBase_Log");
    RightLungBase_Log->SetVisAttributes(att2);
    new G4PVPlacement(0,G4ThreeVector(-60*mm, 0,0*mm),"RightLungBaseWorld",RightLungBase_Log,TorsoPhys,false,0);

    //Left Lung Base
    G4EllipticalTube* LeftLungBase= new G4EllipticalTube("LeftLungBase",37.5*mm,52.5*mm,150*mm);
    G4LogicalVolume* LeftLungBase_Log = new G4LogicalVolume(LeftLungBase,LungInflated,"LeftLungBase_Log");
    LeftLungBase_Log->SetVisAttributes(att2);
    new G4PVPlacement(0,G4ThreeVector(60*mm, 0,0*mm),"LeftLungBaseWorld",LeftLungBase_Log,TorsoPhys,false,0);

    /*
    //Left Lung Top
    G4Ellipsoid* LeftLungTop= new G4Ellipsoid("LeftLungTop",37.5*mm,52.5*mm,20*mm);
    G4LogicalVolume* LeftLungTop_Log = new G4LogicalVolume(LeftLungTop,LungInflated,"LeftLungTop_Log");
    LeftLungTop_Log->SetVisAttributes(att2);
    new G4PVPlacement(0,G4ThreeVector(60*mm,0,120*mm),"LeftLungTopWorld",LeftLungTop_Log,TorsoPhys,false,0);

    //Right Lung Top
    G4Ellipsoid* RightLungTop= new G4Ellipsoid("RightLungTop",37.5*mm,52.5*mm,20*mm);
    G4LogicalVolume* RightLungTop_Log = new G4LogicalVolume(RightLungTop,LungInflated,"RightLungTop_Log");
    RightLungTop_Log->SetVisAttributes(att2);
    new G4PVPlacement(0,G4ThreeVector(-60*mm,0,120*mm),"RightLungTopWorld",RightLungTop_Log,TorsoPhys,false,0);
    */
  }

  //----------------------------------------------------------------------------------------------------------------
  // Abdomen phantom (UdM Collaboration)
  //----------------------------------------------------------------------------------------------------------------

  else if(thePhantom =="Abdomen"){
    rotExt->rotateZ(pi/2.);

    G4Material *Rib2ndICRUreport46             = theMaterial->ConstructMaterial("Rib_2nd_ICRU_report_46",1.41);
    G4Material *AdiposemeanZWoodard1986        = theMaterial->ConstructMaterial("Adipose_mean_Z_Woodard_1986",0.95);
    G4Material *liver1Woodard1986              = theMaterial->ConstructMaterial("liver1_Woodard_1986",1.05);
    G4Material *brainwhiteWoodard1986          = theMaterial->ConstructMaterial("brain_white_Woodard_1986",1.04);
    G4Material *vertebralcolumnC4ICRUreport46  = theMaterial->ConstructMaterial("vertebral_column_C4_ICRU_report_46",1.42);
    G4Material *aortaWoodard1986               = theMaterial->ConstructMaterial("aorta_Woodard_1986",1.05);
    G4Material *bloodwholeWoodard1986          = theMaterial->ConstructMaterial("blood_whole_Woodard_1986",1.06);
    G4Material *bileWoodard1986                = theMaterial->ConstructMaterial("bile_Woodard_1986",1.03);
    G4Material *adrenalglandWoodard1986        = theMaterial->ConstructMaterial("adrenal_gland_Woodard_1986",1.03);
    G4Material *MusclemeanZWoodard1986         = theMaterial->ConstructMaterial("Muscle_mean_Z_Woodard_1986",1.05);
    G4Material *smallintestinewallWoodard1986  = theMaterial->ConstructMaterial("small_intestine_wall_Woodard_1986",1.03);
    G4Material *stomachWoodard1986             = theMaterial->ConstructMaterial("stomach_Woodard_1986",1.05);
    G4Material *kidney1Woodard1986             = theMaterial->ConstructMaterial("kidney1_Woodard_1986",1.05);
    G4Material *spleenWoodard1986              = theMaterial->ConstructMaterial("spleen_Woodard_1986",1.06);
    G4Material *GlandmeanZWoodard1986          = theMaterial->ConstructMaterial("Gland_mean_Z_Woodard_1986",1.02);

    G4VisAttributes* att1    = new G4VisAttributes(true,G4Colour(0,1,1));
    G4VisAttributes* att2    = new G4VisAttributes(true,G4Colour(1,0,1));
    G4VisAttributes* att3    = new G4VisAttributes(true,G4Colour(0.5,1,0));
    //G4VisAttributes* att4    = new G4VisAttributes(true,G4Colour(0,0,1));

    //Outer shape fat
    G4EllipticalTube* OuterShape= new G4EllipticalTube("OuterShape",16*cm,12*cm,10*cm);
    G4LogicalVolume* OuterShape_Log = new G4LogicalVolume(OuterShape,AdiposemeanZWoodard1986,"OuterShape_Log");
    OuterShape_Log->SetVisAttributes(att2);
    G4VPhysicalVolume* OuterShapePhys  = new G4PVPlacement(rotExt,G4ThreeVector(-0.0*cm,0*cm,0*cm),"OuterShapeWorld",OuterShape_Log,cont_phys,false,0);

    //Intercostals Muscle
    G4EllipticalTube* Intercostals      = new G4EllipticalTube("Intercostals",15*cm,11*cm,10*cm);
    G4LogicalVolume* Intercostals_Log   = new G4LogicalVolume(Intercostals,MusclemeanZWoodard1986,"Intercostals_Log");
    Intercostals_Log->SetVisAttributes(att2);
    G4RotationMatrix* Intercostals_rot  = new G4RotationMatrix(0,0,0*pi/180);
    G4ThreeVector IntercostalsShift (0*cm,-0.5*cm,0*cm);
    G4VPhysicalVolume* IntercostalsPhys = new G4PVPlacement(Intercostals_rot,G4ThreeVector(-0.0*cm,-0.5*cm,0*cm),"IntercostalsWorld",Intercostals_Log,OuterShapePhys,false,0);

    //Inner adipose Union Solid
    G4EllipticalTube* InnerAdipose= new G4EllipticalTube("InnerAdipose",12*cm,8.5*cm,10*cm);
    G4EllipticalTube* InnerAdiposeL= new G4EllipticalTube("InnerAdiposeL",5.5*cm,5*cm,10*cm);
    G4EllipticalTube* InnerAdiposeR= new G4EllipticalTube("InnerAdiposeR",5.5*cm,5*cm,10*cm);
    G4EllipticalTube* InnerAdiposeR2= new G4EllipticalTube("InnerAdiposeR2",3*cm,1.5*cm,10*cm);
    G4EllipticalTube* InnerAdiposeL2= new G4EllipticalTube("InnerAdiposeL2",3*cm,1.5*cm,10*cm);
    G4EllipticalTube* AdiposeR6= new G4EllipticalTube("AdiposeR6",2.5*cm,1*cm,10*cm);
    G4EllipticalTube* AdiposeL6= new G4EllipticalTube("AdiposeL6",2.5*cm,1*cm,10*cm);
    G4EllipticalTube* AdiposeR5= new G4EllipticalTube("AdiposeR5",1.5*cm,2.2*cm,10*cm);

    G4RotationMatrix* InnerAdiposeR_rot = new G4RotationMatrix(0,0,-35*pi/180);
    G4RotationMatrix* InnerAdiposeL_rot = new G4RotationMatrix(0,0,35*pi/180);
    G4RotationMatrix* InnerAdiposeR2_rot = new G4RotationMatrix(0,0,80*pi/180);
    G4RotationMatrix* InnerAdiposeL2_rot = new G4RotationMatrix(0,0,-80*pi/180);
    G4RotationMatrix* AdiposeR6_rot = new G4RotationMatrix(0,0,10*pi/180);
    G4RotationMatrix* AdiposeL6_rot = new G4RotationMatrix(0,0,-10*pi/180);
    G4RotationMatrix* AdiposeR5_rot = new G4RotationMatrix(0,0,-20*pi/180);
    G4ThreeVector InnerAdiposeShift (0*cm,1.0*cm,0*cm);

    G4UnionSolid* InnerAdipose_L       = new G4UnionSolid("InnerAdipose+L", InnerAdipose,InnerAdiposeL, InnerAdiposeL_rot, G4ThreeVector(7.0*cm,-2*cm,0*mm)-InnerAdiposeShift);
    G4UnionSolid* InnerAdipose_LR      = new G4UnionSolid("InnerAdipose+LR", InnerAdipose_L,InnerAdiposeR, InnerAdiposeR_rot, G4ThreeVector(-7.0*cm,-2*cm,0*mm)-InnerAdiposeShift);
    G4UnionSolid* InnerAdipose_LR_R2   = new G4UnionSolid("InnerAdipose+LR_R2", InnerAdipose_LR,InnerAdiposeR2, InnerAdiposeR2_rot, G4ThreeVector(-10.8*cm,1.6*cm,0*mm)-InnerAdiposeShift);
    G4UnionSolid* InnerAdipose_LR_LR2  = new G4UnionSolid("InnerAdipose+LR+LR2", InnerAdipose_LR_R2,InnerAdiposeL2, InnerAdiposeL2_rot, G4ThreeVector(10.8*cm,1.6*cm,0*mm)-InnerAdiposeShift);
    G4UnionSolid* InnerAdipose_LR_R6   = new G4UnionSolid("InnerAdipose+LR+R6", InnerAdipose_LR_LR2,AdiposeR6, AdiposeR6_rot, G4ThreeVector(-5.5*cm,-7*cm,0*mm)-InnerAdiposeShift);
    G4UnionSolid* InnerAdipose_LR_LR6  = new G4UnionSolid("InnerAdipose+LR+LR6", InnerAdipose_LR_R6,AdiposeL6, AdiposeL6_rot, G4ThreeVector(5.5*cm,-7*cm,0*mm)-InnerAdiposeShift);
    G4UnionSolid* InnerAdipose_R5      = new G4UnionSolid("InnerAdipose+LR+LR6+R5", InnerAdipose_LR_LR6,AdiposeR5, AdiposeR5_rot, G4ThreeVector(-6.0*cm,-4.5*cm,0*mm)-InnerAdiposeShift);

    G4LogicalVolume* InnerAdipose_Log = new G4LogicalVolume(InnerAdipose_R5,AdiposemeanZWoodard1986,"InnerAdipose_Log");
    G4VPhysicalVolume *InnerAdiposePhys =  new G4PVPlacement(0,G4ThreeVector(0.0*cm,1*cm,0*cm)-IntercostalsShift,"InnerAdiposeWorld",InnerAdipose_Log,IntercostalsPhys,false,0);

    //Liver Union Solid
    G4EllipticalTube* Liver= new G4EllipticalTube("Liver",6.2*cm,4*cm,10*cm);
    G4EllipticalTube* Liverext= new G4EllipticalTube("Liverext",2*cm,0.8*cm,10*cm);
    G4EllipticalTube* Liverext2= new G4EllipticalTube("Liverext2",3*cm,1.5*cm,10*cm);
    G4EllipticalTube* Liverext3= new G4EllipticalTube("Liverext3",2.5*cm,0.8*cm,10*cm);

    G4RotationMatrix* Liverext_rot = new G4RotationMatrix(0,0,(-20+96)*pi/180);
    G4RotationMatrix* Liverext2_rot = new G4RotationMatrix(0,0,(25+96)*pi/180);
    G4RotationMatrix* Liverext3_rot = new G4RotationMatrix(0,0,(40+96)*pi/180);

    G4ThreeVector LiverShift (-9.3*cm,-0.7*cm,0*cm);
    G4UnionSolid* LiverUnion1       = new G4UnionSolid("Liver+ext1",Liver,Liverext, Liverext_rot, G4ThreeVector(G4ThreeVector(-8.5*cm,-6.2*cm,0*cm)-LiverShift).rotateZ(-96*pi/180));
    G4UnionSolid* LiverUnion2       = new G4UnionSolid("Liver+ext1+ext2",LiverUnion1,Liverext2, Liverext2_rot, G4ThreeVector(G4ThreeVector(-7.0*cm,-1.6*cm,0*cm)-LiverShift).rotateZ(-96*pi/180));
    G4UnionSolid* LiverUnion3       = new G4UnionSolid("Liver+ext1+ext2+ext3",LiverUnion2,Liverext3, Liverext3_rot, G4ThreeVector(G4ThreeVector(-10.3*cm,4*cm,0*cm)-LiverShift).rotateZ(-96*pi/180));

    //Liver Subtraction Solid
    G4RotationMatrix* Adipose_rot = new G4RotationMatrix(0,0,(-20-96)*pi/180);
    G4SubtractionSolid* LiverSub1 = new G4SubtractionSolid("Liver-adipose",LiverUnion3,AdiposeR5,Adipose_rot,G4ThreeVector(G4ThreeVector(-6.0*cm,-4.5*cm,0*mm)-LiverShift).rotateZ(-96*pi/180));
    G4LogicalVolume* Liver_Log = new G4LogicalVolume(LiverSub1,liver1Woodard1986,"Liver_Log");
    Liver_Log->SetVisAttributes(att3);
    G4RotationMatrix* Liver_rot = new G4RotationMatrix(0,0,96*pi/180);
    G4VPhysicalVolume *LiverPhys  = new G4PVPlacement(Liver_rot,G4ThreeVector(-9.3*cm,-0.7*cm,0*cm)-InnerAdiposeShift,"LiverWorld",Liver_Log,InnerAdiposePhys,false,0);

    // Vertebral Column Union Solid
    G4EllipticalTube* VertebralColumn      = new G4EllipticalTube("VertebralColumn",2.1*cm,2.1*cm,10*cm);
    G4EllipticalTube* VertebralColumnback  = new G4EllipticalTube("VertebralColumnback",0.4*cm,1.2*cm,10*cm);
    G4EllipticalTube* VertebralColumnRback = new G4EllipticalTube("VertebralColumnRback",0.3*cm,0.9*cm,10*cm);
    G4EllipticalTube* VertebralColumnLback = new G4EllipticalTube("VertebralColumnLback",0.3*cm,0.9*cm,10*cm);

    G4ThreeVector VertebralColumnShift (0*cm,-6.4*cm,0*cm);
    G4RotationMatrix* VertebralColumnLback_rot = new G4RotationMatrix(0,0,75*pi/180);
    G4RotationMatrix* VertebralColumnRback_rot = new G4RotationMatrix(0,0,-75*pi/180);

    G4UnionSolid* VertebralColumn_back   = new G4UnionSolid("VertebralColumn+Back",VertebralColumn,VertebralColumnback, 0, G4ThreeVector(0.0*cm,-9.6*cm,0*mm)-VertebralColumnShift);
    G4UnionSolid* VertebralColumn_backL  = new G4UnionSolid("VertebralColumn+Back+L",VertebralColumn_back,VertebralColumnLback,VertebralColumnLback_rot,G4ThreeVector(-0.9*cm,-8.5*cm,0)-VertebralColumnShift);
    G4UnionSolid* VertebralColumn_backLR = new G4UnionSolid("VertebralColumn+Back+LR",VertebralColumn_backL,VertebralColumnRback,VertebralColumnRback_rot,G4ThreeVector(0.9*cm,-8.5*cm,0)-VertebralColumnShift);

    G4LogicalVolume* VertebralColumn_Log  = new G4LogicalVolume(VertebralColumn_backLR,vertebralcolumnC4ICRUreport46,"VertebralColumn_Log");
    VertebralColumn_Log->SetVisAttributes(att2);
    G4VPhysicalVolume *VertebralColumnPhys = new G4PVPlacement(0,G4ThreeVector(-0.0*cm,-6.4*cm,0*cm)-IntercostalsShift,"VertebralColumnWorld",VertebralColumn_Log,IntercostalsPhys,false,0);

    G4EllipticalTube* Spine= new G4EllipticalTube("Spine",1*cm,1*cm,10*cm);
    G4LogicalVolume* Spine_Log = new G4LogicalVolume(Spine,brainwhiteWoodard1986,"Spine_Log");
    Spine_Log->SetVisAttributes(att2);
    G4RotationMatrix* Spine_rot = new G4RotationMatrix(0,0,0*pi/180);
    new G4PVPlacement(Spine_rot,G4ThreeVector(-0.0*cm,-7.6*cm,0*cm)-VertebralColumnShift,"SpineWorld",Spine_Log,VertebralColumnPhys,false,0);

    G4EllipticalTube* RibR1= new G4EllipticalTube("RibR1",2*cm,0.5*cm,10*cm);
    G4LogicalVolume* RibR1_Log = new G4LogicalVolume(RibR1,Rib2ndICRUreport46,"RibR1_Log");
    RibR1_Log->SetVisAttributes(att2);
    G4RotationMatrix* RibR1_rot = new G4RotationMatrix(0,0,18*pi/180);
    new G4PVPlacement(RibR1_rot,G4ThreeVector(-3.7*cm,-8*cm,0*cm)-IntercostalsShift,"RibR1World",RibR1_Log,IntercostalsPhys,false,0);

    G4EllipticalTube* RibL1= new G4EllipticalTube("RibL1",2*cm,0.5*cm,10*cm);
    G4LogicalVolume* RibL1_Log = new G4LogicalVolume(RibL1,Rib2ndICRUreport46,"RibL1_Log");
    RibL1_Log->SetVisAttributes(att2);
    G4RotationMatrix* RibL1_rot = new G4RotationMatrix(0,0,-18*pi/180);
    new G4PVPlacement(RibL1_rot,G4ThreeVector(3.7*cm,-8*cm,0*cm)-IntercostalsShift,"RibL1World",RibL1_Log,IntercostalsPhys,false,0);

    G4EllipticalTube* RibR2= new G4EllipticalTube("RibR2",1.5*cm,0.5*cm,10*cm);
    G4LogicalVolume* RibR2_Log = new G4LogicalVolume(RibR2,Rib2ndICRUreport46,"RibR2_Log");
    RibR2_Log->SetVisAttributes(att2);
    G4RotationMatrix* RibR2_rot = new G4RotationMatrix(0,0,-25*pi/180);
    new G4PVPlacement(RibR2_rot,G4ThreeVector(-8.8*cm,-7.7*cm,0*cm)-IntercostalsShift,"RibR2World",RibR2_Log,IntercostalsPhys,false,0);

    G4EllipticalTube* RibL2= new G4EllipticalTube("RibL2",1.5*cm,0.5*cm,10*cm);
    G4LogicalVolume* RibL2_Log = new G4LogicalVolume(RibL2,Rib2ndICRUreport46,"RibL2_Log");
    RibL2_Log->SetVisAttributes(att2);
    G4RotationMatrix* RibL2_rot = new G4RotationMatrix(0,0,25*pi/180);
    new G4PVPlacement(RibL2_rot,G4ThreeVector(8.8*cm,-7.7*cm,0*cm)-IntercostalsShift,"RibL2World",RibL2_Log,IntercostalsPhys,false,0);

    G4EllipticalTube* Aorta= new G4EllipticalTube("Aorta",1*cm,1*cm,10*cm);
    G4LogicalVolume* Aorta_Log = new G4LogicalVolume(Aorta,aortaWoodard1986,"Aorta_Log");
    Aorta_Log->SetVisAttributes(att2);
    G4RotationMatrix* Aorta_rot = new G4RotationMatrix(0,0,0*pi/180);
    G4VPhysicalVolume* AortaPhys  =new G4PVPlacement(Aorta_rot,G4ThreeVector(-1.0*cm,-3.5*cm,0*cm)-InnerAdiposeShift,"AortaWorld",Aorta_Log,InnerAdiposePhys,false,0);

    G4EllipticalTube* AortaBlood= new G4EllipticalTube("AortaBlood",0.85*cm,0.85*cm,10*cm);
    G4LogicalVolume* AortaBlood_Log = new G4LogicalVolume(AortaBlood,bloodwholeWoodard1986,"AortaBlood_Log");
    AortaBlood_Log->SetVisAttributes(att2);
    G4RotationMatrix* AortaBlood_rot = new G4RotationMatrix(0,0,0*pi/180);
    new G4PVPlacement(AortaBlood_rot,G4ThreeVector(0*cm,0*cm,0*cm),"AortaBloodWorld",AortaBlood_Log,AortaPhys,false,0);

    G4EllipticalTube* RibR3= new G4EllipticalTube("RibR3",1.5*cm,0.5*cm,10*cm);
    G4LogicalVolume* RibR3_Log = new G4LogicalVolume(RibR3,Rib2ndICRUreport46,"RibR3_Log");
    RibR3_Log->SetVisAttributes(att1);
    G4RotationMatrix* RibR3_rot = new G4RotationMatrix(0,0,-70*pi/180);
    new G4PVPlacement(RibR3_rot,G4ThreeVector(-13.2*cm,-3.5*cm,0*cm)-IntercostalsShift,"RibR3World",RibR3_Log,IntercostalsPhys,false,0);

    G4EllipticalTube* RibL3= new G4EllipticalTube("RibL3",1.5*cm,0.5*cm,10*cm);
    G4LogicalVolume* RibL3_Log = new G4LogicalVolume(RibL3,Rib2ndICRUreport46,"RibL3_Log");
    RibL3_Log->SetVisAttributes(att2);
    G4RotationMatrix* RibL3_rot = new G4RotationMatrix(0,0,70*pi/180);
    new G4PVPlacement(RibL3_rot,G4ThreeVector(13.2*cm,-3.5*cm,0*cm)-IntercostalsShift,"RibL3World",RibL3_Log,IntercostalsPhys,false,0);

    G4EllipticalTube* RibR4= new G4EllipticalTube("RibR4",1*cm,0.4*cm,10*cm);
    G4LogicalVolume* RibR4_Log = new G4LogicalVolume(RibR4,Rib2ndICRUreport46,"RibR4_Log");
    RibR4_Log->SetVisAttributes(att2);
    G4RotationMatrix* RibR4_rot = new G4RotationMatrix(0,0,90*pi/180);
    new G4PVPlacement(RibR4_rot,G4ThreeVector(-13.5*cm,1*cm,0*cm)-IntercostalsShift,"RibR4World",RibR4_Log,IntercostalsPhys,false,0);

    G4EllipticalTube* RibL4= new G4EllipticalTube("RibL4",1*cm,0.4*cm,10*cm);
    G4LogicalVolume* RibL4_Log = new G4LogicalVolume(RibL4,Rib2ndICRUreport46,"RibL4_Log");
    RibL4_Log->SetVisAttributes(att2);
    G4RotationMatrix* RibL4_rot = new G4RotationMatrix(0,0,90*pi/180);
    new G4PVPlacement(RibL4_rot,G4ThreeVector(13.5*cm,1*cm,0*cm)-IntercostalsShift,"RibL4World",RibL4_Log,IntercostalsPhys,false,0);

    G4EllipticalTube* RibR5= new G4EllipticalTube("RibR5",1*cm,0.4*cm,10*cm);
    G4LogicalVolume* RibR5_Log = new G4LogicalVolume(RibR5,Rib2ndICRUreport46,"RibR5_Log");
    RibR5_Log->SetVisAttributes(att2);
    G4RotationMatrix* RibR5_rot = new G4RotationMatrix(0,0,-125*pi/180);
    new G4PVPlacement(RibR5_rot,G4ThreeVector(-12.0*cm,4.5*cm,0*cm)-IntercostalsShift,"RibR5World",RibR5_Log,IntercostalsPhys,false,0);

    G4EllipticalTube* RibL5= new G4EllipticalTube("RibL5",1*cm,0.4*cm,10*cm);
    G4LogicalVolume* RibL5_Log = new G4LogicalVolume(RibL5,Rib2ndICRUreport46,"RibL5_Log");
    RibL5_Log->SetVisAttributes(att2);
    G4RotationMatrix* RibL5_rot = new G4RotationMatrix(0,0,125*pi/180);
    new G4PVPlacement(RibL5_rot,G4ThreeVector(12.0*cm,4.5*cm,0*cm)-IntercostalsShift,"RibL5World",RibL5_Log,IntercostalsPhys,false,0);

    G4EllipticalTube* RibR6= new G4EllipticalTube("RibR6",1*cm,0.4*cm,10*cm);
    G4LogicalVolume* RibR6_Log = new G4LogicalVolume(RibR6,Rib2ndICRUreport46,"RibR6_Log");
    RibR6_Log->SetVisAttributes(att2);
    G4RotationMatrix* RibR6_rot = new G4RotationMatrix(0,0,-165*pi/180);
    new G4PVPlacement(RibR6_rot,G4ThreeVector(-8.5*cm,8*cm,0*cm)-IntercostalsShift,"RibR6World",RibR6_Log,IntercostalsPhys,false,0);

    G4EllipticalTube* RibL6= new G4EllipticalTube("RibL6",1*cm,0.4*cm,10*cm);
    G4LogicalVolume* RibL6_Log = new G4LogicalVolume(RibL6,Rib2ndICRUreport46,"RibL6_Log");
    RibL6_Log->SetVisAttributes(att2);
    G4RotationMatrix* RibL6_rot = new G4RotationMatrix(0,0,165*pi/180);
    new G4PVPlacement(RibL6_rot,G4ThreeVector(8.5*cm,8*cm,0*cm)-IntercostalsShift,"RibL6World",RibL6_Log,IntercostalsPhys,false,0);

    G4EllipticalTube* Spleen= new G4EllipticalTube("Spleen",3*cm,1.5*cm,10*cm);
    G4LogicalVolume* Spleen_Log = new G4LogicalVolume(Spleen,spleenWoodard1986,"Spleen_Log");
    Spleen_Log->SetVisAttributes(att2);
    G4RotationMatrix* Spleen_rot = new G4RotationMatrix(0,0,40*pi/180);
    new G4PVPlacement(Spleen_rot,G4ThreeVector(9.7*cm,-5*cm,0*cm)-IntercostalsShift,"SpleenWorld",Spleen_Log,IntercostalsPhys,false,0);

    G4EllipticalTube* AdrenalGland= new G4EllipticalTube("AdrenalGland",1*cm,0.7*cm,10*cm);
    G4LogicalVolume* AdrenalGland_Log = new G4LogicalVolume(AdrenalGland,adrenalglandWoodard1986,"AdrenalGland_Log");
    AdrenalGland_Log->SetVisAttributes(att2);
    G4RotationMatrix* AdrenalGland_rot = new G4RotationMatrix(0,0,30*pi/180);
    new G4PVPlacement(AdrenalGland_rot,G4ThreeVector(2.3*cm,-2.3*cm,0*cm)-InnerAdiposeShift,"AdrenalGlandWorld",AdrenalGland_Log,InnerAdiposePhys,false,0);

    G4EllipticalTube* Pancreas= new G4EllipticalTube("Pancreas",5*cm,1*cm,10*cm);
    G4LogicalVolume* Pancreas_Log = new G4LogicalVolume(Pancreas,GlandmeanZWoodard1986,"Pancreas_Log");
    Pancreas_Log->SetVisAttributes(att2);
    G4RotationMatrix* Pancreas_rot = new G4RotationMatrix(0,0,-10*pi/180);
    new G4PVPlacement(Pancreas_rot,G4ThreeVector(2.5*cm,0.5*cm,0*cm)-InnerAdiposeShift,"PancreasWorld",Pancreas_Log,InnerAdiposePhys,false,0);

    G4EllipticalTube* Gallbladder= new G4EllipticalTube("Gallbladder",2.5*cm,1.5*cm,10*cm);
    G4LogicalVolume* Gallbladder_Log = new G4LogicalVolume(Gallbladder,bileWoodard1986,"Gallbladder_Log");
    Gallbladder_Log->SetVisAttributes(att2);
    G4RotationMatrix* Gallbladder_rot = new G4RotationMatrix(0,0,(-55-96)*pi/180);
    new G4PVPlacement(Gallbladder_rot,G4ThreeVector(G4ThreeVector(-7.0*cm,2.5*cm,0*cm)-LiverShift).rotateZ(-96*pi/180),"GallbladderWorld",Gallbladder_Log,LiverPhys,false,0);

    G4EllipticalTube* Colon1      = new G4EllipticalTube("Colon1",1.4*cm,1.4*cm,10*cm);
    G4LogicalVolume* Colon1_Log   = new G4LogicalVolume(Colon1,stomachWoodard1986,"Colon1_Log");
    Colon1_Log->SetVisAttributes(att2);
    G4RotationMatrix* Colon1_rot  = new G4RotationMatrix(0,0,0*pi/180);
    G4VPhysicalVolume *Colon1Phys = new G4PVPlacement(Colon1_rot,G4ThreeVector(-2.0*cm,7.1*cm,0*cm)-InnerAdiposeShift,"Colon1World",Colon1_Log,InnerAdiposePhys,false,0);

    G4EllipticalTube* Colon2= new G4EllipticalTube("Colon2",1.3*cm,1.3*cm,10*cm);
    G4LogicalVolume* Colon2_Log = new G4LogicalVolume(Colon2,stomachWoodard1986,"Colon2_Log");
    Colon2_Log->SetVisAttributes(att2);
    G4RotationMatrix* Colon2_rot = new G4RotationMatrix(0,0,0*pi/180);
    G4VPhysicalVolume *Colon2Phys = new G4PVPlacement(Colon2_rot,G4ThreeVector(-0.3*cm,6.8*cm,0*cm)-InnerAdiposeShift,"Colon2World",Colon2_Log,InnerAdiposePhys,false,0);

    G4EllipticalTube* Colon3= new G4EllipticalTube("Colon3",1*cm,1.4*cm,10*cm);
    G4LogicalVolume* Colon3_Log = new G4LogicalVolume(Colon3,stomachWoodard1986,"Colon3_Log");
    Colon3_Log->SetVisAttributes(att2);
    G4RotationMatrix* Colon3_rot = new G4RotationMatrix(0,0,0*pi/180);
    G4VPhysicalVolume *Colon3Phys = new G4PVPlacement(Colon3_rot,G4ThreeVector(1.0*cm,7.2*cm,0*cm)-InnerAdiposeShift,"Colon3World",Colon3_Log,InnerAdiposePhys,false,0);

    G4EllipticalTube* Colon4= new G4EllipticalTube("Colon4",1.3*cm,1*cm,10*cm);
    G4LogicalVolume* Colon4_Log = new G4LogicalVolume(Colon4,stomachWoodard1986,"Colon4_Log");
    Colon4_Log->SetVisAttributes(att2);
    G4RotationMatrix* Colon4_rot = new G4RotationMatrix(0,0,0*pi/180);
    G4VPhysicalVolume *Colon4Phys = new G4PVPlacement(Colon4_rot,G4ThreeVector(2.0*cm,6.6*cm,0*cm)-InnerAdiposeShift,"Colon4World",Colon4_Log,InnerAdiposePhys,false,0);

    G4EllipticalTube* Colon5= new G4EllipticalTube("Colon5",1.3*cm,1.3*cm,10*cm);
    G4LogicalVolume* Colon5_Log = new G4LogicalVolume(Colon5,stomachWoodard1986,"Colon5_Log");
    Colon5_Log->SetVisAttributes(att2);
    G4RotationMatrix* Colon5_rot = new G4RotationMatrix(0,0,0*pi/180);
    G4VPhysicalVolume *Colon5Phys = new G4PVPlacement(Colon5_rot,G4ThreeVector(3.6*cm,6.2*cm,0*cm)-InnerAdiposeShift,"Colon5World",Colon5_Log,InnerAdiposePhys,false,0);

    G4EllipticalTube* ColonAir1= new G4EllipticalTube("ColonAir1",1*cm,1*cm,10*cm);
    G4LogicalVolume* ColonAir1_Log = new G4LogicalVolume(ColonAir1,air,"ColonAir1_Log");
    ColonAir1_Log->SetVisAttributes(att2);
    G4RotationMatrix* ColonAir1_rot = new G4RotationMatrix(0,0,0*pi/180);
    new G4PVPlacement(ColonAir1_rot,G4ThreeVector(),"ColonAir1World",ColonAir1_Log,Colon1Phys,false,0);

    G4EllipticalTube* ColonAir2= new G4EllipticalTube("ColonAir2",0.8*cm,0.8*cm,10*cm);
    G4LogicalVolume* ColonAir2_Log = new G4LogicalVolume(ColonAir2,air,"ColonAir2_Log");
    ColonAir2_Log->SetVisAttributes(att2);
    G4RotationMatrix* ColonAir2_rot = new G4RotationMatrix(0,0,0*pi/180);
    new G4PVPlacement(ColonAir2_rot,G4ThreeVector(),"ColonAir2World",ColonAir2_Log,Colon2Phys,false,0);

    G4EllipticalTube* ColonAir3= new G4EllipticalTube("ColonAir3",0.9*cm,0.9*cm,10*cm);
    G4LogicalVolume* ColonAir3_Log = new G4LogicalVolume(ColonAir3,air,"ColonAir3_Log");
    ColonAir3_Log->SetVisAttributes(att2);
    G4RotationMatrix* ColonAir3_rot = new G4RotationMatrix(0,0,0*pi/180);
    new G4PVPlacement(ColonAir3_rot,G4ThreeVector(),"ColonAir3World",ColonAir3_Log,Colon3Phys,false,0);

    G4EllipticalTube* ColonAir4= new G4EllipticalTube("ColonAir4",0.6*cm,0.6*cm,10*cm);
    G4LogicalVolume* ColonAir4_Log = new G4LogicalVolume(ColonAir4,air,"ColonAir4_Log");
    ColonAir4_Log->SetVisAttributes(att2);
    G4RotationMatrix* ColonAir4_rot = new G4RotationMatrix(0,0,0*pi/180);
    new G4PVPlacement(ColonAir4_rot,G4ThreeVector(),"ColonAir4World",ColonAir4_Log,Colon4Phys,false,0);

    G4EllipticalTube* ColonAir5= new G4EllipticalTube("ColonAir5",0.7*cm,0.7*cm,10*cm);
    G4LogicalVolume* ColonAir5_Log = new G4LogicalVolume(ColonAir5,air,"ColonAir5_Log");
    ColonAir5_Log->SetVisAttributes(att2);
    G4RotationMatrix* ColonAir5_rot = new G4RotationMatrix(0,0,0*pi/180);
    new G4PVPlacement(ColonAir5_rot,G4ThreeVector(),"ColonAir5World",ColonAir5_Log,Colon5Phys,false,0);

    G4EllipticalTube* SmallIntestine= new G4EllipticalTube("SmallIntestine",1.3*cm,2.2*cm,10*cm);
    G4LogicalVolume* SmallIntestine_Log = new G4LogicalVolume(SmallIntestine,smallintestinewallWoodard1986,"SmallIntestine_Log");
    SmallIntestine_Log->SetVisAttributes(att2);
    G4RotationMatrix* SmallIntestine_rot = new G4RotationMatrix(0,0,30*pi/180);
    G4ThreeVector SmallIntestineShift = G4ThreeVector(9.0*cm,2*cm,0*cm);
    G4VPhysicalVolume* SmallIntestinePhys = new G4PVPlacement(SmallIntestine_rot,G4ThreeVector(9.0*cm,2*cm,0*cm)-InnerAdiposeShift,"SmallIntestineWorld",SmallIntestine_Log,InnerAdiposePhys,false,0);

    G4EllipticalTube* SmallIntestine1= new G4EllipticalTube("SmallIntestine1",0.8*cm,1.3*cm,10*cm);
    G4LogicalVolume* SmallIntestine1_Log = new G4LogicalVolume(SmallIntestine1,smallintestinewallWoodard1986,"SmallIntestine1_Log");
    SmallIntestine1_Log->SetVisAttributes(att2);
    G4RotationMatrix* SmallIntestine1_rot = new G4RotationMatrix(0,0,-20*pi/180);
    G4ThreeVector SmallIntestine1Shift = G4ThreeVector(7.3*cm,3.3*cm,0*cm);
    G4VPhysicalVolume* SmallIntestine1Phys = new G4PVPlacement(SmallIntestine1_rot,G4ThreeVector(7.3*cm,3.3*cm,0*cm)-InnerAdiposeShift,"SmallIntestine1World",SmallIntestine1_Log,InnerAdiposePhys,false,0);

    G4EllipticalTube* SmallIntestine2= new G4EllipticalTube("SmallIntestine2",0.7*cm,1.1*cm,10*cm);
    G4LogicalVolume* SmallIntestine2_Log = new G4LogicalVolume(SmallIntestine2,smallintestinewallWoodard1986,"SmallIntestine2_Log");
    SmallIntestine2_Log->SetVisAttributes(att2);
    G4RotationMatrix* SmallIntestine2_rot = new G4RotationMatrix(0,0,35*pi/180);
    new G4PVPlacement(SmallIntestine2_rot,G4ThreeVector(7.8*cm,1.2*cm,0*cm)-InnerAdiposeShift,"SmallIntestine2World",SmallIntestine2_Log,InnerAdiposePhys,false,0);

    G4EllipticalTube* SmallIntestineAir= new G4EllipticalTube("SmallIntestineAir",0.8*cm,0.9*cm,10*cm);
    G4LogicalVolume* SmallIntestineAir_Log = new G4LogicalVolume(SmallIntestineAir,air,"SmallIntestineAir_Log");
    SmallIntestineAir_Log->SetVisAttributes(att1);
    G4Transform3D SmallIntestineAirTransform = G4Transform3D(G4RotationMatrix(0,0,30*pi/180), G4ThreeVector(G4ThreeVector(8.8*cm,2.8*cm,0*cm)-SmallIntestineShift).rotateZ(-30*pi/180) );
    new G4PVPlacement(SmallIntestineAirTransform,"SmallIntestineAirWorld",SmallIntestineAir_Log,SmallIntestinePhys,false,0);

    G4EllipticalTube* SmallIntestineAir1= new G4EllipticalTube("SmallIntestineAir1",0.6*cm,0.5*cm,10*cm);
    G4LogicalVolume* SmallIntestineAir1_Log = new G4LogicalVolume(SmallIntestineAir1,air,"SmallIntestineAir1_Log");
    SmallIntestineAir1_Log->SetVisAttributes(att1);
    G4Transform3D SmallIntestineAir1Transform = G4Transform3D(G4RotationMatrix(0,0,0*pi/180), G4ThreeVector(G4ThreeVector(7.6*cm,3.8*cm,0*cm)-SmallIntestine1Shift).rotateZ(20*pi/180) );
    new G4PVPlacement(SmallIntestineAir1Transform,"SmallIntestineAir1World",SmallIntestineAir1_Log,SmallIntestine1Phys,false,0);

    G4EllipticalTube* KidneyR= new G4EllipticalTube("KidneyR",1.5*cm,2*cm,10*cm);
    G4LogicalVolume* KidneyR_Log = new G4LogicalVolume(KidneyR,kidney1Woodard1986,"KidneyR_Log");
    KidneyR_Log->SetVisAttributes(att2);
    G4RotationMatrix* KidneyR_rot = new G4RotationMatrix(0,0,-20*pi/180);
    new G4PVPlacement(KidneyR_rot,G4ThreeVector(-4.2*cm,-5*cm,0*cm)-InnerAdiposeShift,"KidneyRWorld",KidneyR_Log,InnerAdiposePhys,false,0);

    G4EllipticalTube* KidneyL= new G4EllipticalTube("KidneyL",1.8*cm,2.3*cm,10*cm);
    G4LogicalVolume* KidneyL_Log = new G4LogicalVolume(KidneyL,kidney1Woodard1986,"KidneyL_Log");
    KidneyL_Log->SetVisAttributes(att2);
    G4RotationMatrix* KidneyL_rot = new G4RotationMatrix(0,0,30*pi/180);
    new G4PVPlacement(KidneyL_rot,G4ThreeVector(4.2*cm,-5*cm,0*cm)-InnerAdiposeShift,"KidneyLWorld",KidneyL_Log,InnerAdiposePhys,false,0);

  }
//----------------------------------------------------------------------------------------------------------------
// UDM phantom (UdM Collaboration)
//----------------------------------------------------------------------------------------------------------------

  else if(thePhantom =="HeadUDM"){

    rotExt->rotateZ(pi/2.);
    G4Material *AdiposemeanZWoodard1986        = theMaterial->ConstructMaterial("Adipose_mean_Z_Woodard_1986",0.95);
    G4Material *brainwhiteWoodard1986          = theMaterial->ConstructMaterial("brain_white_Woodard_1986",1.04);
    G4Material *craniumICRUreport46            = theMaterial->ConstructMaterial("Cranium_ICRU_report_46",1.61);
    G4Material *braingreyWoodard1986           = theMaterial->ConstructMaterial("brain_greymatter_Woodard_1986",1.04);
    G4Material *corticalboneICRUreport46       = theMaterial->ConstructMaterial("cortical_bone_ICRU_Report_46",1.92);

    G4VisAttributes* att2    = new G4VisAttributes(true,G4Colour(1,0,1));

    //Outline fat
    G4EllipticalTube* Outline= new G4EllipticalTube("Outline",7*cm,9*cm,10*cm);
    G4LogicalVolume* Outline_Log = new G4LogicalVolume(Outline,AdiposemeanZWoodard1986,"Outline_Log");
    Outline_Log->SetVisAttributes(att2);
    G4VPhysicalVolume* Outline_phys = new G4PVPlacement(rotExt,G4ThreeVector(-0.0*cm,0*cm,0*cm),"OutlineWorld",Outline_Log,cont_phys,false,0);

    // Cranium full union
    G4EllipticalTube* Cranium= new G4EllipticalTube("Cranium",6*cm,6*cm,10*cm);
    G4EllipticalTube* CraniumRight= new G4EllipticalTube("CraniumRight",3*cm,3.5*cm,10*cm);
    G4EllipticalTube* CraniumLeft= new G4EllipticalTube("CraniumLeft",3*cm,3.5*cm,10*cm);
    G4EllipticalTube* CraniumRight4= new G4EllipticalTube("CraniumRight4",0.7*cm,2*cm,10*cm);
    G4EllipticalTube* CraniumLeft4= new G4EllipticalTube("CraniumLeft4",0.7*cm,2*cm,10*cm);

    G4RotationMatrix* Cranium_rot = new G4RotationMatrix(0,0,0*pi/180);
    G4RotationMatrix* CraniumRight_rot = new G4RotationMatrix(0,0,30*pi/180);
    G4RotationMatrix* CraniumLeft_rot = new G4RotationMatrix(0,0,330*pi/180);
    G4RotationMatrix* CraniumRight4_rot = new G4RotationMatrix(0,0,0*pi/180);
    G4RotationMatrix* CraniumLeft4_rot = new G4RotationMatrix(0,0,0*pi/180);

    G4ThreeVector Cranium_shift = G4ThreeVector(-0.0*cm,-2.2*cm,0*cm);

    G4UnionSolid* CraniumUnion1    = new G4UnionSolid("Cranium+CraniumRight",Cranium,CraniumRight, CraniumRight_rot, G4ThreeVector(G4ThreeVector(3.0*cm,1*cm,0*cm)-Cranium_shift));
    G4UnionSolid* CraniumUnion2    = new G4UnionSolid("CraniumUnion1+CraniumLeft",CraniumUnion1,CraniumLeft, CraniumLeft_rot, G4ThreeVector(G4ThreeVector(-3.0*cm,1*cm,0*cm)-Cranium_shift));
    G4UnionSolid* CraniumUnion3    = new G4UnionSolid("CraniumUnion2+CraniumRight4",CraniumUnion2,CraniumRight4, CraniumRight4_rot, G4ThreeVector(G4ThreeVector(5.8*cm,-1*cm,0*cm)-Cranium_shift));
    G4UnionSolid* CraniumUnion4    = new G4UnionSolid("CraniumUnion3+CraniumLeft4",CraniumUnion3,CraniumLeft4, CraniumLeft4_rot, G4ThreeVector(G4ThreeVector(-5.8*cm,-1*cm,0*cm)-Cranium_shift));
    G4LogicalVolume* Cranium_Log = new G4LogicalVolume(CraniumUnion4,craniumICRUreport46,"Cranium_Log");
    Cranium_Log->SetVisAttributes(att2);
    G4VPhysicalVolume* Cranium_phys = new G4PVPlacement(Cranium_rot,G4ThreeVector(-0.0*cm,-2.2*cm,0*cm),"CraniumWorld",Cranium_Log,Outline_phys,false,0);

    // Brain right (grey matter)
    G4EllipticalTube* BrainRight= new G4EllipticalTube("BrainRight",2.2*cm,3*cm,10*cm);
    G4LogicalVolume* BrainRight_Log = new G4LogicalVolume(BrainRight,braingreyWoodard1986,"BrainRight_Log");
    BrainRight_Log->SetVisAttributes(att2);
    G4RotationMatrix* BrainRight_rot = new G4RotationMatrix(0,0,40*pi/180);
    G4VPhysicalVolume* BrainRight_phys = new G4PVPlacement(BrainRight_rot,G4ThreeVector(3.1*cm,1.3*cm,0*cm)-Cranium_shift,"BrainRightWorld",BrainRight_Log,Cranium_phys,false,0);
    G4ThreeVector BrainRight_shift = G4ThreeVector(3.1*cm,1.3*cm,0*cm);

    // Brain left (grey matter)
    G4EllipticalTube* BrainLeft= new G4EllipticalTube("BrainLeft",2.2*cm,3*cm,10*cm);
    G4LogicalVolume* BrainLeft_Log = new G4LogicalVolume(BrainLeft,braingreyWoodard1986,"BrainLeft_Log");
    BrainLeft_Log->SetVisAttributes(att2);
    G4RotationMatrix* BrainLeft_rot = new G4RotationMatrix(0,0,320*pi/180);
    G4VPhysicalVolume* BrainLeft_phys = new G4PVPlacement(BrainLeft_rot,G4ThreeVector(-3.1*cm,1.3*cm,0*cm)-Cranium_shift,"BrainLeftWorld",BrainLeft_Log,Cranium_phys,false,0);
    G4ThreeVector BrainLeft_shift = G4ThreeVector(-3.1*cm,1.3*cm,0*cm);

    // Brain left (white matter)
    G4EllipticalTube* Brain3= new G4EllipticalTube("Brain3",5*cm,3*cm,10*cm);
    G4LogicalVolume* Brain3_Log = new G4LogicalVolume(Brain3,brainwhiteWoodard1986,"Brain3_Log");
    Brain3_Log->SetVisAttributes(att2);
    G4RotationMatrix* Brain3_rot = new G4RotationMatrix(0,0,0*pi/180);
    new G4PVPlacement(Brain3_rot,G4ThreeVector(-0.0*cm,-4*cm,0*cm) - Cranium_shift,"Brain3World",Brain3_Log,Cranium_phys,false,0);

    //  Brain right (white matter)
    G4EllipticalTube* Brain2= new G4EllipticalTube("Brain2",2.8*cm,4.5*cm,10*cm);
    G4LogicalVolume* Brain2_Log = new G4LogicalVolume(Brain2,brainwhiteWoodard1986,"Brain2_Log");
    Brain2_Log->SetVisAttributes(att2);
    G4RotationMatrix* Brain2_rot = new G4RotationMatrix(0,0,30*pi/180 - 40*pi/180);
    G4ThreeVector Brain2_shift = G4ThreeVector( G4ThreeVector(2.3*cm,-1.2*cm,0*cm) - BrainRight_shift).rotateZ(-40*pi/180);
    G4VPhysicalVolume* Brain2_phys = new G4PVPlacement(Brain2_rot, Brain2_shift, "Brain2World",Brain2_Log,BrainRight_phys,false,0);

    // Brain center (white matter)
    G4EllipticalTube* Brain1= new G4EllipticalTube("Brain1",2.8*cm,4.5*cm,10*cm);
    G4LogicalVolume* Brain1_Log = new G4LogicalVolume(Brain1,brainwhiteWoodard1986,"Brain1_Log");
    Brain1_Log->SetVisAttributes(att2);
    G4RotationMatrix* Brain1_rot  = new G4RotationMatrix(0,0,(330-320)*pi/180);
    G4ThreeVector Brain1_shift = G4ThreeVector( G4ThreeVector(-2.3*cm,-1.2*cm,0*cm) - BrainLeft_shift).rotateZ(-320*pi/180);
    G4VPhysicalVolume* Brain1_phys = new G4PVPlacement(Brain1_rot, Brain1_shift, "Brain1World",Brain1_Log,BrainLeft_phys,false,0);

    //Cranium Right 3
    G4EllipticalTube* CraniumRight3= new G4EllipticalTube("CraniumRight3",0.5*cm,1*cm,10*cm);
    G4EllipticalTube* CraniumRight2= new G4EllipticalTube("CraniumRight2",0.4*cm,2.5*cm,10*cm);
    G4RotationMatrix* CraniumRight2_rot = new G4RotationMatrix(0,0,(60-80)*pi/180);
    G4UnionSolid* CraniumRight3_Union1 = new G4UnionSolid("CraniumRight3+CraniumRight2",CraniumRight3,CraniumRight2,CraniumRight2_rot,G4ThreeVector(G4ThreeVector(3.6*cm,0.1*cm,0*cm)
															    -G4ThreeVector(4.9*cm,-0.2*cm,0*cm)).rotateZ(-80*pi/180));
    G4LogicalVolume* CraniumRight3_Log = new G4LogicalVolume(CraniumRight3_Union1,craniumICRUreport46,"CraniumRight3_Log");
    CraniumRight3_Log->SetVisAttributes(att2);
    G4RotationMatrix* CraniumRight3_rot = new G4RotationMatrix(0,0,(80-30)*pi/180);
    G4ThreeVector CraniumRight3_shift   = G4ThreeVector( G4ThreeVector(4.9*cm,-0.2*cm,0*cm) - G4ThreeVector(2.3*cm,-1.2*cm,0*cm)).rotateZ(-30*pi/180) ;
    G4VPhysicalVolume* CraniumRight3_phys =  new G4PVPlacement(CraniumRight3_rot,CraniumRight3_shift,"CraniumRight3World",CraniumRight3_Log,Brain2_phys,false,0);

    //Cranium Left 3
    G4EllipticalTube* CraniumLeft3= new G4EllipticalTube("CraniumLeft3",0.5*cm,1*cm,10*cm);
    G4EllipticalTube* CraniumLeft2= new G4EllipticalTube("CraniumLeft2",0.4*cm,2.5*cm,10*cm);
    G4RotationMatrix* CraniumLeft2_rot = new G4RotationMatrix(0,0,(300-280)*pi/180);
    G4UnionSolid* CraniumLeft3_Union1 = new G4UnionSolid("CraniumLeft3+CraniumLeft2",CraniumLeft3,CraniumLeft2,CraniumLeft2_rot,G4ThreeVector(G4ThreeVector(-3.6*cm,0.1*cm,0*cm)
														     -G4ThreeVector(-4.9*cm,-0.2*cm,0*cm)).rotateZ(-280*pi/180));
    G4LogicalVolume* CraniumLeft3_Log = new G4LogicalVolume(CraniumLeft3_Union1,craniumICRUreport46,"CraniumLeft3_Log");
    CraniumLeft3_Log->SetVisAttributes(att2);
    G4RotationMatrix* CraniumLeft3_rot = new G4RotationMatrix(0,0,(280-330)*pi/180);
    G4ThreeVector     CraniumLeft3_shift = G4ThreeVector(G4ThreeVector(-4.9*cm,-0.2*cm,0*cm) - G4ThreeVector(-2.3*cm,-1.2*cm,0*cm)).rotateZ((-330)*pi/180);
    G4VPhysicalVolume* CraniumLeft3_phys = new G4PVPlacement(CraniumLeft3_rot, CraniumLeft3_shift,"CraniumLeft3World",CraniumLeft3_Log,Brain1_phys,false,0);

    //Ear Cavity Right
    G4EllipticalTube* EarCavityRight= new G4EllipticalTube("EarCavityRight",0.4*cm,0.6*cm,10*cm);
    G4LogicalVolume* EarCavityRight_Log = new G4LogicalVolume(EarCavityRight,air,"EarCavityRight_Log");
    EarCavityRight_Log->SetVisAttributes(att2);
    G4RotationMatrix* EarCavityRight_rot = new G4RotationMatrix(0,0,(60-80)*pi/180);
    new G4PVPlacement(EarCavityRight_rot,G4ThreeVector(G4ThreeVector(5.7*cm,-0.7*cm,0*cm) - G4ThreeVector(4.9*cm,-0.2*cm,0*cm)).rotateZ(-80*pi/180) ,"EarCavityRightWorld",EarCavityRight_Log,
		                                                                                                                                                   CraniumRight3_phys,false,0);
    //Ear Cavity Left
    G4EllipticalTube* EarCavityLeft= new G4EllipticalTube("EarCavityLeft",0.4*cm,0.6*cm,10*cm);
    G4LogicalVolume* EarCavityLeft_Log = new G4LogicalVolume(EarCavityLeft,air,"EarCavityLeft_Log");
    EarCavityLeft_Log->SetVisAttributes(att2);
    G4RotationMatrix* EarCavityLeft_rot = new G4RotationMatrix(0,0,(300-280)*pi/180);
    new G4PVPlacement(EarCavityLeft_rot,G4ThreeVector(G4ThreeVector(-5.7*cm,-0.7*cm,0*cm) - G4ThreeVector(-4.9*cm,-0.2*cm,0*cm)).rotateZ(-280*pi/180),"EarCavityLeftWorld",EarCavityLeft_Log,
		                                                                                                                                                   CraniumLeft3_phys,false,0);
    //Eye Bone Right
    G4EllipticalTube* EyeBoneRight= new G4EllipticalTube("EyeBoneRight",0.5*cm,2.5*cm,10*cm);
    G4LogicalVolume* EyeBoneRight_Log = new G4LogicalVolume(EyeBoneRight,craniumICRUreport46,"EyeBoneRight_Log");
    EyeBoneRight_Log->SetVisAttributes(att2);
    G4RotationMatrix* EyeBoneRight_rot = new G4RotationMatrix(0,0, (305-40)*pi/180);
    new G4PVPlacement(EyeBoneRight_rot,G4ThreeVector(G4ThreeVector(2.8*cm,4.6*cm,0*cm) - BrainRight_shift).rotateZ(-40*pi/180),"EyeBoneRightWorld",EyeBoneRight_Log, BrainRight_phys,false,0);

    //Eye Bone Left
    G4EllipticalTube* EyeBoneLeft= new G4EllipticalTube("EyeBoneLeft",0.5*cm,2.5*cm,10*cm);
    G4LogicalVolume* EyeBoneLeft_Log = new G4LogicalVolume(EyeBoneLeft,craniumICRUreport46,"EyeBoneLeft_Log");
    EyeBoneLeft_Log->SetVisAttributes(att2);
    G4RotationMatrix* EyeBoneLeft_rot = new G4RotationMatrix(0,0,(55-320)*pi/180);
    new G4PVPlacement(EyeBoneLeft_rot,G4ThreeVector(G4ThreeVector(-2.8*cm,4.6*cm,0*cm) - BrainLeft_shift).rotateZ(-320*pi/180),"EyeBoneLeftWorld",EyeBoneLeft_Log,BrainLeft_phys,false,0);

    //Nose Fat
    G4EllipticalTube* NoseOutline= new G4EllipticalTube("NoseOutline",1.3*cm,3.5*cm,10*cm);
    G4LogicalVolume* NoseOutline_Log = new G4LogicalVolume(NoseOutline,AdiposemeanZWoodard1986,"NoseOutline_Log");
    NoseOutline_Log->SetVisAttributes(att2);
    G4RotationMatrix* NoseOutline_rot = new G4RotationMatrix(0,0,0*pi/180);
    G4VPhysicalVolume* NoseOutline_phys = new G4PVPlacement(NoseOutline_rot,G4ThreeVector(-0.0*cm,6*cm,0*cm),"NoseOutlineWorld",NoseOutline_Log,Outline_phys,false,0);

    //Nose Bone
    G4EllipticalTube* NoseBone= new G4EllipticalTube("NoseBone",1.1*cm,3.3*cm,10*cm);
    G4LogicalVolume* NoseBone_Log = new G4LogicalVolume(NoseBone,corticalboneICRUreport46,"NoseBone_Log");
    NoseBone_Log->SetVisAttributes(att2);
    G4RotationMatrix* NoseBone_rot = new G4RotationMatrix(0,0,0*pi/180);
    G4VPhysicalVolume* NoseBone_phys = new G4PVPlacement(NoseBone_rot,G4ThreeVector(-0.0*cm,0*cm,0*cm),"NoseBoneWorld",NoseBone_Log,NoseOutline_phys,false,0);

    //Nose Air
    G4EllipticalTube* NoseAir= new G4EllipticalTube("NoseAir",0.9*cm,3*cm,10*cm);
    G4LogicalVolume* NoseAir_Log = new G4LogicalVolume(NoseAir,air,"NoseAir_Log");
    NoseAir_Log->SetVisAttributes(att2);
    G4RotationMatrix* NoseAir_rot = new G4RotationMatrix(0,0,0*pi/180);
    new G4PVPlacement(NoseAir_rot,G4ThreeVector(-0.0*cm,0*cm,0*cm),"NoseAirWorld",NoseAir_Log,NoseBone_phys,false,0);

    G4EllipticalTube* EyeR= new G4EllipticalTube("EyeR",1.5*cm,1.5*cm,10*cm);
    G4LogicalVolume* EyeR_Log = new G4LogicalVolume(EyeR,water,"EyeR_Log");
    EyeR_Log->SetVisAttributes(att2);
    G4RotationMatrix* EyeR_rot = new G4RotationMatrix(0,0,0*pi/180);
    new G4PVPlacement(EyeR_rot,G4ThreeVector(3.0*cm,7*cm,0*cm),"EyeRWorld",EyeR_Log,Outline_phys,false,0);

    G4EllipticalTube* EyeL= new G4EllipticalTube("EyeL",1.5*cm,1.5*cm,10*cm);
    G4LogicalVolume* EyeL_Log = new G4LogicalVolume(EyeL,water,"EyeL_Log");
    EyeL_Log->SetVisAttributes(att2);
    G4RotationMatrix* EyeL_rot = new G4RotationMatrix(0,0,0*pi/180);
    new G4PVPlacement(EyeL_rot,G4ThreeVector(-3.0*cm,7*cm,0*cm),"EyeLWorld",EyeL_Log,Outline_phys,false,0);

  }


  else if(thePhantom =="XCAT"){ // XCAT Phantom
    //vector<G4Material*> theMaterialList;
    // Define a vector of material
    vector< G4String > MaterialName  = {"Air","Water","Skin","Breast","MuscleTissue","Brain","Liver","GallbladderPituitaryGlandTracheaThymusTonsilsUreters"
					,"LungInflated","Oesophagus","Cartilage","Stomach", "Pancreas","Spleen","Ribs2and6","CorticalBone",
					"vertebral_column_C4_ICRU_report_46", "red_marrow_Woodard_1986","Blood",
					"UrinaryBladder","Gland_mean_Z_Woodard_1986","SmallIntestine","SmallIntestine","Gland_mean_Z_Woodard_1986","Heart",
					"Cartilage","Lymph","Uterus","Thyroid","SoftTissue"};

    vector<G4float> Threshold = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29};

    vector<G4float > Density       = {0.001,1.0, 1.09, 0.99, 1.05, 1.045, 1.08, 1.01, 0.26, 1.10,1.10, 1.05, 1.03, 1.06, 1.41,1.84,1.42,
				      1.03,1.01, 1.03,1.02,1.03,1.03,1.02, 1.055, 1.10, 1.0, 1.03, 1.01,1.01};

    for (unsigned int k=0;k<MaterialName.size();k++) theMaterialList.push_back( theMaterial->ConstructMaterial(MaterialName.at(k),Density.at(k)) );

    size_t* materialIDs          = new size_t[NbinsX*NbinsY*NbinsZ];
    unsigned int p               = 0;

    for (G4int k=0;k<NbinsZ;k++) {
    for (G4int j=0;j<NbinsY;j++) {
    for (G4int i=0;i<NbinsX;i++) {
      G4int huUnit = hu->GetBinContent(i+1,j+1,k+1);
      materialIDs[p] = huUnit;
      p++;
    }
    }
    }
    // Voxel volume
    G4Box* phantomVoxel_vol              = new G4Box("phantomVoxel_vol",VoxelHalfX,VoxelHalfY,VoxelHalfZ);
    G4LogicalVolume *phantomVoxel_log    = new G4LogicalVolume(phantomVoxel_vol,air,"phantomVoxel_log",0,0,0);

    param = new G4PhantomParameterisation();
    param->SetVoxelDimensions(VoxelHalfX,VoxelHalfY,VoxelHalfZ);
    param->SetNoVoxel(NbinsX,NbinsY,NbinsZ);
    param->SetMaterials(theMaterialList);
    param->SetMaterialIndices(materialIDs);
    param->SetSkipEqualMaterials(true);
    param->BuildContainerSolid(cont_phys);
    param->CheckVoxelsFillContainer(cont_vol->GetXHalfLength(),cont_vol->GetYHalfLength(),cont_vol->GetZHalfLength());
    phantomPhys = new G4PVParameterised("phantomPhys",phantomVoxel_log,cont_log,kZAxis,param->GetNoVoxel(),param);
    phantomPhys->SetRegularStructureId(1);
  }
  else{
    cout<<"The phantom name cannot be found."<<endl;
  }

  return physWorld;
}
