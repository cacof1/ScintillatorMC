#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "Analysis.hh"
#include "OrganicMaterial.hh"
#include "G4NistManager.hh"
#include "DetectorConstruction.hh"
#include "SteppingAction.hh"
#include "G4EmCalculator.hh"
#include "G4IonTable.hh"
#include "G4RunManager.hh"
#include <TROOT.h>
#include "G4ExceptionHandler.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4EmCalculator.hh"
#include "G4Proton.hh"
#include "G4UnitsTable.hh"
#include "pCTconfig.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UImanager.hh"
#include <stdlib.h>
#ifdef VIS
#include "G4VisExecutive.hh"

#endif
#include "G4UIExecutive.hh"
#include <iostream>
#include <fstream>
using namespace std;
void calcRSP(PrimaryGeneratorAction *);
void calcStoppingPower(PrimaryGeneratorAction *, DetectorConstruction*);
int main(int ,char** argv) {
  gROOT->ProcessLine("#include <vector>");
  const string configFile = G4String(argv[1]);
  pCTconfig cfg(configFile.data()); // Create a class instance for parsing the configuration file             
  cout<<"Config file "<<configFile<<endl; 

  G4int nProtons = 1000;
  cfg.addItem("nProtons", nProtons);

  G4float Energy = 200; // MeV
  cfg.addItem("Energy", Energy);   

  G4String Model  = "XCAT";       // []
  cfg.addItem("Model", Model);     

  G4String CTPath     = "None";
  cfg.addItem("CTPath", CTPath);   

  G4float Angle  = 0; // Degrees
  cfg.addItem("Angle", Angle);   

  G4float Thickness  = 30; // cm
  cfg.addItem("Thickness", Thickness);   

  G4int   thread  = atoi(argv[2]); // []
  cfg.addItem("thread", thread);   

  G4int  ANumber  = 1; // Atomic Number
  cfg.addItem("ANumber", ANumber);
  
  G4int NPB       = 100; // []
  cfg.addItem("NPB", NPB);
  
  G4float sigmaX_pos = 5;// mm
  cfg.addItem("sigmaX_pos", sigmaX_pos);

  G4float sigmaY_pos = 5; // mm
  cfg.addItem("sigmaY_pos", sigmaY_pos);
    
  G4float sigma_AngX = 25; // mRad
  cfg.addItem("sigma_AngX", sigma_AngX);
  
  G4float sigma_AngY = 25;// mRad
  cfg.addItem("sigma_AngY", sigma_AngY);

  G4float fieldSizeY = 300;// mm
  cfg.addItem("fieldSizeY", fieldSizeY);

  G4float centerY = 0;// mm
  cfg.addItem("centerY", centerY);  

  G4float fieldSizeZ = 300;// mm
  cfg.addItem("fieldSizeZ", fieldSizeZ);      

  G4float centerZ = 0;// mm
  cfg.addItem("centerZ", centerZ);      

  G4int saveYXProj = 0; // []
  cfg.addItem("saveYXProj", saveYXProj);

  G4int saveYZProj = 0; // []
  cfg.addItem("saveYZProj", saveYZProj);

  G4int saveZXProj = 0; // []
  cfg.addItem("saveZXProj", saveZXProj);

  G4int saveTTree = 0; // []
  cfg.addItem("saveTTree", saveTTree);    

  G4int particleCount = 0; // []
  cfg.addItem("particleCount", particleCount);  

  G4String SourceType  = "Beam";       // []
  cfg.addItem("SourceType", SourceType);     

  
  cfg.Configure();
  CLHEP::RanecuEngine *theRanGenerator = new CLHEP::RanecuEngine;  
  theRanGenerator->setSeed(thread);
  CLHEP::HepRandom::setTheEngine(theRanGenerator);
  G4RunManager* runManager   = new G4RunManager;  
  runManager->SetUserInitialization(new PhysicsList());
  DetectorConstruction* myDC = new DetectorConstruction();
  PrimaryGeneratorAction *theGenerator =  new PrimaryGeneratorAction();
  Analysis* theAnalysis      = new Analysis();


  runManager->SetUserAction(theGenerator);
  //runManager->SetUserAction( new SteppingAction());
  runManager->SetUserInitialization(myDC);
  runManager->SetVerboseLevel(0);
  runManager->Initialize();

  //G4UImanager * UImanager = G4UImanager::GetUIpointer();
  //UImanager->ApplyCommand("/run/setCut 1 mm");
  //UImanager->ApplyCommand("/process/inactivate msc all");
  //UImanager->ApplyCommand("/process/inactivate had all");
  //UImanager->ApplyCommand("/process/eLoss/fluct false");
  //UImanager->ApplyCommand("/process/eLoss/CSDARange true");
  //UImanager->ApplyCommand("/run/verbose 1");
  //UImanager->ApplyCommand("/event/verbose 1");
  //UImanager->ApplyCommand("/tracking/verbose 2");
  //G4UIExecutive * ui = new G4UIExecutive(argc,argv);  

  #ifdef VIS
  G4UImanager * UImanager = G4UImanager::GetUIpointer();
  G4VisManager* visManager = new G4VisExecutive;
  visManager->SetVerboseLevel(0);
  visManager->Initialize();
  //G4UImanager * UImanager = G4UImanager::GetUIpointer();
  G4cout << " UI session starts ..." << G4endl;
  G4UIExecutive* ui = new G4UIExecutive(2, argv);
  UImanager->ApplyCommand("/control/execute vis.mac");
  ui->SessionStart();
  #endif

  int NProton_tot = cfg.item_int["NPB"]*cfg.item_int["NPB"]*cfg.item_int["nProtons"];
  runManager->BeamOn(NProton_tot);
  theAnalysis->SaveAndClose();
  calcRSP(theGenerator);
  //calcStoppingPower(theGenerator, myDC);
  //delete visManager;
  return 0;
  delete runManager;
  }

void calcRSP(PrimaryGeneratorAction* theGenerator){
  G4ParticleDefinition* particle = theGenerator->particle;
  G4EmCalculator* emCal = new G4EmCalculator;
  ofstream myfile;
  myfile.open ("RSP.txt");
  myfile<<"Density RelativeElectronDensity IValue RSP Name"<<endl; 
  G4MaterialTable *theMaterialTable = G4Material::GetMaterialTable();
  G4float waterelectrondensity = 0;
  for(size_t i =0;i<theMaterialTable->size();i++){
      
    G4int I = 0;
    G4float tot =0;
    G4Material* water = theMaterialTable->at(0);
    waterelectrondensity = water->GetElectronDensity()/(g/cm3);
    for(int j=1;j<50000;j++){
      G4float dedx_w = emCal->ComputeElectronicDEDX( float(j)/10*MeV,particle,water);
      G4float dedx_b = emCal->ComputeElectronicDEDX( float(j)/10*MeV,particle,theMaterialTable->at(i));
      tot +=dedx_b/dedx_w;
      I+=1;
    }
    G4float RSP = tot/I;
    G4float electrondensity = theMaterialTable->at(i)->GetElectronDensity()/(g/cm3);
    myfile<<theMaterialTable->at(i)->GetDensity()/(g/cm3)<<" "<<electrondensity/waterelectrondensity<<" "<<theMaterialTable->at(i)->GetIonisation()->GetMeanExcitationEnergy()/eV<<" "<<RSP<<" "<<theMaterialTable->at(i)->GetName()<<endl;
    }
  myfile.close();
  
}

void calcStoppingPower(PrimaryGeneratorAction* theGenerator, DetectorConstruction* myDC){
  G4ParticleDefinition* particle = theGenerator->particle;
  G4EmCalculator* emCal = new G4EmCalculator;

  ofstream myfile;
  myfile.open ("Water_Geant4.dat");
  //G4MaterialTable *theMaterialTable = G4Material::GetMaterialTable();
  G4Material* water = myDC->water;//theMaterialTable->at(0);
  
  cout<<emCal->GetCSDARange(105.43*MeV,particle, myDC->water)*mm<<endl;
  for(int j=1;j<50000;j++){
    G4float dedx_w = emCal->ComputeElectronicDEDX( float(j)/10*MeV,particle,water);
    myfile<<float(j)/10*MeV<<" "<<dedx_w*MeV/mm<<" "<<endl;
  }
  myfile.close();

}


