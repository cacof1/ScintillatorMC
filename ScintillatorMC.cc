#include "PhysicsList.hh"
//#include "EmPhysics_pCT.hh"
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

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UImanager.hh"

#ifdef VIS
#include "G4VisExecutive.hh"

#endif
#include "G4UIExecutive.hh"
#include <iostream>
#include <fstream>
using namespace std;
void calcRSP(PrimaryGeneratorAction *);
void calcStoppingPower(PrimaryGeneratorAction *, DetectorConstruction*);
int main(int argc,char** argv) {
  gROOT->ProcessLine("#include <vector>");


  if(argc<=7) {
    cout<<"Please input the following arguments : NParticules Energy Model Angle Thick Thread ANumber"<<endl;
    return 0;
  }
  G4int nProtons  = atoi(argv[1]); 
  G4double Energy = atof(argv[2]); // MeV
  G4String Model  = argv[3];
  G4double angle  = atof(argv[4]); // Degrees
  G4double thick  = atof(argv[5]); // cm
  G4int   thread  = atoi(argv[6]); 
  G4int  ANumber  = atoi(argv[7]); //Atomic Number
  G4String CT     = "";
  if(argc==9) CT     = argv[8];
  CLHEP::RanecuEngine *theRanGenerator = new CLHEP::RanecuEngine;  
  theRanGenerator->setSeed(thread);
  CLHEP::HepRandom::setTheEngine(theRanGenerator);
  G4String paraWorldName = "";
  G4RunManager* runManager   = new G4RunManager;  
  runManager->SetUserInitialization(new PhysicsList(paraWorldName));
  DetectorConstruction* myDC = new DetectorConstruction(Model,angle,thick,CT);
  cout<<Energy<<" "<<ANumber<<" "<<nProtons<<endl;
  PrimaryGeneratorAction *theGenerator =  new PrimaryGeneratorAction(Energy,ANumber, nProtons);
  cout<<__LINE__<<endl;
  Analysis* theAnalysis      = new Analysis(thread,angle,Model);
  cout<<__LINE__<<endl;
  //G4ExceptionHandler* handler = new G4ExceptionHandler();
  runManager->SetUserAction(theGenerator);
  runManager->SetUserAction( new SteppingAction() );
  runManager->SetUserInitialization( myDC );
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
  G4UIExecutive* ui = new G4UIExecutive(argc, argv);
  UImanager->ApplyCommand("/control/execute vis.mac");
  ui->SessionStart();
  #endif

  runManager->BeamOn( nProtons );  
  theAnalysis->Save();
  //calcRSP(theGenerator);
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
  G4double waterelectrondensity = 0;
  for(size_t i =0;i<theMaterialTable->size();i++){
      
    G4int I = 0;
    G4double tot =0;
    G4Material* water = theMaterialTable->at(0);
    waterelectrondensity = water->GetElectronDensity()/(g/cm3);
    for(int j=1;j<50000;j++){
      G4double dedx_w = emCal->ComputeElectronicDEDX( double(j)/10*MeV,particle,water);
      G4double dedx_b = emCal->ComputeElectronicDEDX( double(j)/10*MeV,particle,theMaterialTable->at(i));
      tot +=dedx_b/dedx_w;
      I+=1;
    }
    G4double RSP = tot/I;
    G4double electrondensity = theMaterialTable->at(i)->GetElectronDensity()/(g/cm3);
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
    G4double dedx_w = emCal->ComputeElectronicDEDX( double(j)/10*MeV,particle,water);
    myfile<<double(j)/10*MeV<<" "<<dedx_w*MeV/mm<<" "<<endl;
  }
  myfile.close();

}


