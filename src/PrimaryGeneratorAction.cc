#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "G4EmCalculator.hh"
#include "G4NistManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4UImanager.hh"
#include "G4Proton.hh"
#include "G4Alpha.hh"
#include "G4Deuteron.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include "G4Proton.hh"
#include "G4IonTable.hh"
#include  <fstream>
#include  <sstream>
#include <math.h>
#include "G4NavigationHistory.hh"
#include "G4TouchableHistory.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TransportationManager.hh"
#include "TFile.h"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Analysis.hh"

using namespace std;

PrimaryGeneratorAction* PrimaryGeneratorAction::theGenerator = NULL;

PrimaryGeneratorAction::PrimaryGeneratorAction()
{ 
  theConfig = pCTconfig::GetInstance();
  NPBY = theConfig->item_int["NPB"];
  NPBZ = theConfig->item_int["NPB"];

  A        = theConfig->item_int["ANumber"];
  Nprotons = theConfig->item_int["nProtons"];
  ENER     = theConfig->item_float["Energy"];
  
  theGenerator    = this;
  theDetector     = DetectorConstruction::GetInstance();
  particleSource  = new G4GeneralParticleSource();

  //Generic beam (no phase space)
  if(A==1) particle = G4Proton::Proton();//G4IonTable::GetIonTable()->GetIon(1,1,0); // proton
  else if(A==2) particle = G4Deuteron::Deuteron();
  else if(A==4) particle = G4Alpha::Alpha();
  else particle = G4IonTable::GetIonTable()->GetIon(int(A/2),A,0); // rest
  particleSource->SetParticleDefinition(particle);
  

  // Mono-energetic
  eneDist = particleSource->GetCurrentSource()->GetEneDist();
  eneDist->SetEnergyDisType("Mono");
  eneDist->SetMonoEnergy(ENER*A*MeV);

  //Energy Spectrum
  //eneDist = particleSource->GetCurrentSource()->GetEneDist();
  //eneDist->SetEnergyDisType("Arb");
  //eneDist->ArbEnergyHistoFile("macro/source/protonKEspectrumAtZ-0.1cm_71_163.9MeV.dat");
  //eneDist->ArbInterpolate("Lin");
  
  ////////////////////////////////////////////////////////////////////////////////////
  //Source 1
  //Position

  particleSource->SetCurrentSourceIntensity(0.75);
  posDist = particleSource->GetCurrentSource()->GetPosDist();
  posDist->SetPosDisType(theConfig->item_str["SourceType"]);
  posDist->SetPosDisShape("Plane");
  posDist->SetBeamSigmaInX(theConfig->item_float["sigmaX_pos"]*mm);//6.393*mm;
  posDist->SetBeamSigmaInY(theConfig->item_float["sigmaY_pos"]*mm);//6.5524*mm;

  //Angular Specturm
  angDist = particleSource->GetCurrentSource()->GetAngDist();
  angDist->SetAngDistType("beam2d");
  angDist->SetBeamSigmaInAngX(theConfig->item_float["sigma_AngX"]/1000.);
  angDist->SetBeamSigmaInAngY(theConfig->item_float["sigma_AngY"]/1000.);    // Rad -> mRad


  particleSource->SetParticleDefinition(particle);
  ////////////////////////////////////////////////////////////////////////////////////
  //Source 2
  //Position

  particleSource->AddaSource(0.25);
  posDist = particleSource->GetCurrentSource()->GetPosDist();
  posDist->SetPosDisType(theConfig->item_str["SourceType"]);  
  posDist->SetPosDisShape("Plane");
  posDist->SetBeamSigmaInX(theConfig->item_float["sigmaX_pos"]*mm);//6.393*mm;
  posDist->SetBeamSigmaInY(theConfig->item_float["sigmaY_pos"]*mm);//6.5524*mm;
  
  //Angular Specturm
  angDist = particleSource->GetCurrentSource()->GetAngDist();
  angDist->SetAngDistType("planar");
  angDist->SetBeamSigmaInAngX(theConfig->item_float["sigma_AngX"]/1000.);
  angDist->SetBeamSigmaInAngY(theConfig->item_float["sigma_AngY"]/1000.);    // Rad -> mRad

  particleSource->SetParticleDefinition(particle); 
  //Pencil beam parameters
  nProtonsGenerated     = 0;
  nProtonsPerPencilBeam = 0;
  idPBGlobal            = 0;// Start at pencil beam 0
			   
  fieldSizeY = theConfig->item_float["fieldSizeY"]; //2*(theDetector->ScintHalfY);//300 mm
  fieldSizeZ = theConfig->item_float["fieldSizeZ"]; //2*(theDetector->ScintHalfZ);//300 mm
  
  if(theConfig->item_int["NPB"]==1){ // For the calibration
    PencilBeamPosY.push_back(0);
    PencilBeamPosZ.push_back(0);
  }
  else{
    PencilBeamPosY = linspace(-fieldSizeY/2 + theConfig->item_float["centerY"], fieldSizeY/2 + theConfig->item_float["centerY"], NPBY); // mm
    PencilBeamPosZ = linspace(-fieldSizeZ/2 + theConfig->item_float["centerZ"], fieldSizeZ/2 + theConfig->item_float["centerZ"], NPBZ); // mm
  }

  //Initial Setup
  x0           = -1*theDetector->PhantomHalfX - theDetector->midX -1*mm -10*mm;
  y0           = PencilBeamPosY[0]; 
  z0           = PencilBeamPosZ[0];
  particleSource->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(x0,y0,z0));
  angleY = atan2(-y0,theConfig->item_float["isoMagnetY"]*cm);
  angleZ = atan2(-z0,theConfig->item_float["isoMagnetZ"]*cm);
  
  pz = cos(angleY) * cos(angleZ);
  py = sin(angleY) * cos(angleZ);
  px = sin(angleZ);

  Momentum = G4ThreeVector(px, py, pz).unit();    
  particleSource->SetCurrentSourceto(0);
  particleSource->GetCurrentSource()->GetAngDist()->DefineAngRefAxes("angref1", Momentum);
  particleSource->GetCurrentSource()->GetPosDist()->SetPosRot1(Momentum);
}
  
PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  theGenerator = NULL;
  delete particleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  theAnalysis  = Analysis::GetInstance(); // singleton problems, might need to de-singleton it
  nProtonsGenerated++; nProtonsPerPencilBeam++;

  if(nProtonsGenerated%20000==0) cout << nProtonsGenerated << "/"<<Nprotons*NPBY*NPBZ<< endl;
  if(nProtonsPerPencilBeam >= theConfig->item_int["nProtons"] ){

    idPBGlobal +=1;
    nProtonsPerPencilBeam = 0;
    idPBY        = (idPBGlobal - (idPBGlobal%NPBY))/NPBY;
    idPBZ        = idPBGlobal%NPBY;  

    x0           = -1*theDetector->PhantomHalfX - theDetector->midX -1*mm -10*mm;    
    y0           = PencilBeamPosY[idPBY]; 
    z0           = PencilBeamPosZ[idPBZ];
    angleY = atan2(-y0,theConfig->item_float["isoMagnetY"]*cm);
    angleZ = atan2(-z0,theConfig->item_float["isoMagnetZ"]*cm);
    
    pz = cos(angleY) * cos(angleZ);
    py = sin(angleY) * cos(angleZ);
    px = sin(angleZ);
    Momentum = G4ThreeVector(px, py, pz).unit();

    
    particleSource->SetCurrentSourceto(0);
    particleSource->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(x0,y0,z0));
    particleSource->GetCurrentSource()->GetAngDist()->DefineAngRefAxes("angref1", Momentum);
    particleSource->GetCurrentSource()->GetPosDist()->SetPosRot1(Momentum);

    particleSource->SetCurrentSourceto(1);
    particleSource->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(x0,y0,z0));
    particleSource->GetCurrentSource()->GetAngDist()->DefineAngRefAxes("angref1", Momentum);
    particleSource->GetCurrentSource()->GetPosDist()->SetPosRot1(Momentum);
    Einit = particleSource->GetCurrentSource()->GetParticleEnergy();


    theAnalysis->SaveAndReset();      
  }
  particleSource->GeneratePrimaryVertex(anEvent);
  
}

vector<G4double> PrimaryGeneratorAction::linspace(double start, double stop, int NStep) {
  std::vector<double> array;

  if (NStep == 0) return array; 
  if (NStep == 1)
    {
      array.push_back(start);
      return array;
    }
  double delta = (stop - start) / (NStep);
  for(int i=0; i < NStep; ++i)  array.push_back(start + delta * i);
  array.push_back(stop);
  return array;
}



