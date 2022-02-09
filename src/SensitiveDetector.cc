#include "SensitiveDetector.hh"
#include "Analysis.hh"

SensitiveDetector::SensitiveDetector(G4String name):G4VSensitiveDetector(name),theName(name)
{
  theAnalysis = Analysis::GetInstance();  
}

G4bool SensitiveDetector::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  if ( aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary){
    if(theName == "FrontTracker" || theName == "RearTracker") theAnalysis->RearFrontDetector(aStep, theName);
  }
  //if( theName == "Scintillator" && aStep->GetPreStepPoint()->GetStepStatus() != fGeomBoundary)  theAnalysis->FillScintillatorDose(aStep);
  if( theName == "Scintillator")  theAnalysis->FillScintillatorDose(aStep);
  return true;
}
