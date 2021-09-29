#include "SensitiveDetectorTracker.hh"
#include "Analysis.hh"

SensitiveDetectorTracker::SensitiveDetectorTracker(G4String name):G4VSensitiveDetector(name),theName(name)
{
  theAnalysis = Analysis::GetInstance();  
}

G4bool SensitiveDetectorTracker::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  if ( aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary){
    theAnalysis->RearFrontDetector(aStep, theName);
  }
  return true;
}
