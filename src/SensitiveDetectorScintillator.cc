#include "SensitiveDetectorScintillator.hh"
#include "Analysis.hh"

SensitiveDetectorScintillator::SensitiveDetectorScintillator(G4String name):G4VSensitiveDetector(name),theName(name)
{
  theAnalysis  = Analysis::GetInstance();  
}

G4bool SensitiveDetectorScintillator::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  theAnalysis->FillScintillatorDose(aStep);
  return true;
}
