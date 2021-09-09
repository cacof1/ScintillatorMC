#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "G4EmConfigurator.hh"
#include "globals.hh"

class G4VPhysicsConstructor;
class HadrontherapyStepMax;
class HadrontherapyPhysicsListMessenger;

class PhysicsList: public G4VModularPhysicsList
{
public:
  PhysicsList();
  virtual ~PhysicsList();

  void ConstructParticle();

  void SetCuts();
  void SetCutForGamma(G4double);
  void SetCutForElectron(G4double);
  void SetCutForPositron(G4double);
  
  void AddPhysicsList(const G4String& name);
  void ConstructProcess();

  void AddStepMax();
  HadrontherapyStepMax* GetStepMaxProcess() {return stepMaxProcess;};
  void AddPackage(const G4String& name);

private:
  void AddParallelWorldProcess();
  G4EmConfigurator em_config;

  G4double cutForGamma;
  G4double cutForElectron;
  G4double cutForPositron;

  G4bool helIsRegisted;
  G4bool bicIsRegisted;
  G4bool biciIsRegisted;
  G4bool locIonIonInelasticIsRegistered;
  G4bool radioactiveDecayIsRegisted;

  G4String                             emName;
  G4VPhysicsConstructor*               emPhysicsList;
  G4VPhysicsConstructor*               decPhysicsList;
  std::vector<G4VPhysicsConstructor*>  hadronPhys;

  HadrontherapyStepMax* stepMaxProcess;

};

#endif
