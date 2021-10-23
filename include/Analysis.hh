#ifndef Analysis_hh
#define Analysis_hh
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "TTree.h"
#include "TMap.h"
#include "TH1F.h"
#include <vector>
#include <iostream>
#include "PrimaryGeneratorAction.hh"
#include "pCTconfig.hh"
using namespace std;
class G4Step;
class TFile ;
class TTree ;
class TH3F  ;
class TH2F  ;
class TH1F  ;
class TProfile2D;
class TMap  ;
class PrimaryGeneratorAction;
class DetectorConstruction;
class SteppingAction;
class Analysis
{
public:

  Analysis();//G4int,G4double, G4String);
  ~Analysis();
  static inline Analysis* GetInstance() { return theAnalysis; }
  void analyseHit(G4Step*,G4String);
  TTree  *t; // Phasespace
  TTree  *t2; // Pencil beam dose and LET
  TFile *f1;

  pCTconfig* theConfig = pCTconfig::GetInstance();
  TH3F* Edep_Tot;
  TH3F* L_Tot;  
  TH3F* Entries_Tot;
  TH3F* LET_Tot;

  TH2F* XYProj_Tot;
  TH2F* XZProj_Tot;
  TH2F* YZProj_Tot;

  TH2F* XYProj_Tot_Q;
  TH2F* XZProj_Tot_Q;
  TH2F* YZProj_Tot_Q;
  void SaveAndClose();
  void SaveAndReset();  
  void RearFrontDetector(G4Step* aStep, G4String theName);
  void FillScintillatorDose(G4Step*  aStep);
  double findWET(double, double);

  std::string  proc_name, part_name;
  G4double TotEnergyDeposit = 0.;

  int NbinsX, NbinsY, NbinsZ, binGlobal;
  int idPB;
  float Xmin, Xmax, Ymin, Ymax, Zmin, Zmax;
private:
  pair<map<int,pair<float,pair<float,int>>>::iterator,bool> ret;
  map<int,pair<float,pair<float,int>>>::iterator it;
  static Analysis* theAnalysis;
  PrimaryGeneratorAction *theGenerator;
  DetectorConstruction   *theDetector ;
  SteppingAction         *theSteppingAction;
  float x_scint, y_scint ,z_scint;
  float Estop_scint;
  float L,dX; // quenched light emission
  float kB,A; // quenching parameters
  float  x0=0,y0=0,z0=0,px0=1,py0=0,pz0=0;

  float  x1,y1,z1,px1,py1,pz1;
  float  theta_y1, theta_z1;
  float  Einit,Estop, EnergyMid,LET;
  int N;
  vector<double> Energy;
  vector<double> dEdXBins;

  G4int  idPBY, idPBZ;

  // Non Quenched
  /*
  TH1D* PDD[theGenerator->NPBY*theGenerator->NPBZ];  
  TH2F* YXProj[theGenerator->NPBY*theGenerator->NPBZ];
  TH2F* ZXProj[theGenerator->NPBY*theGenerator->NPBZ];
  TH2F* YZProj[theGenerator->NPBY*theGenerator->NPBZ];
  */
 
  //Quenched
  /*TH1D* PDD_Q[theGenerator->NPBY*theGenerator->NPBZ];  
  TH2F* YXProj_Q[theGenerator->NPBY*theGenerator->NPBZ];
  TH2F* ZXProj_Q[theGenerator->NPBY*theGenerator->NPBZ];
  TH2F* YZProj_Q[theGenerator->NPBY*theGenerator->NPBZ];*/

  //TH1F* PDD_Q;
  TH2F* YXProj_Q;
  TH2F* ZXProj_Q;
  TH2F* YZProj_Q;

  //Single Event radiograph
  TProfile2D* Front;
  TProfile2D* Back;

  //std::vector<  map<int,pair<float,float>> >  Edep_PB = vector< map<int,pair<float,float>> >(); //theGenerator->NPBY*theGenerator->NPBZ];
  //map<int,pair<float,pair<float,int>>> Edep_PB[theGenerator->NPBY*theGenerator->NPBZ];
};
#endif
