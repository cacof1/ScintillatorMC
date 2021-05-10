#ifndef Analysis_hh
#define Analysis_hh
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "TTree.h"
#include "TMap.h"
#include "PrimaryGeneratorAction.hh"

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

  Analysis(G4int,G4double, G4String);
  ~Analysis();
  static inline Analysis* GetInstance() { return theAnalysis; }
  void analyseHit(G4Step*,G4String);
  TTree  *t; // Phasespace
  TTree  *t2; // Pencil beam dose and LET
  TFile *f1;

  //map<int,pair<float,float>> Edep_PB[NPBY*NPBZ];
  map<int,pair<float,pair<float,int>>> Edep_PB[NPBY*NPBZ];
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
  void Save();
  void RearFrontDetector(G4Step* aStep, G4String theName);
  void FillScintillatorDose(G4Step*  aStep);
  double findWET(double, double);
  vector<double>* tracks_X;
  vector<double>* tracks_Y;
  vector<double>* tracks_Z;
  vector<double>* Eloss;
  vector<double>* Length;
  vector<double>* tracks_E;
  vector<double>* Radlen;
  vector<TString>*  mat_name;
  vector<TString>*  vol_name;
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
  TH1D* PDD[NPBY*NPBZ];  
  TH2F* YXProj[NPBY*NPBZ];
  TH2F* ZXProj[NPBY*NPBZ];
  TH2F* YZProj[NPBY*NPBZ];
  */

  //Quenched
  TH1F* PDD_Q[NPBY*NPBZ];  
  TH2F* YXProj_Q[NPBY*NPBZ];  
  TH2F* ZXProj_Q[NPBY*NPBZ];
  TH2F* YZProj_Q[NPBY*NPBZ];
  
  //Single Event radiograph
  TProfile2D* Front;
  TProfile2D* Back;

};
#endif
