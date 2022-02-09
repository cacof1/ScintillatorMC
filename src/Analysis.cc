#include "vector"
#include "G4ios.hh"
#include "TTree.h"
#include "TFile.h"
#include "TH3D.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TProfile2D.h"
#include "TMap.h"
#include "globals.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4Event.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "SteppingAction.hh"
#include "G4VProcess.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Analysis.hh"

Analysis* Analysis::theAnalysis=NULL;
Analysis::~Analysis(){theAnalysis=NULL;}

Analysis::Analysis()
{
  theAnalysis       = this;
  theGenerator      = PrimaryGeneratorAction::GetInstance();
  theDetector       = DetectorConstruction::GetInstance();
  theConfig         = pCTconfig::GetInstance();


  f1 = new TFile(Form("%s_%.0f_%.1f_%d_%d.root",theConfig->item_str["Model"].data(),theConfig->item_float["Energy"],theConfig->item_float["angle"],theConfig->item_int["thread"],theConfig->item_int["ANumber"]),"recreate");
  f1->mkdir("PDDList");  f1->mkdir("YZProj");  f1->mkdir("YXProj");  f1->mkdir("ZXProj"); // Normal
  f1->mkdir("PDDList_Q");  f1->mkdir("ZXProj_Q");  f1->mkdir("YZProj_Q");  f1->mkdir("YXProj_Q"); // Quenched

  NbinsX = 300; NbinsY = 300; NbinsZ = 300;

  //Full histogram
  Edep_Tot    = new TH3F("Edep_Tot", "Edep_Tot", NbinsX, -theDetector->ScintHalfX, theDetector->ScintHalfX,
		                                 NbinsY, -theDetector->ScintHalfY, theDetector->ScintHalfY,
		                                 NbinsZ, -theDetector->ScintHalfZ, theDetector->ScintHalfZ);


  L_Tot       = new TH3F("L_Tot", "L_Tot", NbinsX, -theDetector->ScintHalfX, theDetector->ScintHalfX,
		                                 NbinsY, -theDetector->ScintHalfY, theDetector->ScintHalfY,
		                                 NbinsZ, -theDetector->ScintHalfZ, theDetector->ScintHalfZ);


  LET_Tot    = new TH3F("LET_Tot", "LET_Tot",    NbinsX, -theDetector->ScintHalfX, theDetector->ScintHalfX,
		                                 NbinsY, -theDetector->ScintHalfY, theDetector->ScintHalfY,
	 	                                 NbinsZ, -theDetector->ScintHalfZ, theDetector->ScintHalfZ);

  Entries_Tot = new TH3F("Entries_Tot", "Entries_Tot", NbinsX, -theDetector->ScintHalfX, theDetector->ScintHalfX,
 		                                       NbinsY, -theDetector->ScintHalfY, theDetector->ScintHalfY,
                                                       NbinsZ, -theDetector->ScintHalfZ, theDetector->ScintHalfZ);
  // Full projection cumulative
  /*XYProj_Tot  = new TH2F("XYProj_Tot", "XYProj_Tot", NbinsX, -theDetector->ScintHalfX, theDetector->ScintHalfX,
			       NbinsY, -theDetector->ScintHalfY, theDetector->ScintHalfY);

  XZProj_Tot  = new TH2F("XZProj_Tot", "XZProj_Tot", NbinsX, -theDetector->ScintHalfX, theDetector->ScintHalfX,
			       NbinsZ, -theDetector->ScintHalfZ, theDetector->ScintHalfZ);

  YZProj_Tot  = new TH2F("YZProj_Tot", "YZProj_Tot", NbinsY, -theDetector->ScintHalfY, theDetector->ScintHalfY,
			       NbinsZ, -theDetector->ScintHalfZ, theDetector->ScintHalfZ);

  XYProj_Tot_Q  = new TH2F("XYProj_Tot_Q", "XYProj_Tot_Q", NbinsX, -theDetector->ScintHalfX, theDetector->ScintHalfX,
				 NbinsY, -theDetector->ScintHalfY, theDetector->ScintHalfY);

  XZProj_Tot_Q  = new TH2F("XZProj_Tot_Q", "XZProj_Tot_Q", NbinsX, -theDetector->ScintHalfX, theDetector->ScintHalfX,
				 NbinsZ, -theDetector->ScintHalfZ, theDetector->ScintHalfZ);

  YZProj_Tot_Q  = new TH2F("YZProj_Tot_Q", "YZProj_Tot_Q", NbinsY, -theDetector->ScintHalfY, theDetector->ScintHalfY,
				 NbinsZ, -theDetector->ScintHalfZ, theDetector->ScintHalfZ);
  */
  //Single event radiograph
  Front         = new TProfile2D("Front", "Front", NbinsY, -theDetector->ScintHalfY, theDetector->ScintHalfY, NbinsZ, -theDetector->ScintHalfZ, theDetector->ScintHalfZ, "s");
  Back          = new TProfile2D("Back", "Back", NbinsY, -theDetector->ScintHalfY, theDetector->ScintHalfY, NbinsZ, -theDetector->ScintHalfZ, theDetector->ScintHalfZ, "s");

  std::string line;
  std::ifstream SPWater ("Water_Geant4_P.dat");
  double data[3];
  while(getline(SPWater, line)) {
    stringstream ss(line);
    for(int i=0;i<3;i++) ss >> data[i];
    Energy.push_back(data[0]);
    dEdXBins.push_back(data[1]);
  }

  if(theConfig->item_int["saveTTree"] ){ //Phasespace
    t = new TTree("phase","PS");
    t->Branch("x0",&x0,"x0/F");
    t->Branch("y0",&y0,"y0/F");
    t->Branch("z0",&z0,"z0/F");

    t->Branch("px0",&px0,"px0/F");
    t->Branch("py0",&py0,"py0/F");
    t->Branch("pz0",&pz0,"pz0/F");
    t->Branch("Einit",&Einit,"Einit/F");

    t->Branch("x1",&x1,"x1/F");
    t->Branch("y1",&y1,"y1/F");
    t->Branch("z1",&z1,"z1/F");

    t->Branch("px1",&px1,"px1/F");
    t->Branch("py1",&py1,"py1/F");
    t->Branch("pz1",&pz1,"pz1/F");
    t->Branch("Estop",&Estop,"Estop/F");

    t->Branch("proc_name",&proc_name);
    t->Branch("part_name",&part_name);
    t->Branch("idPBY",&idPBY,"idPBY/I");
    t->Branch("idPBZ",&idPBZ,"idPBZ/I");
  }
  /*PDD    = new TH1F(Form("PDD_Q_%d",0),Form("PDD_Q_%d",0),NbinsX, -theDetector->ScintHalfX, theDetector->ScintHalfX);
  ZXProj = new TH2F(Form("ZXProj_%d",0), Form("ZXProj_%d",0), NbinsZ, -theDetector->ScintHalfZ, theDetector->ScintHalfZ,NbinsX, -theDetector->ScintHalfX, theDetector->ScintHalfX);
  YXProj = new TH2F(Form("YXProj_%d",0), Form("YXProj_%d",0), NbinsY, -theDetector->ScintHalfY, theDetector->ScintHalfY,NbinsX, -theDetector->ScintHalfX, theDetector->ScintHalfX);
  YZProj = new TH2F(Form("YZProj_%d",0), Form("YZProj_%d",0), NbinsY, -theDetector->ScintHalfY, theDetector->ScintHalfY,NbinsZ, -theDetector->ScintHalfZ, theDetector->ScintHalfZ);  */

  //PDD_Q    = new TH1F(Form("PDD_Q_%d",0),Form("PDD_Q_%d",0),NbinsX, -theDetector->ScintHalfX, theDetector->ScintHalfX);
  ZXProj_Q = new TH2F(Form("ZXProj_Q_%d",0), Form("ZXProj_Q_%d",0), NbinsZ, -theDetector->ScintHalfZ, theDetector->ScintHalfZ,NbinsX, -theDetector->ScintHalfX, theDetector->ScintHalfX);
  YXProj_Q = new TH2F(Form("YXProj_Q_%d",0), Form("YXProj_Q_%d",0), NbinsY, -theDetector->ScintHalfY, theDetector->ScintHalfY,NbinsX, -theDetector->ScintHalfX, theDetector->ScintHalfX);
  YZProj_Q = new TH2F(Form("YZProj_Q_%d",0), Form("YZProj_Q_%d",0), NbinsY, -theDetector->ScintHalfY, theDetector->ScintHalfY,NbinsZ, -theDetector->ScintHalfZ, theDetector->ScintHalfZ);

}

void Analysis::RearFrontDetector(G4Step* aStep, G4String theName){
  f1->cd();
  if(theName=="FrontTracker"){
    x0  = aStep->GetPreStepPoint()->GetPosition()[0];
    y0  = aStep->GetPreStepPoint()->GetPosition()[1];
    z0  = aStep->GetPreStepPoint()->GetPosition()[2];

    px0 = aStep->GetPreStepPoint()->GetMomentumDirection()[0];
    py0 = aStep->GetPreStepPoint()->GetMomentumDirection()[1];
    pz0 = aStep->GetPreStepPoint()->GetMomentumDirection()[2];
    Einit = theConfig->item_float["Energy"];

  }
  else if(theName=="RearTracker"){
    x1  = aStep->GetPreStepPoint()->GetPosition()[0];
    y1  = aStep->GetPreStepPoint()->GetPosition()[1];
    z1  = aStep->GetPreStepPoint()->GetPosition()[2];

    px1 = aStep->GetPreStepPoint()->GetMomentumDirection()[0];
    py1 = aStep->GetPreStepPoint()->GetMomentumDirection()[1];
    pz1 = aStep->GetPreStepPoint()->GetMomentumDirection()[2];
    Estop = aStep->GetPreStepPoint()->GetKineticEnergy();
    if(aStep->GetTrack()->GetCreatorProcess()!=0) proc_name = aStep->GetTrack()->GetCreatorProcess()->GetProcessName();
    else proc_name  = "primary";
    part_name  = aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName();
    idPBY = theGenerator->idPBY;
    idPBZ = theGenerator->idPBZ;
    if(part_name.compare("gamma")!=0 && part_name.compare("neutron")!=0 && part_name.compare("e-")!=0){
      double WET = findWET(Einit, Estop);
      Front->Fill(y0,z0,WET);
      Back->Fill(y1,z1,WET);
      if(theConfig->item_int["saveTTree"]) t->Fill();
    }
  }
}

void Analysis::FillScintillatorDose(G4Step* aStep)
{
  part_name  = aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName();
  if(part_name.compare("gamma")!=0 && part_name.compare("neutron")!=0 && part_name.compare("e-")!=0){
    x_scint      = aStep->GetPostStepPoint()->GetPosition()[0] - (theDetector->ScintPosX); // substract the reference 0
    y_scint      = aStep->GetPostStepPoint()->GetPosition()[1];
    z_scint      = aStep->GetPostStepPoint()->GetPosition()[2];
    Estop_scint  = aStep->GetTotalEnergyDeposit();
    dX           = aStep->GetStepLength();

    if(aStep->GetStepLength()>0) LET = (aStep->GetTotalEnergyDeposit() - aStep->GetNonIonizingEnergyDeposit())/dX;
    else LET = 0;

    //Quenching factor
    A    = 1.0;
    kB   = 1.59E-2; // cm/MeV
    L    = dX*A*LET/(1+ kB*LET);
    // For particles count
    if(theConfig->item_int["particleCount"]){
      if(aStep->GetTrack()->GetTrackStatus() == fStopAndKill){ // End of Track
	if(aStep->GetTrack()->GetCreatorProcess()==0){ // Primary Particle
	  //PDD_Q->Fill(x_scint,1); // -- Non Quenched
	  //Distal projection beam by beam
	  if(theConfig->item_int["saveYXProj"]){ YXProj_Q->Fill(y_scint,x_scint,1);}
	  if(theConfig->item_int["saveZXProj"]){ ZXProj_Q->Fill(z_scint,x_scint,1);}
	  if(theConfig->item_int["saveYZProj"]){ YZProj_Q->Fill(z_scint,y_scint,1);}
	}
      }
    }
    else{
      //Lateral/Distal projection beam by beam 2-D (light)
      //PDD_Q->Fill(x_scint,L); // Quenched
      if(theConfig->item_int["saveYXProj"]){ YXProj_Q->Fill(y_scint,x_scint,L);}
      if(theConfig->item_int["saveZXProj"]){ ZXProj_Q->Fill(z_scint,x_scint,L);}
      if(theConfig->item_int["saveYZProj"]){ YZProj_Q->Fill(z_scint,y_scint,L);}
    }

    //Lateral/Distal projection beam by beam 2-D (energy)
    //PDD->Fill(x_scint,Estop_scint);
    //YXProj->Fill(y_scint,x_scint,Estop_scint);
    //ZXProj->Fill(z_scint,x_scint,Estop_scint);
    //YZProj->Fill(z_scint,x_scint,Estop_scint);

    //Cumulative Projections
    /*
    XYProj_Tot->Fill(x_scint,y_scint,Estop_scint);
    XZProj_Tot->Fill(x_scint,z_scint,Estop_scint);
    YZProj_Tot->Fill(y_scint,z_scint,Estop_scint);

    XYProj_Tot_Q->Fill(x_scint,y_scint,L);
    XZProj_Tot_Q->Fill(x_scint,z_scint,L);
    YZProj_Tot_Q->Fill(y_scint,z_scint,L);
    */

    // 3-D distribution Pencil beam by pencil beam
    //binGlobal = Edep_Tot->FindBin(x_scint,y_scint,z_scint); // Global bin associated with this voxel x-y-z
    // pair<float,int> dict1 = pair<float,int>(LET, 1);
    //pair<float,pair<float,int>> dict = pair<float,pair<float,int>>(Estop_scint,dict1); //find my key
    //ret = Edep_PB.insert(pair<int,pair<float,pair<float,int>>>(binGlobal,dict));
    /*if ( !ret.second ) {
      Edep_PB[binGlobal].first         += Estop_scint;
      Edep_PB[binGlobal].second.first  += LET ;
      Edep_PB[binGlobal].second.second += 1 ;
      }*/

    LET_Tot->Fill(x_scint,y_scint,z_scint,LET);
    Edep_Tot->Fill(x_scint, y_scint, z_scint, Estop_scint);
    L_Tot->Fill(x_scint, y_scint, z_scint, L);
    Entries_Tot->Fill(x_scint, y_scint, z_scint);
  }
}

void Analysis::SaveAndReset(){ // Save current PB histograms, close it and open the next one

  // PDD Pencil Beam By Pencil Beam (Energy)
  /*
  f1->cd("PDDList");
  PDD->Write("",TObject::kOverwrite);
  f1->cd();

  // 2-D Projection Pencil Beam By Pencil Beam (Energy)
  f1->cd("YZProj");
  YZProj->Write("",TObject::kOverwrite);
  f1->cd();

  // 2-D Projection Pencil Beam By Pencil Beam (Energy)
  f1->cd("YXProj");
  YXProj->Write("",TObject::kOverwrite);
  f1->cd();

  // 2-D Projection Pencil Beam By Pencil Beam (Energy)
  f1->cd("ZXProj");
  ZXProj->Write("",TObject::kOverwrite);
  f1->cd();


  // PDD Pencil Beam By Pencil Beam (Quenched Light)
  f1->cd("PDDList_Q");
  PDD_Q->Write("",TObject::kOverwrite);
  f1->cd();
  */

  //2-D Projection Pencil Beam By Pencil Beam (Quenched Light)
  if(theConfig->item_int["saveZXProj"]){
  f1->cd("ZXProj_Q");
  ZXProj_Q->Write("",TObject::kOverwrite);
  f1->cd();
  }

  // 2-D Projection Pencil Beam By Pencil Beam (Quenched Light)
  if(theConfig->item_int["saveYZProj"] ){
  f1->cd("YZProj_Q");
  YZProj_Q->Write("",TObject::kOverwrite);
  f1->cd();
  }

  // 2-D Projection Pencil Beam By Pencil Beam (Quenched Light)
  if(theConfig->item_int["saveYXProj"] ){
  f1->cd("YXProj_Q");
  YXProj_Q->Write("",TObject::kOverwrite);
  f1->cd();
  }

  // 3-D histogram Pencil Beam By Pencil Beam -- Sparse
  /*
  t2 = new TTree("PencilBeam","PB");
  t2->Branch("idPB",&idPB,"idPB/I"); // Pencil Beam ID
  t2->Branch("binglobal",&binGlobal,"binGlobal/I"); // Global Bin (voxel)
  t2->Branch("Estop_scint",&Estop_scint,"Estop_scint/F"); // Dose
  t2->Branch("N",&N,"N/I"); // LET
  t2->Branch("LET",&LET,"LET/F"); // LET
  for(int i =0; i<theGenerator->NPBY*theGenerator->NPBZ; i++){
    idPB = i;
    for(it = Edep_PB[i].begin(); it != Edep_PB[i].end(); it++) {
      binGlobal = it->first;
      Estop_scint  = it->second.first;
      LET          = it->second.second.first;
      N            = it->second.second.second;
      t2->Fill();
    }
  }
  t2->Write("",TObject::kOverwrite);
  */


  delete ZXProj_Q;delete YXProj_Q;delete YZProj_Q;


  /*PDD    = new TH1F(Form("PDD_Q_%d",theGenerator->idPBGlobal),Form("PDD_Q_%d",theGenerator->idPBGlobal),NbinsX, -theDetector->ScintHalfX, theDetector->ScintHalfX);
  ZXProj = new TH2F(Form("ZXProj_%d",theGenerator->idPBGlobal), Form("ZXProj_%d",theGenerator->idPBGlobal), NbinsZ, -theDetector->ScintHalfZ, theDetector->ScintHalfZ,NbinsX, -theDetector->ScintHalfX, theDetector->ScintHalfX);
  YXProj = new TH2F(Form("YXProj_%d",theGenerator->idPBGlobal), Form("YXProj_%d",theGenerator->idPBGlobal), NbinsY, -theDetector->ScintHalfY, theDetector->ScintHalfY,NbinsX, -theDetector->ScintHalfX, theDetector->ScintHalfX);
  YZProj = new TH2F(Form("YZProj_%d",theGenerator->idPBGlobal), Form("YZProj_%d",theGenerator->idPBGlobal), NbinsY, -theDetector->ScintHalfY, theDetector->ScintHalfY,NbinsZ, -theDetector->ScintHalfZ, theDetector->ScintHalfZ);  */

  //PDD_Q    = new TH1F(Form("PDD_Q_%d",theGenerator->idPBGlobal),Form("PDD_Q_%d",theGenerator->idPBGlobal),NbinsX, -theDetector->ScintHalfX, theDetector->ScintHalfX);

  ZXProj_Q = new TH2F(Form("ZXProj_Q_%d",theGenerator->idPBGlobal), Form("ZXProj_Q_%d",theGenerator->idPBGlobal), NbinsZ, -theDetector->ScintHalfZ, theDetector->ScintHalfZ,NbinsX, -theDetector->ScintHalfX, theDetector->ScintHalfX);
  YXProj_Q = new TH2F(Form("YXProj_Q_%d",theGenerator->idPBGlobal), Form("YXProj_Q_%d",theGenerator->idPBGlobal), NbinsY, -theDetector->ScintHalfY, theDetector->ScintHalfY,NbinsX, -theDetector->ScintHalfX, theDetector->ScintHalfX);
  YZProj_Q = new TH2F(Form("YZProj_Q_%d",theGenerator->idPBGlobal), Form("YZProj_Q_%d",theGenerator->idPBGlobal), NbinsY, -theDetector->ScintHalfY, theDetector->ScintHalfY,NbinsZ, -theDetector->ScintHalfZ, theDetector->ScintHalfZ);
}
void Analysis::SaveAndClose(){
  cout<<"Recording"<<endl;
  f1->cd();

  // Single Event
  if(theConfig->item_int["saveTTree"]) t->Write("",TObject::kOverwrite);

  // Full 3-D Histogram
  Edep_Tot->Write("",TObject::kOverwrite);
  Entries_Tot->Write("",TObject::kOverwrite);
  LET_Tot->Write("",TObject::kOverwrite);
  L_Tot->Write("",TObject::kOverwrite);

  /*
  // Integrated 2-D projection
  XYProj_Tot->Write("",TObject::kOverwrite);
  XZProj_Tot->Write("",TObject::kOverwrite);
  YZProj_Tot->Write("",TObject::kOverwrite);
  XYProj_Tot_Q->Write("",TObject::kOverwrite);
  XZProj_Tot_Q->Write("",TObject::kOverwrite);
  YZProj_Tot_Q->Write("",TObject::kOverwrite);
  */

  // Single event
  Front->Write("",TObject::kOverwrite);
  Back->Write("",TObject::kOverwrite);

  //Pencil Beam Position
  f1->cd();
  TTree* t3 = new TTree("PencilBeamPos","PB");
  float PBPosZ, PBPosY;
  t3->Branch("PBPosZ",&PBPosZ,"PBPosZ/F"); // Pencil Beam ID
  t3->Branch("PBPosY",&PBPosY,"PBPosY/F"); // Pencil Beam ID
  for(int i =0; i<theGenerator->NPBY*theGenerator->NPBZ; i++){
    int idY        = (i - (i%theGenerator->NPBY))/theGenerator->NPBY;
    int idZ        =  i%theGenerator->NPBY;
    PBPosZ = theGenerator->PencilBeamPosZ[idZ];
    PBPosY = theGenerator->PencilBeamPosY[idY];
    t3->Fill();
  }
  t3->Write("",TObject::kOverwrite);

  // Header containing relevant info
  TTree* t4 = new TTree("Header","");
  t4->Branch("NPB",&theConfig->item_int["NPB"],"NPB/I"); // Pencil Beam ID
  t4->Branch("NProtons",&theConfig->item_int["NProtons"],"NProtons/I"); // Pencil Beam ID
  t4->Branch("Energy",&theConfig->item_float["Energy"],"Energy/F"); // Pencil Beam ID
  t4->Branch("Thickness",&theConfig->item_float["Thickness"],"Thickness/I"); // Pencil Beam ID
  t4->Branch("ANumber",&theConfig->item_int["ANumber"],"ANumber/I"); // Pencil Beam ID
  t4->Branch("Angle",&theConfig->item_float["Angle"],"Angle/I"); // Pencil Beam ID
  t4->Branch("sigmaX_pos",&theConfig->item_float["sigmaX_pos"],"sigmaX_pos/F"); // Pencil Beam ID
  t4->Branch("sigmaY_pos",&theConfig->item_float["sigmaY_pos"],"sigmaY_pos/F"); // Pencil Beam ID
  t4->Branch("fieldSizeY",&theConfig->item_float["fieldSizeY"],"fieldSizeY/F"); // Pencil Beam ID
  t4->Branch("fieldSizeZ",&theConfig->item_float["fieldSizeZ"],"fieldSizeZ/F"); // Pencil Beam ID
  t4->Branch("sigma_AngX",&theConfig->item_float["sigma_AngX"],"sigma_AngX/F"); // Pencil Beam ID
  t4->Branch("sigma_AngY",&theConfig->item_float["sigma_AngY"],"sigma_AngY/F"); // Pencil Beam ID
  t4->Branch("centerY",&theConfig->item_float["centerY"],"centerY/F"); // Pencil Beam ID
  t4->Branch("centerZ",&theConfig->item_float["centerZ"],"centerZ/F"); // Pencil Beam ID

  //Strings
  t4->Branch("SourceType",&theConfig->item_str["SourceType"]); // Pencil Beam ID
  t4->Branch("Model",&theConfig->item_str["Model"]); // Pencil Beam ID
  t4->Branch("Phase",&theConfig->item_str["CTPath"]); // Pencil Beam ID

  t4->Fill();
  t4->Write("",TObject::kOverwrite);
  f1->Close();

}

double Analysis::findWET(double Ein,double Eout){
  int it_Einit = lower_bound(Energy.begin(), Energy.end(), Ein)-Energy.begin();
  int it_Estop = lower_bound(Energy.begin(), Energy.end(), Eout)-Energy.begin();
  double WET = 0 ;
  for(int i=it_Estop;i<it_Einit;i++){
    WET += 0.1/dEdXBins[i];
  }
  return WET;
}
