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

Analysis::Analysis(G4int thread, G4double angle,G4String theName){
  theAnalysis       = this;
  theGenerator      = PrimaryGeneratorAction::GetInstance();
  theDetector       = DetectorConstruction::GetInstance();
  f1 = new TFile(Form("%s_%.0f_%.1f_%d_%d.root",theName.data(),theGenerator->ENER,angle,thread,theGenerator->A),"recreate");

  NbinsX = 150; NbinsY = 150; NbinsZ = 150;

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
  XYProj_Tot  = new TH2F("XYProj_Tot", "XYProj_Tot", NbinsX, -theDetector->ScintHalfX, theDetector->ScintHalfX, 
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

  //Phasespace
  /*t = new TTree("phase","PS");
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
  */
  //Pencil beam by pencil beam event
  for(int i =0;i<theGenerator->NPBY*theGenerator->NPBZ;i++){

    // Non Quenched
    /*
    PDD[i]      = new TH1D(Form("PDD_%d",i),Form("PDD_%d",i),NbinsX, -theDetector->ScintHalfX, theDetector->ScintHalfX);
    ZXProj[i]   = new TH2F(Form("ZXProj_%d",i), Form("ZXProj_%d",i), NbinsZ, -theDetector->ScintHalfZ, theDetector->ScintHalfZ,
    				 NbinsX, -theDetector->ScintHalfX, theDetector->ScintHalfX);
    YXProj[i]   = new TH2F(Form("YXProj_%d",i), Form("YXProj_%d",i), NbinsY, -theDetector->ScintHalfY, theDetector->ScintHalfY,
    				 NbinsX, -theDetector->ScintHalfX, theDetector->ScintHalfX);
    YZProj[i]   = new TH2F(Form("YZProj_%d",i), Form("YZProj_%d",i), NbinsY, -theDetector->ScintHalfY, theDetector->ScintHalfY,
    				 NbinsX, -theDetector->ScintHalfX, theDetector->ScintHalfX);
    */
    //Quenched

    /*PDD_Q[i]      = new TH1F(Form("PDD_Q_%d",i),Form("PDD_Q_%d",i),NbinsX, -theDetector->ScintHalfX, theDetector->ScintHalfX);
    ZXProj_Q[i]   = new TH2F(Form("ZXProj_Q_%d",i), Form("ZXProj_Q_%d",i), NbinsZ, -theDetector->ScintHalfZ, theDetector->ScintHalfZ,
			     NbinsX, -theDetector->ScintHalfX, theDetector->ScintHalfX);
    YXProj_Q[i]   = new TH2F(Form("YXProj_Q_%d",i), Form("YXProj_Q_%d",i), NbinsY, -theDetector->ScintHalfY, theDetector->ScintHalfY,
			     NbinsX, -theDetector->ScintHalfX, theDetector->ScintHalfX);
    YZProj_Q[i]   = new TH2F(Form("YZProj_Q_%d",i), Form("YZProj_Q_%d",i), NbinsY, -theDetector->ScintHalfY, theDetector->ScintHalfY,
			     NbinsX, -theDetector->ScintHalfX, theDetector->ScintHalfX);
    */
    PDD_Q.push_back(new TH1F(Form("PDD_Q_%d",i),Form("PDD_Q_%d",i),NbinsX, -theDetector->ScintHalfX, theDetector->ScintHalfX));
    ZXProj_Q.push_back(new TH2F(Form("ZXProj_Q_%d",i), Form("ZXProj_Q_%d",i), NbinsZ, -theDetector->ScintHalfZ, theDetector->ScintHalfZ,
				NbinsX, -theDetector->ScintHalfX, theDetector->ScintHalfX));
    YXProj_Q.push_back(new TH2F(Form("YXProj_Q_%d",i), Form("YXProj_Q_%d",i), NbinsY, -theDetector->ScintHalfY, theDetector->ScintHalfY,
				NbinsX, -theDetector->ScintHalfX, theDetector->ScintHalfX));
    YZProj_Q.push_back(new TH2F(Form("YZProj_Q_%d",i), Form("YZProj_Q_%d",i), NbinsY, -theDetector->ScintHalfY, theDetector->ScintHalfY,
				NbinsX, -theDetector->ScintHalfX, theDetector->ScintHalfX));

  }  
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
    Einit = theGenerator->Einit; 

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
    }
  }
}

void Analysis::FillScintillatorDose(G4Step* aStep)
{
  part_name  = aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName();
  if(part_name.compare("gamma")!=0 && part_name.compare("neutron")!=0 && part_name.compare("e-")!=0){
    theDetector = DetectorConstruction::GetInstance();
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
    /*
    if(aStep->GetTrack()->GetTrackStatus() == fStopAndKill){ // End of Track
      if(aStep->GetTrack()->GetCreatorProcess()==0){ // Primary Particle
	// For the pencil beams  
	//PDD[theGenerator->idPBGlobal]->Fill(x_scint,1); // -- Non Quenched
	//Distal projection beam by beam
	//YXProj[theGenerator->idPBGlobal]->Fill(y_scint,x_scint,1);
	//ZXProj[theGenerator->idPBGlobal]->Fill(z_scint,x_scint,1);
	//YZProj[theGenerator->idPBGlobal]->Fill(z_scint,x_scint,1);	
	}
      }
    */
    
    //Lateral/Distal projection beam by beam 2-D (energy)
    //PDD[theGenerator->idPBGlobal]->Fill(x_scint,Estop_scint); 
    //YXProj[theGenerator->idPBGlobal]->Fill(y_scint,x_scint,Estop_scint);
    //ZXProj[theGenerator->idPBGlobal]->Fill(z_scint,x_scint,Estop_scint);
    //YZProj[theGenerator->idPBGlobal]->Fill(z_scint,x_scint,Estop_scint);    

    //Lateral/Distal projection beam by beam 2-D (light)
    PDD_Q[theGenerator->idPBGlobal]->Fill(x_scint,L); // Quenched    
    YXProj_Q[theGenerator->idPBGlobal]->Fill(y_scint,x_scint,L);
    ZXProj_Q[theGenerator->idPBGlobal]->Fill(z_scint,x_scint,L);
    YZProj_Q[theGenerator->idPBGlobal]->Fill(z_scint,y_scint,L);    
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
    //ret = Edep_PB[theGenerator->idPBGlobal].insert(pair<int,pair<float,pair<float,int>>>(binGlobal,dict));
    /*if ( !ret.second ) {
      Edep_PB[theGenerator->idPBGlobal][binGlobal].first         += Estop_scint;
      Edep_PB[theGenerator->idPBGlobal][binGlobal].second.first  += LET ;
      Edep_PB[theGenerator->idPBGlobal][binGlobal].second.second += 1 ;    
      }*/

    
    LET_Tot->Fill(x_scint,y_scint,z_scint,LET);     
    Edep_Tot->Fill(x_scint, y_scint, z_scint, Estop_scint);
    L_Tot->Fill(x_scint, y_scint, z_scint, L);    
    Entries_Tot->Fill(x_scint, y_scint, z_scint);
  }
}
void Analysis::Save(){
  cout<<"Recording"<<endl;
  f1->cd();
  // Single Event
  //t->Write("",TObject::kOverwrite);  

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

  // PDD Pencil Beam By Pencil Beam (Energy)
  /*
  f1->mkdir("PDDList");
  f1->cd("PDDList");
  for(int i =0;i<theGenerator->NPBY*theGenerator->NPBZ;i++) PDD[i]->Write("",TObject::kOverwrite);
  f1->cd();  

  // 2-D Projection Pencil Beam By Pencil Beam (Energy)
  f1->mkdir("YZProj");
  f1->cd("YZProj");
  for(int i =0;i<theGenerator->NPBY*theGenerator->NPBZ;i++) YZProj[i]->Write("",TObject::kOverwrite);
  f1->cd();

  // 2-D Projection Pencil Beam By Pencil Beam (Energy)
  f1->mkdir("YXProj");
  f1->cd("YXProj");
  for(int i =0;i<theGenerator->NPBY*theGenerator->NPBZ;i++) YXProj[i]->Write("",TObject::kOverwrite);
  f1->cd();

  // 2-D Projection Pencil Beam By Pencil Beam (Energy)
  f1->mkdir("ZXProj");
  f1->cd("ZXProj");
  for(int i =0;i<theGenerator->NPBY*theGenerator->NPBZ;i++) ZXProj[i]->Write("",TObject::kOverwrite);
  f1->cd();  
  */

  // PDD Pencil Beam By Pencil Beam (Quenched Light)
  f1->mkdir("PDDList_Q");
  f1->cd("PDDList_Q");
  for(int i =0;i<theGenerator->NPBY*theGenerator->NPBZ;i++) PDD_Q[i]->Write("",TObject::kOverwrite);
  f1->cd();  
  
  //2-D Projection Pencil Beam By Pencil Beam (Quenched Light)
  f1->mkdir("ZXProj_Q");
  f1->cd("ZXProj_Q");
  for(int i =0;i<theGenerator->NPBY*theGenerator->NPBZ;i++) ZXProj_Q[i]->Write("",TObject::kOverwrite);
  f1->cd();
  
  // 2-D Projection Pencil Beam By Pencil Beam (Quenched Light)

  f1->mkdir("YZProj_Q");
  f1->cd("YZProj_Q");
  for(int i =0;i<theGenerator->NPBY*theGenerator->NPBZ;i++) YZProj_Q[i]->Write("",TObject::kOverwrite);
  f1->cd();

  
  // 2-D Projection Pencil Beam By Pencil Beam (Quenched Light)
  f1->mkdir("YXProj_Q");
  f1->cd("YXProj_Q");
  for(int i =0;i<theGenerator->NPBY*theGenerator->NPBZ;i++) YXProj_Q[i]->Write("",TObject::kOverwrite);
  f1->cd();

  
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
  int NGlobal = theGenerator->NPBY*theGenerator->NPBZ;
  TTree* t4 = new TTree("Header","");
  t4->Branch("NPB",&NGlobal,"NPB/I"); // Pencil Beam ID
  t4->Branch("PB_Sigma_Z",&theGenerator->PencilBeamStdZ,"PB_Sigma_Z/D"); // Pencil Beam ID   
  t4->Branch("PB_Sigma_Y",&theGenerator->PencilBeamStdY,"PB_Sigma_Y/D"); // Pencil Beam ID
  t4->Branch("PB_Sigma_Ang",&theGenerator->PencilBeamStdAng,"PB_Sigma_Ang/D"); // Pencil Beam ID  
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
