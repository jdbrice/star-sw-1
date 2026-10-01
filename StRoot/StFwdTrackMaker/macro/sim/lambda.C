//usr/bin/env root4star -l -b -q $0'('$1', '$2')'; exit $?
// macro to instantiate the Geant3 from within
// STAR  C++  framework and get the starsim prompt
// To use it do
//  root4star starsim.C

class St_geant_Maker;
St_geant_Maker *geant_maker = 0;

class StarGenEvent;
StarGenEvent   *event       = 0;

class StarPrimaryMaker;
StarPrimaryMaker *_primary = 0;

class StarKinematics;
StarKinematics *kinematics = 0;


TH1F* hM = 0;
float numParticles = 1;

// ----------------------------------------------------------------------------
void geometry( TString tag, Bool_t agml=true )
{
  TString cmd = "DETP GEOM "; cmd += tag + " field=-5.0";
  if ( !geant_maker ) geant_maker = (St_geant_Maker *)chain->GetMaker("geant");
  geant_maker -> LoadGeometry(cmd);
  //  if ( agml ) command("gexec $STAR_LIB/libxgeometry.so");
}
// ----------------------------------------------------------------------------
void command( TString cmd )
{
  if ( !geant_maker ) geant_maker = (St_geant_Maker *)chain->GetMaker("geant");
  geant_maker -> Do( cmd );
}
// ----------------------------------------------------------------------------
void trig( Int_t n=1 )
{

  
  for ( Int_t i=0; i<n; i++ ) {

    // Clear the chain from the previous event
    chain->Clear();

    //(Momentum, Energy units are Gev/C, GeV)
    Double_t masses[2] = { 0.13957039, 0.93827208943};

    TGenPhaseSpace genEvent;
    TLorentzVector W;
    // W.SetPtEtaPhiM( 0.0, 100.0, 0, 3.096 );
    double px = gRandom->Uniform(-1.0, 1.0);
    double py = gRandom->Uniform(-1.0, 1.0);
    double pz = gRandom->Uniform(0, 50.0);
    W.SetXYZM( px, py, pz, 1.115683 );
    printf("W: px = %f, py = %f, pz = %f, m = %f\n", W.Px(), W.Py(), W.Pz(), W.M() );
    genEvent.SetDecay(W, 2, masses);

    TLorentzVector lv;
    for ( int j = 0; j < numParticles; j++ ){
      Double_t weight = genEvent.Generate();
      TLorentzVector lvPion = *(genEvent.GetDecay(0));
      TLorentzVector lvProton = *(genEvent.GetDecay(1));
      lv = lvPion + lvProton;

      StarGenParticle *pion;
      pion = kinematics->AddParticle( "pi-" );
      pion->SetPx(lvPion.Px());
      pion->SetPy(lvPion.Py());
      pion->SetPz(lvPion.Pz());
      pion->SetMass( masses[0] );

      // doesnt seem to do anything!
      // pion->SetVx(1.0);
      // pion->SetVy(1.0);
      // pion->SetVz(10.0);

      StarGenParticle *proton;
      proton = kinematics->AddParticle( "p" );
      proton->SetPx(lvProton.Px());
      proton->SetPy(lvProton.Py());
      proton->SetPz(lvProton.Pz());
      proton->SetMass( masses[1] );

      // doesnt seem to do anything!
      // proton->SetVx(1.0);
      // proton->SetVy(1.0);
      // proton->SetVz(10.0);

      hM->Fill( lv.M() );

      cout << "pion eta = " << lvPion.Eta() << endl;
      cout << "proton eta = " << lvProton.Eta() << endl;
    }

    
		// kinematics->Kine( numParticles, nameParticle.Data(), 10.2, 12.0, 2.5, 4.00  );

    // Generate the event
    chain->Make();

    // Print the event
    cout << "Event: " << i << " generated with " << numParticles << " Lambda." << endl;
  }
}
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
void Kinematics()
{
  
  //  gSystem->Load( "libStarGeneratorPoolPythia6_4_23.so" );
  gSystem->Load( "libKinematics.so");
  kinematics = new StarKinematics();
    
  _primary->AddGenerator(kinematics);
}
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
void lambda( Int_t nevents=10000, Int_t rngSeed=0 )
{ 

  hM = new TH1F("hM",";M_{p#pi};counts [10MeV]", 50, 1.0, 1.5 );
  cout << "Generating: " << nevents << " events with seed: " << rngSeed << endl;
  gSystem->Load( "libStarRoot.so" );

  gROOT->ProcessLine(".L bfc.C");
  {
    TString simple = "y2024 geant gstar usexgeom agml ";
    bfc(0, simple );
  }

  gSystem->Load( "libVMC.so");

  gSystem->Load( "StarGeneratorUtil.so" );
  gSystem->Load( "StarGeneratorEvent.so" );
  gSystem->Load( "StarGeneratorBase.so" );

  gSystem->Load( "libMathMore.so"   );  
  gSystem->Load( "xgeometry.so"     );

  

  // Setup RNG seed and map all ROOT TRandom here
  StarRandom::seed( rngSeed );
  StarRandom::capture();

  //
  // Create the primary event generator and insert it
  // before the geant maker
  //
  //  StarPrimaryMaker *
  _primary = new StarPrimaryMaker();
  {
    _primary -> SetFileName( "lambda.root");
    chain -> AddBefore( "geant", _primary );
  }

  Kinematics();

  //
  // Initialize primary event generator and all sub makers
  //
  _primary -> Init();
  _primary->SetSigma( 0.001, 0.001, 0.001 );
  _primary->SetVertex(1.0, 1.0, 50.0 );

  //
  // Setup geometry and set starsim to use agusread for input
  //
  geometry("y2024");   // DETP GEOM ... field=: fz header records the field (reversed full, as in 2022+ data)
  command("gkine -4 0");
  command("gfile o fwd_lambda.fzd");
 
  //
  // Trigger on nevents
  //
  trig( nevents );

  TFile * f = new TFile( "fwd_lambda_gen.root", "RECREATE" );
  f->cd();
  hM->Write();
  f->Write();

  command("call agexit");  // Make sure that STARSIM exits properly

}
// ----------------------------------------------------------------------------

