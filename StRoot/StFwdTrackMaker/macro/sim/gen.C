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


TH1F* hNumHits = 0;
TString nameParticle = "mu+";
int numParticles = 1;
float minPt = 0.0;
float maxPt = 1.0;
float minEta = 2.5;
float maxEta = 4.00;
float minPhi = 0.0;
float maxPhi = 2.0 * TMath::Pi();

float vtxX = 0.0;
float vtxY = 0.0;
float vtxZ = 0.0;

float vtxSigmaX = 0.0001;
float vtxSigmaY = 0.0001;
float vtxSigmaZ = 0.0001;

TString fzdFilename = "sim.fzd";
TString primaryName = "sim.root";

// Geometry tag and DB lane for the GENERATOR. These were hardcoded as
// "y2024" + "sdt20211016" in two places each; they are globals now so a caller
// can put generation and reconstruction on the same lane, which is the only
// way the two agree.
//
// The default is y2022c, not y2024, on purpose: y2022c is the tag PR28 wired to
// FSTMv2m, i.e. fstmconfig=2, the as-installed FST wedge mounting, which is
// what StFstDb::fstmConfig() returns. Generating on y2024 and reconstructing
// with a y2024 cache trips StFwdTrackMaker's mounting-consistency check --
// "FST geometry and StFstDb disagree about the wedge mounting" -- because the
// y2024 AGML is still on the other phase. Do not change this back without
// checking that error is absent from the reco log.
//
// dev2022sm is NOT an option here: it is the SIMPLIFIED geometry and produces
// no FST g2t hits at all, so a sample generated on it has nothing to track.
TString genGeomTag = "y2022c";
TString genSdt     = "sdt20211112.000004";

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
void trig_event( Int_t i )
{
  kinematics->Kine( numParticles, nameParticle.Data(), minPt, maxPt, minEta, maxEta, minPhi, maxPhi );
}
// ----------------------------------------------------------------------------
void trig( Int_t n=1 )
{
  for ( Int_t i=0; i<n; i++ ) {
    // Clear the chain from the previous event
    chain->Clear();
    trig_event( i );
    // Generate the event
    chain->Make();
  }
}
// ----------------------------------------------------------------------------
void Kinematics()
{

  //  gSystem->Load( "libStarGeneratorPoolPythia6_4_23.so" );
  gSystem->Load( "libKinematics.so");
  kinematics = new StarKinematics();
  _primary->AddGenerator(kinematics);
}
// ----------------------------------------------------------------------------
void gen( Int_t nevents=1000, Int_t rngSeed=12352342 )
{

  cout << "Generating: " << nevents << " events with seed: " << rngSeed << endl;
  gSystem->Load( "libStarRoot.so" );
  // The macro path used to be a hardcoded list ending in a dozen
  // /afs/rhic.bnl.gov/... directories. /afs is gone -- inside the container it
  // is an empty bind mount -- and because SetMacroPath REPLACES the path
  // rather than extending it, this threw away the working $STAR entry too.
  // Every gun built on gen.C therefore died at the first line that matters:
  //   Error in <TRint::ProcessLine>: macro bfc.C not found in path ...
  //   Error: Function bfc(0,simple) is not defined in current scope gen.C
  // Build it from $STAR instead, and keep whatever was already there.
  gROOT->SetMacroPath( Form(".:./StRoot/macros:%s/StRoot/macros:%s",
                            gSystem->Getenv("STAR"), gROOT->GetMacroPath()) );

  gROOT->ProcessLine(".L bfc.C");
  {
    TString simple = genSdt + " " + genGeomTag + " geant gstar usexgeom agml ";
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
    _primary -> SetFileName( primaryName );
    chain -> AddBefore( "geant", _primary );
  }

  Kinematics();

  //
  // Initialize primary event generator and all sub makers
  //
  _primary -> Init();
  _primary->SetSigma( vtxSigmaX, vtxSigmaY, vtxSigmaZ ); // 1mm x 1mm x 1mm smearing at the vertex
  _primary->SetVertex(vtxX, vtxY, vtxZ );

  //
  // Setup geometry and set starsim to use agusread for input
  //
  geometry(genGeomTag);   // DETP GEOM ... field=: fz header records the field (reversed full, as in 2022+ data)
  command("gkine -4 0");
  command( TString::Format("gfile o %s", fzdFilename.Data()) );


  // command( "DCAY 0" );
  // command( "ANNI 0" );
  // command( "BREM 0" );
  // command( "COMP 0" );
  // command( "HADR 0" );
  // command( "MUNU 0" );
  // command( "PAIR 0" );
  // command( "PFIS 0" );
  // command( "PHOT 0" );
  // command( "RAYL 0" );
  // command( "LOSS 4" );
  // command( "DRAY 0" );
  // command( "MULS 0" );
  // command( "STRA 0" );
  // command( "physi"  );

  //
  // Trigger on nevents
  //
  // StarMagField::setConstBz(true);
  trig( nevents );

  command("call agexit");  // Make sure that STARSIM exits properly

}
// ----------------------------------------------------------------------------
