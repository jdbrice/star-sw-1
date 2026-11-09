//usr/bin/env root4star -l -b -q $0'('$1', '$2')'; exit $?
//
// lambda.C -- Lambda -> p pi- into the FORWARD acceptance, as input for forward
//             tracking and GenFit/RAVE secondary-vertex finding.
//
//   root4star -l -b -q 'lambda.C(10)'                  10 events, "gun" mode
//   root4star -l -b -q 'lambda.C(10,20261109,"ps")'    the old phase-space mode
//
// Output: fwd_lambda_<mode>.fzd   (feed this to sim.C)
//         fwd_lambda_<mode>.root  the StarPrimaryMaker generator record
//         fwd_lambda_gen.root     hM, the generated p-pi invariant mass ("ps")
//
// TWO MODES, and the difference is the whole point:
//
//   "gun" (default)  one Lambda0 per event, thrown flat in the forward
//                    acceptance and handed to GEANT3, which transports it and
//                    DECAYS it itself. The decay vertex is therefore REAL and
//                    DISPLACED: c*tau = 7.89 cm, so <L> = (p/m)*7.89 cm, i.e.
//                    ~7 cm per GeV of Lambda momentum. A 10 GeV Lambda flies
//                    ~70 cm before decaying and its daughters start well
//                    downstream of the primary vertex. This is the sample a
//                    secondary-vertex study needs.
//                    GEANT uses the real branching ratios, ~64% p pi- and
//                    ~36% n pi0, so only about two events in three give the
//                    two charged forward tracks RAVE can make a vertex from.
//
//   "ps"             what this macro used to do: TGenPhaseSpace decays the
//                    Lambda HERE and the proton and the pi- are pushed as
//                    PRIMARIES, so both start AT the primary vertex and the
//                    decay length is exactly zero. Useful as a control -- it
//                    shows what RAVE does with a pair that genuinely comes
//                    from one point -- but it is not a Lambda sample, and it
//                    cannot tell you anything about displaced vertices.
//
// WHAT CHANGED 2026-11-09, and why the old version could not be used as it was
//
//   1. The parent was thrown as px,py in (-1,1) and pz in (0,50) uniformly.
//      That is a box in momentum, not an acceptance: pz near zero puts the
//      Lambda at eta ~ 0, nowhere near the forward detectors. Now it is thrown
//      flat in eta over [minEta,maxEta] = [2.5,4.0], which is what the FST and
//      the sTGC actually cover.
//   2. The daughters never had their ENERGY set. StarPrimaryMaker reads
//      GetEnergy() and hands it straight to StarParticleStack::PushTrack, so
//      every particle this macro produced went into GEANT with E = 0 while
//      carrying non-zero momentum. StarKinematics::Kine() sets the energy;
//      this hand-rolled AddParticle path did not.
//   3. The primary vertex was (1,1,50) cm with 10 um sigma. 50 cm downstream
//      is not where the beams cross; it is now (0,0,0), settable.
//   4. bfc got no "sdt", so generation read whatever DB lane today's date
//      selects, while reco (sim.C) names its lane explicitly. The two must
//      MATCH or GEANT and the tracker use different geometry. It is now an
//      argument, defaulting to the same sdt20240622 that sim.C uses.
//
// Note that eta is thrown, not rapidity, and the energy is then built as
// sqrt(pt^2 + pz^2 + m^2). Do not be tempted to replace this with
// StarKinematics::Kine(): without the "rapidity" attribute set it treats its
// y argument as pseudorapidity but computes E = pt*cosh(eta), which is the
// MASSLESS energy -- harmless for a muon gun, wrong by ~60 MeV for a Lambda,
// and below the mass itself at low pt.

class St_geant_Maker;
St_geant_Maker *geant_maker = 0;

class StarGenEvent;
StarGenEvent   *event       = 0;

class StarPrimaryMaker;
StarPrimaryMaker *_primary = 0;

class StarKinematics;
StarKinematics *kinematics = 0;

TH1F *hM = 0;   // generated p-pi mass, "ps" mode only
TH1F *hEta = 0; // generated Lambda eta
TH1F *hP = 0;   // generated Lambda momentum
TH1F *hL = 0;   // mean decay length beta*gamma*c*tau implied by that momentum

// --- what to throw -----------------------------------------------------------
int   nLambda = 1;        // Lambdas per event
float minPt   = 0.3;      // GeV/c
float maxPt   = 1.5;
float minEta  = 2.5;      // the forward acceptance
float maxEta  = 4.0;
float minPhi  = 0.0;
float maxPhi  = 2.0 * TMath::Pi();

// --- where to put the primary vertex -----------------------------------------
float vtxX = 0.0, vtxY = 0.0, vtxZ = 0.0;             // cm
float vtxSigmaX = 0.01, vtxSigmaY = 0.01, vtxSigmaZ = 0.01;  // cm, 100 um

const double kMLambda = 1.115683;
const double kMProton = 0.93827208943;
const double kMPion   = 0.13957039;
const double kCTau     = 7.89;   // cm, Lambda

TString gMode = "gun";
TString gFzd  = "fwd_lambda_gun.fzd";

// ----------------------------------------------------------------------------
void geometry( TString tag, Bool_t agml=true )
{
  TString cmd = "DETP GEOM "; cmd += tag + " field=-5.0";
  if ( !geant_maker ) geant_maker = (St_geant_Maker *)chain->GetMaker("geant");
  geant_maker -> LoadGeometry(cmd);
}
// ----------------------------------------------------------------------------
void command( TString cmd )
{
  if ( !geant_maker ) geant_maker = (St_geant_Maker *)chain->GetMaker("geant");
  geant_maker -> Do( cmd );
}
// ----------------------------------------------------------------------------
// Throw one Lambda flat in pt, eta, phi and return it as a TLorentzVector.
TLorentzVector throwLambda()
{
  double pt  = gRandom->Uniform( minPt,  maxPt  );
  double eta = gRandom->Uniform( minEta, maxEta );
  double phi = gRandom->Uniform( minPhi, maxPhi );

  double px = pt * TMath::Cos( phi );
  double py = pt * TMath::Sin( phi );
  double pz = pt * TMath::SinH( eta );
  double E  = TMath::Sqrt( px*px + py*py + pz*pz + kMLambda*kMLambda );

  TLorentzVector W;
  W.SetPxPyPzE( px, py, pz, E );
  return W;
}
// ----------------------------------------------------------------------------
void trig( Int_t n = 1 )
{
  for ( Int_t iev = 0; iev < n; iev++ ) {

    chain->Clear();

    for ( int il = 0; il < nLambda; il++ ) {

      TLorentzVector W = throwLambda();

      // beta*gamma*c*tau, the mean decay length this Lambda will fly. GEANT
      // samples an exponential around it; this is only the scale.
      double bgctau = ( W.P() / kMLambda ) * kCTau;

      hEta->Fill( W.Eta() );
      hP  ->Fill( W.P()   );
      hL  ->Fill( bgctau  );

      printf( "ev %3d  Lambda pt=%6.3f eta=%5.3f p=%7.3f GeV/c  <L>=%7.2f cm\n",
              iev, W.Pt(), W.Eta(), W.P(), bgctau );

      if ( gMode == "gun" ) {

        // Hand the Lambda itself to GEANT and let IT decay the particle, so
        // the daughters are born at a real, displaced decay vertex.
        StarGenParticle *lam = kinematics->AddParticle( "Lambda0" );
        lam->SetPx( W.Px() );
        lam->SetPy( W.Py() );
        lam->SetPz( W.Pz() );
        lam->SetEnergy( W.E() );
        lam->SetMass( kMLambda );

      } else {

        // "ps": decay it here and push the daughters as primaries. Both start
        // at the primary vertex -- zero decay length, by construction.
        Double_t masses[2] = { kMPion, kMProton };
        TGenPhaseSpace genEvent;
        genEvent.SetDecay( W, 2, masses );
        genEvent.Generate();

        TLorentzVector lvPion   = *(genEvent.GetDecay(0));
        TLorentzVector lvProton = *(genEvent.GetDecay(1));
        hM->Fill( (lvPion + lvProton).M() );

        StarGenParticle *pion = kinematics->AddParticle( "pi-" );
        pion->SetPx( lvPion.Px() );
        pion->SetPy( lvPion.Py() );
        pion->SetPz( lvPion.Pz() );
        pion->SetEnergy( lvPion.E() );
        pion->SetMass( kMPion );

        StarGenParticle *proton = kinematics->AddParticle( "p" );
        proton->SetPx( lvProton.Px() );
        proton->SetPy( lvProton.Py() );
        proton->SetPz( lvProton.Pz() );
        proton->SetEnergy( lvProton.E() );
        proton->SetMass( kMProton );

        printf( "          -> pi- eta=%5.3f pt=%6.3f   p eta=%5.3f pt=%6.3f\n",
                lvPion.Eta(), lvPion.Pt(), lvProton.Eta(), lvProton.Pt() );
      }
    }

    chain->Make();
  }
}
// ----------------------------------------------------------------------------
void Kinematics()
{
  gSystem->Load( "libKinematics.so" );
  kinematics = new StarKinematics();
  _primary->AddGenerator( kinematics );
}
// ----------------------------------------------------------------------------
void lambda( Int_t nevents = 10, Int_t rngSeed = 20261109,
             TString mode = "gun", TString sdt = "sdt20240622",
             TString geomtag = "y2024" )
{
  gMode = mode;
  if ( gMode != "gun" && gMode != "ps" ) {
    printf( "lambda.C: mode must be \"gun\" or \"ps\", not \"%s\"\n", gMode.Data() );
    return;
  }
  gFzd = TString::Format( "fwd_lambda_%s.fzd", gMode.Data() );

  printf( "=== lambda.C ===\n" );
  printf( "  mode     %s\n", gMode.Data() );
  printf( "  events   %d   (%d Lambda per event)\n", nevents, nLambda );
  printf( "  seed     %d\n", rngSeed );
  printf( "  geometry %s  %s\n", geomtag.Data(), sdt.Data() );
  printf( "  pt       %.2f - %.2f GeV/c\n", minPt, maxPt );
  printf( "  eta      %.2f - %.2f\n", minEta, maxEta );
  printf( "  vertex   (%.2f, %.2f, %.2f) sigma (%.4f, %.4f, %.4f) cm\n",
          vtxX, vtxY, vtxZ, vtxSigmaX, vtxSigmaY, vtxSigmaZ );
  printf( "  output   %s\n", gFzd.Data() );
  printf( "\n  REMINDER: sim.C must be run with the SAME %s, and with an\n",
          sdt.Data() );
  printf( "  fGeom cache built for the same tag and lane, or the tracker\n" );
  printf( "  applies different geometry than GEANT used.\n\n" );

  hM   = new TH1F( "hM",   ";M_{p#pi};counts [10MeV]", 50, 1.0, 1.5 );
  hEta = new TH1F( "hEta", ";#eta_{#Lambda};counts",   50, 2.0, 4.5 );
  hP   = new TH1F( "hP",   ";p_{#Lambda} [GeV/c];counts", 50, 0, 50 );
  hL   = new TH1F( "hL",   ";#beta#gammac#tau [cm];counts", 50, 0, 400 );

  gSystem->Load( "libStarRoot.so" );

  gROOT->ProcessLine(".L bfc.C");
  {
    // The sdt MUST be explicit and MUST match reco. Without it the generator
    // silently uses today's date as the DB timestamp.
    TString simple = sdt + " " + geomtag + " geant gstar usexgeom agml ";
    bfc( 0, simple );
  }

  gSystem->Load( "libVMC.so" );
  gSystem->Load( "StarGeneratorUtil.so" );
  gSystem->Load( "StarGeneratorEvent.so" );
  gSystem->Load( "StarGeneratorBase.so" );
  gSystem->Load( "libMathMore.so" );
  gSystem->Load( "xgeometry.so" );

  StarRandom::seed( rngSeed );
  StarRandom::capture();

  _primary = new StarPrimaryMaker();
  {
    _primary -> SetFileName( TString::Format("fwd_lambda_%s.root", gMode.Data()) );
    chain -> AddBefore( "geant", _primary );
  }

  Kinematics();

  _primary -> Init();
  _primary -> SetSigma( vtxSigmaX, vtxSigmaY, vtxSigmaZ );
  _primary -> SetVertex( vtxX, vtxY, vtxZ );

  // DETP GEOM ... field=-5.0: the fz header records the field, reversed full,
  // as in the 2022+ data.
  geometry( geomtag );
  command( "gkine -4 0" );
  command( TString::Format( "gfile o %s", gFzd.Data() ) );

  // Decays are ON -- that is what makes the "gun" mode work. Do NOT add
  // "DCAY 0" here: it would stop the Lambda from ever decaying.

  trig( nevents );

  TFile *f = new TFile( "fwd_lambda_gen.root", "RECREATE" );
  f->cd();
  hM->Write();
  hEta->Write();
  hP->Write();
  hL->Write();
  f->Write();
  // Deliberately NOT f->Close(). StarPrimaryMaker writes its "stats" object at
  // Finish, into whatever gDirectory is current; closing this file leaves that
  // as Rint and the write fails with a TROOT::WriteTObject error. Leaving the
  // file open costs nothing -- Write() has already flushed it.

  printf( "\n=== lambda.C done: %d events -> %s ===\n", nevents, gFzd.Data() );
  printf( "  <p_Lambda>  = %.2f GeV/c\n", hP->GetMean() );
  printf( "  <bg c tau>  = %.2f cm  (mean decay length)\n", hL->GetMean() );

  command( "call agexit" );  // make sure STARSIM exits properly
}
// ----------------------------------------------------------------------------
