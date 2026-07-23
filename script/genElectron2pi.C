// genElectron2pi.C
//
// Electron particle gun with phi released to the full 2pi range (unlike
// genElectron45.C's single narrow-phi quadrant gun), and a configurable
// number of electrons per event -- for testing matching under more
// realistic occupancy (numParticles>1 puts multiple tracks/disk/event, closer
// to real multi-track data than the single-track 45deg sample).
// Same gen.C-inlining rationale as genElectron45.C (its own gen() wrapper
// breaks the macro path for bfc.C).
//
// Usage:
//   root4star -b -q 'script/genElectron2pi.C(4000, 1, 1)'   // 1 e-/event
//   root4star -b -q 'script/genElectron2pi.C(4000, 1, 8)'   // 8 e-/event

void genElectron2pi(int nevents = 4000, int seed = 1, int numPerEvent = 1) {
    gROOT->LoadMacro("StRoot/StFwdTrackMaker/macro/sim/gen.C");

    nameParticle = "e-";
    numParticles = numPerEvent;
    minPt = 1.499; maxPt = 1.501; // StarRandom::flat() needs a strictly nonzero range
    minEta = 2.5;  maxEta = 4.0;
    minPhi = 0.0;  maxPhi = 2.0 * TMath::Pi();

    vtxX = 0.0; vtxY = 0.0; vtxZ = 0.0;
    vtxSigmaX = 0.0001; vtxSigmaY = 0.0001; vtxSigmaZ = 0.0001;

    fzdFilename = Form("ele2pi_np%d.vz0.run%d.fzd", numPerEvent, seed);
    primaryName = Form("ele2pi_np%d.vz0.run%d.root", numPerEvent, seed);

    printf("genElectron2pi: nevents=%d seed=%d numPerEvent=%d -> phi=[0,2pi], file=%s\n",
           nevents, seed, numPerEvent, fzdFilename.Data());

    // ---- inlined gen(nevents, seed), minus the macro-path-breaking line ----
    gSystem->Load("libStarRoot.so");
    gROOT->LoadMacro("bfc.C");
    {
        TString simple = "sdt20211016 y2024 geant gstar usexgeom agml ";
        bfc(0, simple);
    }

    gSystem->Load("libVMC.so");
    gSystem->Load("StarGeneratorUtil.so");
    gSystem->Load("StarGeneratorEvent.so");
    gSystem->Load("StarGeneratorBase.so");
    gSystem->Load("libMathMore.so");
    gSystem->Load("xgeometry.so");

    StarRandom::seed(seed);
    StarRandom::capture();

    _primary = new StarPrimaryMaker();
    _primary->SetFileName(primaryName);
    chain->AddBefore("geant", _primary);

    Kinematics();

    _primary->Init();
    _primary->SetSigma(vtxSigmaX, vtxSigmaY, vtxSigmaZ);
    _primary->SetVertex(vtxX, vtxY, vtxZ);

    command("gkine -4 0");
    command(TString::Format("gfile o %s", fzdFilename.Data()));

    trig(nevents);

    command("call agexit");
}
