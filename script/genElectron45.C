// genElectron45.C
//
// Single-electron particle gun at a chosen phi, narrow enough to stay inside
// one FTT quadrant -- same idea as ~/fcstrk11/star-sw-fwd/jpsi/electron_45.html.
// Loads StFwdTrackMaker/macro/sim/gen.C for its helper functions (geometry,
// command, trig, Kinematics) and global tunables, but does NOT call its
// gen() wrapper directly: gen()'s own gROOT->SetMacroPath() call overwrites
// the macro path with a stale, AFS-only list that no longer contains
// wherever this STAR environment's bfc.C actually lives, breaking the very
// next line's ".L bfc.C" in-function. Inlines the rest of gen()'s body
// instead, using the correct (STAR-configured, untouched) macro path.
//
// Usage:
//   root4star -b -q 'script/genElectron45.C(20, 1, 45.0)'

void genElectron45(int nevents = 20, int seed = 1, float phiDeg = 45.0) {
    gROOT->LoadMacro("StRoot/StFwdTrackMaker/macro/sim/gen.C");

    nameParticle = "e-";
    numParticles = 1;
    minPt = 1.499; maxPt = 1.501; // StarRandom::flat() needs a strictly nonzero range
    minEta = 2.5;  maxEta = 4.0;

    float phiRad = phiDeg * TMath::Pi() / 180.0;
    minPhi = phiRad - 0.001;
    maxPhi = phiRad + 0.001;

    vtxX = 0.0; vtxY = 0.0; vtxZ = 0.0;
    vtxSigmaX = 0.0001; vtxSigmaY = 0.0001; vtxSigmaZ = 0.0001;

    fzdFilename = Form("ele45.vz0.run%d.fzd", seed);
    primaryName = Form("ele45.vz0.run%d.root", seed);

    printf("genElectron45: nevents=%d seed=%d phiDeg=%.1f -> phi=[%.5f,%.5f] rad, file=%s\n",
           nevents, seed, phiDeg, minPhi, maxPhi, fzdFilename.Data());

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
