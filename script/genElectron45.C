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
//   root4star -b -q 'script/genElectron45.C(500, 1, 95.0)'   // wedge-11 edge
//
// GEOMETRY. The tag and sdt are arguments, because the default used to be y2024
// and that is the WRONG geometry for anything about FST wedges: y2022*, y2023*,
// y2024* and y2025 all use FSTM config FSTMv1 -> FstmGeo.xml, which has ZERO
// <Misalign> tags and no halfid logic (measured: all 108 FTUS nodes sit at
// exactly (0,0)). Only the M variants carry the survey --
//     FSTMv1m  -> FstmGeoM.xml   tags dev2022m, dev2026
//     FSTMv1sm -> FstmGeoSM.xml  tag  dev2022sm
// -- and the iwedgem.ge.6 half-support-structure fix (d4378cbde9) was applied to
// those two files only. dev2022sm is also what every Run-22 fGeom reco cache is
// built from, so it is the default here: generated truth then sits on the same
// misaligned geometry the reconstruction reads back.
//
// The "misalign" bfc option is what makes AgPosition actually read the survey
// tables; without it the M variant builds ideal anyway. sdt selects WHICH tables
// (latest beginTime <= sdt wins), so it has to match the cache the reco side
// uses -- 20211112.000004 is the post-strip-frame-fix Run-22 entry.
//
// Pass tag="y2024", sdt="20240622", misalign=0 to reproduce the old ideal-FST
// behaviour.
void genElectron45(int nevents = 20, int seed = 1, float phiDeg = 45.0,
                   const char* tag = "dev2022sm", const char* sdt = "20211112.000004",
                   int misalign = 1) {
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

    // The output name carries BOTH phiDeg and the geometry tag. It used to be
    // hardcoded "ele45" whatever was asked for, so a second angle silently
    // overwrote the first -- and now that the tag is settable, an ideal-FST
    // sample and a misaligned one would collide the same way. Two different
    // geometries under one filename is the worst of the three outcomes, so the
    // tag goes in even though it renames the old ele45.vz0.* (which was y2024,
    // i.e. ideal FST, and should not be confused with a dev2022sm sample).
    fzdFilename = Form("ele%g.%s.vz0.run%d.fzd", phiDeg, tag, seed);
    primaryName = Form("ele%g.%s.vz0.run%d.root", phiDeg, tag, seed);

    printf("genElectron45: nevents=%d seed=%d phiDeg=%.1f -> phi=[%.5f,%.5f] rad, file=%s\n"
           "genElectron45: tag=%s sdt=%s misalign=%d\n",
           nevents, seed, phiDeg, minPhi, maxPhi, fzdFilename.Data(), tag, sdt, misalign);

    // ---- inlined gen(nevents, seed), minus the macro-path-breaking line ----
    gSystem->Load("libStarRoot.so");
    gROOT->LoadMacro("bfc.C");
    {
        TString simple = Form("sdt%s %s geant gstar usexgeom agml %s",
                              sdt, tag, misalign ? "misalign " : "");
        printf("genElectron45: bfc chain = [%s]\n", simple.Data());
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

    geometry(tag);   // gen.C helper, DETP GEOM <tag> field=-5.0: fz header records the reversed field
    command("gkine -4 0");
    command(TString::Format("gfile o %s", fzdFilename.Data()));

    trig(nevents);

    command("call agexit");
}
