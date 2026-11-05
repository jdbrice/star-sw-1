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
//
// GEOMETRY: same arguments and the same reason as script/genElectron45.C -- the
// old hardcoded y2024 uses FSTM config FSTMv1 -> FstmGeo.xml, which has ZERO
// <Misalign> tags and no halfid logic, so every sample generated before
// 2026-10-24 sat on IDEAL FST. Only the M variants carry the survey, and
// dev2022sm is what the Run-22 fGeom reco caches are built from. The "misalign"
// bfc option is what makes AgPosition read the tables at all, and sdt picks
// which tables. See genElectron45.C's header for the full variant map.
void genElectron2pi(int nevents = 4000, int seed = 1, int numPerEvent = 1,
                    const char* tag = "dev2022sm", const char* sdt = "20211112.000004",
                    int misalign = 1) {
    gROOT->LoadMacro("StRoot/StFwdTrackMaker/macro/sim/gen.C");

    nameParticle = "e-";
    numParticles = numPerEvent;
    minPt = 1.499; maxPt = 1.501; // StarRandom::flat() needs a strictly nonzero range
    minEta = 2.5;  maxEta = 4.0;
    minPhi = 0.0;  maxPhi = 2.0 * TMath::Pi();

    vtxX = 0.0; vtxY = 0.0; vtxZ = 0.0;
    vtxSigmaX = 0.0001; vtxSigmaY = 0.0001; vtxSigmaZ = 0.0001;

    // Name carries the geometry tag: an ideal-FST sample and a misaligned one
    // under one filename is the worst outcome. Renames the old ele2pi_np*.vz0.*
    // (which were y2024, i.e. ideal FST).
    fzdFilename = Form("ele2pi_np%d.%s.vz0.run%d.fzd", numPerEvent, tag, seed);
    primaryName = Form("ele2pi_np%d.%s.vz0.run%d.root", numPerEvent, tag, seed);

    printf("genElectron2pi: nevents=%d seed=%d numPerEvent=%d -> phi=[0,2pi], file=%s\n"
           "genElectron2pi: tag=%s sdt=%s misalign=%d\n",
           nevents, seed, numPerEvent, fzdFilename.Data(), tag, sdt, misalign);

    // ---- inlined gen(nevents, seed), minus the macro-path-breaking line ----
    gSystem->Load("libStarRoot.so");
    gROOT->LoadMacro("bfc.C");
    {
        TString simple = Form("sdt%s %s geant gstar usexgeom agml %s",
                              sdt, tag, misalign ? "misalign " : "");
        printf("genElectron2pi: bfc chain = [%s]\n", simple.Data());
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

    geometry(tag);   // gen.C helper, DETP GEOM ... field=-5.0: fz header records the reversed field
    command("gkine -4 0");
    command(TString::Format("gfile o %s", fzdFilename.Data()));

    trig(nevents);

    command("call agexit");
}
