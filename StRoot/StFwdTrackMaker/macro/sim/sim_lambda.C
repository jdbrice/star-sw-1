//usr/bin/env root4star -l -b -q $0'('$1', '$2')'; exit $?
//
// sim_lambda.C -- run the forward chain over a Lambda sample made by lambda.C,
//                 with the geometry lane and the geometry cache set to match
//                 the generator. Run it from the TOP of the checkout.
//
//   root4star -l -b -q 'StRoot/StFwdTrackMaker/macro/sim/sim_lambda.C(10)'
//
// This exists because sim.C's two geometry knobs are file-scope globals, which
// means the sequence has to be: LoadMacro (which does NOT execute), set the
// globals, and only then call sim(). Getting that order wrong leaves the
// defaults in place, and the defaults are the Run-22 lane and the fGeom.root
// symlink -- so a Run-24 sample reconstructs against Run-22 geometry without a
// single complaint in the log. See the _geomCache comment in sim.C.
//
// The RAVE vertex count per event comes out of the log as
//   FwdVertexFinding: <n> RAVE vertices from <m> global tracks [id=.. z=.. nTrk=..]
// one line per event, from FwdTracker::findFwdVertices.

void sim_lambda( int n = 10,
                 const char *fzd  = "fwd_lambda_gun.fzd",
                 const char *sdt  = "sdt20240622",
                 const char *geom = "fGeom.y2024.misalign0.sdt20240622.root" )
{
    gROOT->LoadMacro("StRoot/StFwdTrackMaker/macro/sim/sim.C");

    // Set sim.C's globals through the interpreter, so they resolve against
    // sim.C's own scope rather than this macro's.
    gROOT->ProcessLine( Form("_sdt = \"%s\";", sdt) );
    gROOT->ProcessLine( Form("_geomCache = \"%s\";", geom) );

    printf("sim_lambda: %s, %d events, lane %s, cache %s\n", fzd, n, sdt, geom);

    gROOT->ProcessLine( Form("sim(\"%s\", %d);", fzd, n) );
}
