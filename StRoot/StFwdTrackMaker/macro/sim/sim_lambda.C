//usr/bin/env root4star -l -b -q $0'('$1', '$2')'; exit $?
//
// sim_lambda.C -- DEPRECATED, and kept only so commands already in circulation
//                 keep working. Nothing in it was ever Lambda-specific.
//
// It existed because sim.C's geometry knobs were file-scope globals, so the
// caller had to do LoadMacro (which does NOT execute) -> set the globals ->
// call sim(), in that order. Both knobs are now trailing ARGUMENTS of sim(),
// so there is nothing left to wrap. Use sim.C directly:
//
//   root4star -l -b -q 'StRoot/StFwdTrackMaker/macro/sim/sim.C("fwd_lambda_gun.fzd",10,"sdt20211112.000004","fGeom.y2022c.misalign0.sdt20211112.000004.root")'
//
// Mind the argument ORDER when you switch: sim.C takes (file, n, ...) while
// this took (n, file, ...). That is reason enough to move rather than keep
// both spellings alive.
//
// Still true however you call it, and the whole point of the knobs: the DB lane
// and the geometry cache MUST match what the sample was GENERATED with. The
// generator defaults to y2022c (macro/sim/gen.C, genGeomTag) deliberately --
// that is the tag PR28 wired to the wedge mounting StFstDb::fstmConfig()
// reports -- so pair it with an fGeom.y2022c.* cache, never the fGeom.root
// symlink, which points at Run-22 DATA geometry.
//
// Delete this file once no one is passing the old argument order.

void sim_lambda( int n = 10,
                 const char *fzd  = "fwd_lambda_gun.fzd",
                 const char *sdt  = "sdt20240622",
                 const char *geom = "fGeom.y2024.misalign0.sdt20240622.root" )
{
    printf("sim_lambda.C is DEPRECATED -- sim.C now takes the lane and the cache\n"
           "   directly, in the order (file, n, sdt, geomCache):\n"
           "   sim.C(\"%s\",%d,\"%s\",\"%s\")\n"
           "   Forwarding this call.\n", fzd, n, sdt, geom);

    gROOT->LoadMacro("StRoot/StFwdTrackMaker/macro/sim/sim.C");
    gROOT->ProcessLine( Form("sim(\"%s\", %d, \"%s\", \"%s\");", fzd, n, sdt, geom) );
}
