// testGate.C -- exercise StFstRawHitMaker's occupancy gate on real AuAu DAQ.
//
//   root4star -b -q 'testGate.C(5,-1,"<daq>")'    gate OFF (default), reference
//   root4star -b -q 'testGate.C(5,1000,"<daq>")'  gate ON above 1000 FST raw hits
//
// bfc(0,...) initialises the chain and returns before the event loop (bfc.C
// line 226), which is where we reach in and set the cut -- the same place a
// production steering macro or a BFC .call attribute would.
//
// What to check in the log:
//   gate OFF -> no "OCCUPANCY GATE" lines, FST/FTT/FwdTrack all run
//   gate ON  -> "FST OCCUPANCY GATE: event dropped" per gated event, naming the
//               makers switched off; FCS still runs; events still written.

void testGate(int nev = 5, int maxRawHits = -1,
              const char *daq = "st_physics_24192017_raw_0700047.daq")
{
   const char *opts = "DbV20261006 P2023a StiCA fst fstMuRawHit ftt fcs FwdTrack BEmcChkStat";

   gROOT->LoadMacro("$STAR/StRoot/macros/bfc.C");
   bfc(0, opts, daq);                       // init only

   StMaker *top = StMaker::GetChain();
   if (!top) { printf("testGate: NO CHAIN\n"); return; }

   StMaker *mk = top->GetMakerInheritsFrom("StFstRawHitMaker");
   if (!mk) { printf("testGate: NO StFstRawHitMaker IN CHAIN\n"); return; }

   if (maxRawHits > 0) {
      // Set through the interpreter so this macro needs no compile-time
      // dependency on the maker's header.
      gROOT->ProcessLine(Form("((StFstRawHitMaker*)%p)->setMaxRawHits(%d);",
                              (void *) mk, maxRawHits));
      printf("testGate: maxRawHits set to %d on %s\n", maxRawHits, mk->GetName());
   } else {
      printf("testGate: gate left OFF (default)\n");
   }

   // EventLoop lives on StChain, not StMaker, and bfc.C forces SetRunNumber(-2)
   // so InitRun fires at least once. Drive bfc.C's own 'chain' global so the
   // types are resolved by the interpreter rather than here.
   gROOT->ProcessLine("{ StEvtHddr *hd = (StEvtHddr*)chain->GetDataSet(\"EvtHddr\");"
                      "  if (hd) hd->SetRunNumber(-2); }");
   gROOT->ProcessLine(Form("chain->EventLoop(1,%d,0);", nev));
   gROOT->ProcessLine("chain->Finish();");
}
