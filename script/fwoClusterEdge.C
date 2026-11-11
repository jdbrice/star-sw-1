// fwoClusterEdge.C -- cluster shape against distance from the wedge edge.
//
// WHY THIS EXISTS
//   The per-strip occupancy says the edge strip (d = 0) carries MORE hits than
//   the bulk -- outer +12%, inner +4% -- while its coincidence fraction with
//   the overlapping neighbour is only about half the plateau. The proposed
//   explanation was edge-truncated clusters: a track just outside the silicon
//   deposits on the last strip only, so it is reconstructed at d = 0 with no
//   partner possible. That picture makes a sharp prediction -- those clusters
//   should be BOTH narrower and lower in charge than bulk ones -- and
//   fstWedgeOverlap.C fills hClPhi/hClQ (outer) and hClPhiI/hClQI (inner) to
//   test it. Run this on the hadd-ed campaign output.
//
//   The answer on 60.4M events (fwo_20261110b) was neither of the two outcomes
//   that were anticipated:
//     * d = 0 carries the SAME charge as the bulk (ratio 0.999), so the
//       truncation picture is wrong -- a truncated cluster must lose charge.
//     * d = 0 is nonetheless 3% narrower, and that is pure geometry: a cluster
//       centred on the last strip cannot extend outward, so n >= 4 at d = 0 is
//       0.0001 against 0.003 in the bulk. Clipping of the SIZE without loss of
//       charge is what a one-sided window does.
//     * the surprise is at d = 1, 2, 3: a 20-80x excess of wide (n >= 5)
//       clusters -- 1.5%, 6.1%, 3.9% over a 0.08% bulk -- with the charge
//       excess to match (+9%, +21%, +18%). It is confined to three strips just
//       inside the wedge edge and is gone by d = 5. That is inside the overlap
//       region the width measurement uses, so it has to be understood before
//       the width is final.
//
// Usage:
//   root4star -l -b -q 'script/fwoClusterEdge.C("<hadd-ed fwo file>")'

void fwoClusterEdgeOne(TFile *f, const char *hname, const char *what)
{
   TH2F *h = (TH2F*)f->Get(hname);
   if (!h) { printf("  %s MISSING\n", hname); return; }
   printf("\n%s  --  %s\n", hname, what);
   printf("   d    entries      mean       rms   mean/bulk\n");

   // bulk = d 6..11, far enough from the edge that neither the clipping nor
   // the wide-cluster excess reaches it.
   double bulk = 0; int nb = 0;
   for (int d = 6; d <= 11; d++) {
      int ix = h->GetXaxis()->FindBin((double)d);
      TH1D *pb = h->ProjectionY("fwoCEpb", ix, ix);
      if (pb->GetEntries() > 0) { bulk += pb->GetMean(); nb++; }
      delete pb;
   }
   if (nb) bulk /= nb;

   for (int d = 0; d <= 11; d++) {
      int ix = h->GetXaxis()->FindBin((double)d);
      TH1D *p1 = h->ProjectionY("fwoCEp1", ix, ix);
      if (p1->GetEntries() < 1) { delete p1; continue; }
      printf("  %2d %10.0f %9.3f %9.3f %9.4f\n",
             d, p1->GetEntries(), p1->GetMean(), p1->GetRMS(),
             bulk > 0 ? p1->GetMean()/bulk : 0.0);
      delete p1;
   }
   printf("   bulk (d = 6..11) mean = %.3f\n", bulk);
}

void fwoClusterEdgeSpec(TFile *f, const char *hname, const char *what)
{
   TH2F *h = (TH2F*)f->Get(hname);
   if (!h) { printf("  %s MISSING\n", hname); return; }
   printf("\n%s size spectrum  --  %s  (fraction of that d's clusters)\n", hname, what);
   printf("   d      n=1     n=2     n=3     n=4    n>=5\n");
   for (int d = 0; d <= 11; d++) {
      int ix = h->GetXaxis()->FindBin((double)d);
      TH1D *ps = h->ProjectionY("fwoCEps", ix, ix);
      double tot = ps->Integral(1, 8);
      if (tot <= 0) { delete ps; continue; }
      printf("  %2d  %7.4f %7.4f %7.4f %7.4f %7.4f\n", d,
             ps->GetBinContent(1)/tot, ps->GetBinContent(2)/tot,
             ps->GetBinContent(3)/tot, ps->GetBinContent(4)/tot,
             ps->Integral(5, 8)/tot);
      delete ps;
   }
}

void fwoClusterEdge(const char *fn = "fwo_pairR.root")
{
   TFile *f = TFile::Open(fn);
   if (!f || f->IsZombie()) { printf("cannot open %s\n", fn); return; }

   printf("CLUSTER SHAPE vs DISTANCE FROM THE WEDGE EDGE\n");
   printf("The truncation hypothesis predicts d = 0 clusters are both narrower\n");
   printf("AND lower in charge. Read the charge ratio first: it is the half\n");
   printf("that cannot be faked by the one-sided window at the last strip.\n");

   fwoClusterEdgeOne(f, "hClPhi",  "OUTER cluster phi size");
   fwoClusterEdgeOne(f, "hClQ",    "OUTER cluster charge");
   fwoClusterEdgeOne(f, "hClPhiI", "INNER cluster phi size");
   fwoClusterEdgeOne(f, "hClQI",   "INNER cluster charge");

   fwoClusterEdgeSpec(f, "hClPhi",  "OUTER");
   fwoClusterEdgeSpec(f, "hClPhiI", "INNER");
}
