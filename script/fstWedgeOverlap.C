// fstWedgeOverlap.C
//
// Does wedge W's outer silicon really overlap wedge W+-1's, and which strips of
// the neighbour does it meet? Answered from data alone, with NO tracking, NO
// geometry and NO DB: every quantity used is a hardware index out of
// StMuFstHit (disk, wedge, sensor, meanPhiStrip). Nothing here depends on the
// XML, on fstmConfig, or on StFstDb, so no geometry setting can bias the
// answer -- which is the point: the result is what the detector says, to be
// compared against what our mapping predicts.
//
// THE GEOMETRY BEING TESTED (numbers inlined; StFstConsts.h consts are often
// not visible to CINT)
//   12 wedges per disk, 3 disks (36 wedges, kFstNumWedges=36 -- not 4 disks).
//   Strip pitch in phi = 30/128 deg = 0.234375 deg, the same for all sensors.
//   Inner sensor (sensor 0, rStrip 0-3): 128 strips, phiStrip 0-127, spans
//     exactly 30 deg -> abuts its neighbours, no overlap. NOT used here.
//   Outer sensors (sensor 1 rStrip 4-7 phiStrip 0-63, sensor 2 phiStrip
//     64-127): 128 strips at the same pitch = 30 deg of silicon PLUS the 1 deg
//     dead gap on the wedge centreline -> 31 deg of silhouette for a 30 deg
//     wedge pitch. AGML carries exactly this: the two outer FTUS span
//     344.5-359.5 and 0.5-15.5 in the wedge frame.
//   So each wedge overhangs its nominal 30 deg boundary by 0.5 deg on BOTH
//   sides, and since both neighbours overhang towards each other the shared
//   band is 1.0 deg wide and BOTH wedges have silicon across all of it:
//   1.0 / 0.234375 = 4.27 strips from each wedge. Hence strips 0-3 and 124-127.
//
// WHAT THE HISTOGRAMS ARE
//   One TH2F per wedge (36) plus summaries. For every outer hit in wedge W:
//     x = its phiStrip, 0-127 (all 128, so the non-edge strips are the
//         accidental-coincidence control in the same plot)
//     y = which neighbour-edge also had a hit THIS EVENT, 4 bins:
//         0 : wedge W-1 strips 0-3        1 : wedge W-1 strips 124-127
//         2 : wedge W+1 strips 0-3        3 : wedge W+1 strips 124-127
//   W+-1 is index arithmetic within the disk, mod 12. NOTE THE DIRECTION:
//   wedge index increases as phi DECREASES (wedge 1 spans 60-90 deg, wedge 2
//   30-60, wedge 3 0-30 ...), so W+1 is the LOWER-phi neighbour. Getting this
//   backwards inverts the conclusion, so the axis labels say it too.
//
// WHAT THE TWO HYPOTHESES PREDICT
//   Our mapping (StFstDb::stripPhi, which production uses and which does not
//   depend on fstmConfig) puts the SAME-numbered edges of adjacent wedges
//   against each other, because kFstzDirct = {1,-1,1,-1,...} alternates the
//   strip direction wedge to wedge. Evaluating it for the edge strips gives
//   pairings like w1 s0 <-> w12 s3 and w1 s127 <-> w2 s124, all within 0.063
//   deg -- the two strip grids in the overlap band are interleaved, offset by
//   1.0 - 4*0.234375 = 0.0625 deg. The predicted cell therefore ALTERNATES
//   with wedge parity, and disk 2 is reversed again (kFstzFilp = {1,-1,1}):
//
//     disks 1,3   odd  m : strips 0-3 -> y=0,  strips 124-127 -> y=3
//                 even m : strips 0-3 -> y=2,  strips 124-127 -> y=1
//     disk  2     odd  m : strips 0-3 -> y=2,  strips 124-127 -> y=1
//                 even m : strips 0-3 -> y=0,  strips 124-127 -> y=3
//
//   If the wedges are physically mirrored relative to that convention, every
//   wedge's strip order reverses, so my strips 0-3 meet the OTHER neighbour:
//   the mirrored hypothesis is exactly y0<->y2 and y1<->y3. That is a
//   two-cell-vs-two-cell choice with nothing in between, which is why this
//   test settles the wedge flip without tracking.
//
//   Because the predicted cell alternates with parity, the all-wedge average
//   histogram necessarily shows BOTH cells at ~50/50 and says nothing on its
//   own. The parity-split summaries are the ones to read.
//
// The input name MUST contain ".lis", ".files" or ".MuDst.root", or end in "/":
// StMuChainMaker::make keys off that and silently refuses anything else, so
// run24_mudst.txt has to be handed over as a copy called e.g. zf_run24.list.
//
// Usage: root4star -l -b -q 'fstWedgeOverlap.C("zf_run24.list",20,-1,"out.root")'

void loadLibsFWO();

int   gFwoInnN[36];        // inner hits per wedge this event
int   gFwoInnS[36];        // the inner strip, when there is exactly one
int   gFwoOutN[36];        // outer hits per wedge this event
int   gFwoOutS[36];
int   gFwoOccI[4608];      // inner-sensor occupancy, same indexing
long  gFwoCoinI[3][2][2][4];  // inner edge coincidences [disk][parity][end][y]
long  gFwoEdgeI[3][2][2];
long  gFwoProfI[12][2];       // inner: strip from edge, [dist][0=pred 1=mirror]
long  gFwoProfIR[4][12][2];   // inner, split by rStrip 0-3: [r][dist][pred/mirror].
                              // THE SHARPER OF THE TWO r TESTS. The inner has NO 1 deg
                              // gap, so its width is 128*eps - 128*d/r with no 4.27-strip
                              // pedestal: the three hypotheses give 0.99 (pure pitch),
                              // 1.92 (the outer fit) and 3.50 (pure shift) strips at
                              // rStrip 0, a factor 3.5 apart, where across the whole
                              // outer range they differ by 5%. Small absolute width is
                              // not the same as small discriminating power.
long  gFwoProfR[4][12][2];    // outer, split by rStrip 4-7: [r][dist][pred/mirror].
                              // A radial displacement of the strip fan scales the
                              // apparent pitch by (1-d/R), so the overlap width would
                              // shrink with r (5.51 -> 5.11 strips across rStrip 4-7
                              // for d = -0.176 cm); a pitch error gives the same width
                              // in all four rows. This separates them without needing
                              // the inner sensor at all.
int   gFwoOcc[4608];   // [disk-1][m-1][strip] flattened: (d*12 + m)*128 + s
long  gFwoCoin[3][2][2][4];   // edge hits : [disk-1][parity of m][my edge 0=s0-3 1=s124-127][y]
long  gFwoCtrl[3][2][4];   // control     : strips 60-67, same bookkeeping
long  gFwoEdge[3][2][2];      // edge hits seen, same split
long  gFwoCtrlN[3][2];     // control hits seen

int fwoAnyI(int disk, int m, int lo, int hi)
{
   int base = ((disk - 1) * 12 + (m - 1)) * 128;
   for (int jb = lo; jb <= hi; jb++) if (gFwoOccI[base + jb]) return 1;
   return 0;
}

int fwoAny(int disk, int m, int lo, int hi)
{
   // any hit in strips [lo,hi] of (disk, m)?  disk 1-3, m 1-12
   int base = ((disk - 1) * 12 + (m - 1)) * 128;
   for (int ja = lo; ja <= hi; ja++) if (gFwoOcc[base + ja]) return 1;
   return 0;
}

void fstWedgeOverlap(const Char_t *fileList = "run24_mudst.txt",
                     int nFiles = 20, int nEvents = -1,
                     const Char_t *outFile = "fstWedgeOverlap.root")
{
   loadLibsFWO();

   StChain *chain = new StChain("StChain");
   StMuDstMaker *muDstMaker = new StMuDstMaker(0, 0, "", fileList, "MuDst.root", nFiles);
   TChain *muChain = muDstMaker->chain();
   if (!muChain) {
      printf("\nERROR: no chain. StMuChainMaker::make only accepts an input name\n"
             "       containing \".lis\", \".files\" or \".MuDst.root\", or ending in \"/\".\n"
             "       Got \"%s\" -- copy the list to something like zf_run24.list.\n", fileList);
      return;
   }
   long nAvail = muChain->GetEntries();
   printf("\nfiles requested %d, events in chain %ld\n", nFiles, nAvail);

   if (chain->Init()) { printf("ERROR: chain Init failed\n"); return; }

   // ---- books ----
   TH2F *hw[3][12];
   for (int ia = 0; ia < 3; ia++) {
      for (int ib = 0; ib < 12; ib++) {
         hw[ia][ib] = new TH2F(Form("hD%dW%02d", ia + 1, ib + 1),
                               Form("disk %d wedge %d (m=%d, %s);outer phiStrip;neighbour edge hit",
                                    ia + 1, ib + 1, ib + 1, (ib % 2) ? "even m" : "odd m"),
                               128, 0, 128, 4, 0, 4);
      }
   }
   TH2F *hAll  = new TH2F("hAll",  "all 36 wedges (parity washes out, see hOdd/hEven);outer phiStrip;neighbour edge hit", 128, 0, 128, 4, 0, 4);
   TH2F *hOdd  = new TH2F("hOdd",  "odd m (1,3,5,7,9,11), all disks;outer phiStrip;neighbour edge hit",  128, 0, 128, 4, 0, 4);
   TH2F *hEven = new TH2F("hEven", "even m (2,4,6,8,10,12), all disks;outer phiStrip;neighbour edge hit", 128, 0, 128, 4, 0, 4);
   TH2F *hDisk[3];
   for (int ic = 0; ic < 3; ic++)
      hDisk[ic] = new TH2F(Form("hDisk%d", ic + 1), Form("disk %d, all wedges;outer phiStrip;neighbour edge hit", ic + 1),
                           128, 0, 128, 4, 0, 4);

   // front/back: the mounting sense, orientation = kFstzFilp[disk] * kFstzDirct[m],
   // i.e. odd m on disks 1,3 and even m on disk 2 are one sense, the rest the other.
   // This is the grouping the pattern actually alternates with, so it is the plot
   // to read; hOdd/hEven only coincide with it on disks 1 and 3.
   TH2F *hFront = new TH2F("hFront", "front wedges (orientation +1: odd m disks 1,3 / even m disk 2);outer phiStrip;neighbour edge hit", 128, 0, 128, 4, 0, 4);
   TH2F *hBack  = new TH2F("hBack",  "back wedges (orientation -1);outer phiStrip;neighbour edge hit", 128, 0, 128, 4, 0, 4);

   // Inner vs outer in the same wedge. THIS COMES OUT FLAT, and it has to:
   // the inner sensor covers r 5.0-16.5 (rStrip 0-3) and the outer 16.5-28.0 (rStrip 4-7), so one
   // track crosses a wedge at ONE radius and gives an inner hit OR an outer
   // hit, never both. Same-wedge inner+outer pairs are therefore different
   // particles and carry no phi relation. Kept because the flatness is worth
   // seeing rather than rediscovering: measured 15.7% in the sum window against
   // 11.5% in the diff window, i.e. nothing.
   //   To actually test the inner-vs-outer sense you need two DISKS: for a
   // vertex-pointing track r scales with z (r3/r1 = 178.781/151.75 = 1.178), so
   // a track with r1 in 14.01-16.5 is INNER on disk 1 and OUTER on disk 3 at
   // the same phi. Pair those (and expect sum = const for the 1-3 pair, which
   // share sgn, but difference = const whenever disk 2 is involved, since
   // kFstzFilp inverts it). Narrow band, so it needs much more data.
   //   Note that the relative inner/outer sense is NOT in dispute anyway: the
   // channel-map note (note_fst_channel_map_20260928.txt, from the cosmic
   // channel map scan) and stripPhi agree that they run opposite ways. What the
   // note cannot fix, and what the overlap above measures, is the ABSOLUTE
   // sense -- a bench sheet marked "face up"/"face down" does not say which way
   // the module points once installed.
   TH2F *hIO   = new TH2F("hIO", "same wedge, same event: inner vs outer strip;inner phiStrip;outer phiStrip", 128, 0, 128, 128, 0, 128);
   TH1F *hIOs  = new TH1F("hIOsum",  "inner + outer strip (flat if uncorrelated, peak if opposite senses);s_inner + s_outer;pairs", 256, 0, 256);
   TH1F *hIOd  = new TH1F("hIOdiff", "inner - outer strip (peak if same sense);s_inner - s_outer;pairs", 257, -128.5, 128.5);

   // does the DB mapping index match the note's APV<->phi assignment?
   // note: outer APV5 phi 0-31, APV6 32-63, APV7 64-95, APV8 96-127;
   //       inner APV1 0-31, APV2 32-63, APV3 64-95, APV4 96-127.
   TH2F *hAPV = new TH2F("hAPVvsPhiStrip", "APV vs phiStrip (DB map check);phiStrip;getApv()", 128, 0, 128, 18, -1, 17);

   // strip-to-strip in the overlap band: does my s0 really meet the neighbour's
   // s3?  x = my strip (0-3 then 124-127), y = neighbour strip, same 8 slots
   TH2F *hFine = new TH2F("hFine", "edge strip vs neighbour edge strip, all wedges;my edge strip;neighbour edge strip", 8, 0, 8, 8, 0, 8);
   // the same for the INNER sensor. The inner has no 1 deg gap, so on nominal
   // geometry there is nothing to see here at all -- any structure is the
   // overlap the pitch excess produces. Same 8 slots as the outer so the two
   // plots can be read side by side; the inner band is narrower, so expect the
   // correlation to sit in the corner rather than fill the block.
   TH2F *hFineI = new TH2F("hFineI", "INNER edge strip vs neighbour edge strip, all wedges;my edge strip;neighbour edge strip", 8, 0, 8, 8, 0, 8);
   // PAIR MAP -- the centroid method, and the sharper of the two.
   //
   // Fold both ends to a distance from the wedge edge (d = c for the low end,
   // 127-c for the high end -- akio's 2*127 - i - j) and histogram my distance
   // against the partner's. The two strip grids run OPPOSITE ways through the
   // overlap band, so every pair at the same physical angle lies on
   //     d_me + d_nb = const,
   // and that constant is the overlap width less one: W = <s> + 1. The SUM is
   // therefore a direct measure of the width, and its centroid is far sharper
   // than an edge read off a falling yield profile -- the same reason the
   // pairing centroid of section 5 beat the yield method there.
   //
   // 12 slots a side, not the 8 of hFine: with W ~ 5.3 my strip 0 pairs with
   // the neighbour's strip ~4, and an 8-slot window (0-3) loses that cell
   // entirely, sampling the s = 4 ridge at only three of its five cells.
   //
   // Strips 0 and 127 are visibly less efficient than their neighbours, which
   // biases <s> downward, and a PAIR carries the suppression TWICE (mine and
   // the partner's). hPhiStripPerSensor is the per-strip occupancy to correct
   // with, on the assumption that illumination is uniform across the wedge.
   TH2F *hPair  = new TH2F("hPair",  "OUTER pair map;my distance from the wedge edge;partner distance", 12, -0.5, 11.5, 12, -0.5, 11.5);
   TH2F *hPairI = new TH2F("hPairI", "INNER pair map;my distance from the wedge edge;partner distance", 12, -0.5, 11.5, 12, -0.5, 11.5);
   // EDGE-TRUNCATED CLUSTERS, the test of why strip 0 misbehaves.
   //
   // The per-strip occupancy says strip 0 has MORE hits than the bulk (outer
   // +12%, inner +4%), not fewer, while its coincidence fraction is about half
   // the plateau. The natural reading is that the excess is clusters whose
   // charge is cut off by the silicon edge: the reconstructed meanPhiStrip
   // piles up on the last strip while the true position is outside the sensor,
   // so there is no partner to find.
   //
   // If that is right, hits at d=0 should be NARROWER in phi and LOWER in
   // charge than the bulk, because part of the cluster is missing. If instead
   // they look exactly like bulk clusters, the excess is something else and
   // cluster shape gives no handle -- which matters, because with most FST
   // hits single-strip there may be no cut with any separation to offer.
   // These two plots answer that directly rather than by assertion.
   TH2F *hClPhi  = new TH2F("hClPhi",  "OUTER cluster phi size vs distance from edge;distance from the wedge edge;nRawHitsPhi", 12, -0.5, 11.5, 8, 0.5, 8.5);
   TH2F *hClPhiI = new TH2F("hClPhiI", "INNER cluster phi size vs distance from edge;distance from the wedge edge;nRawHitsPhi", 12, -0.5, 11.5, 8, 0.5, 8.5);
   TH2F *hClQ    = new TH2F("hClQ",    "OUTER cluster charge vs distance from edge;distance from the wedge edge;charge", 12, -0.5, 11.5, 100, 0, 4000);
   TH2F *hClQI   = new TH2F("hClQI",   "INNER cluster charge vs distance from edge;distance from the wedge edge;charge", 12, -0.5, 11.5, 100, 0, 4000);

   // diagnostics that must be checked before believing any of it
   TH2F *hPS = new TH2F("hPhiStripPerSensor", "phiStrip by sensor;meanPhiStrip;sensor", 132, -2, 130, 4, -1, 3);
   TH2F *hRS = new TH2F("hRStripPerSensor",   "rStrip by sensor;meanRStrip;sensor",      12, -2,  10, 4, -1, 3);
   TH1F *hNH = new TH1F("hNHits", "FST hits per event;hits;events", 100, 0, 100);

   const Char_t *yl[4] = {"W-1 (higher phi) s0-3", "W-1 (higher phi) s124-127",
                          "W+1 (lower phi) s0-3",  "W+1 (lower phi) s124-127"};
   for (int id = 0; id < 4; id++) {
      hAll->GetYaxis()->SetBinLabel(id + 1, yl[id]);
      hOdd->GetYaxis()->SetBinLabel(id + 1, yl[id]);
      hEven->GetYaxis()->SetBinLabel(id + 1, yl[id]);
      hFront->GetYaxis()->SetBinLabel(id + 1, yl[id]);
      hBack->GetYaxis()->SetBinLabel(id + 1, yl[id]);
      for (int ie = 0; ie < 3; ie++) hDisk[ie]->GetYaxis()->SetBinLabel(id + 1, yl[id]);
      for (int ig = 0; ig < 3; ig++) for (int ih = 0; ih < 12; ih++)
         hw[ig][ih]->GetYaxis()->SetBinLabel(id + 1, yl[id]);
   }
   const Char_t *el[8] = {"0", "1", "2", "3", "124", "125", "126", "127"};
   for (int ii = 0; ii < 8; ii++) {
      hFine->GetXaxis()->SetBinLabel(ii + 1, el[ii]);
      hFine->GetYaxis()->SetBinLabel(ii + 1, el[ii]);
      hFineI->GetXaxis()->SetBinLabel(ii + 1, el[ii]);
      hFineI->GetYaxis()->SetBinLabel(ii + 1, el[ii]);
   }

   for (int ij = 0; ij < 3; ij++) for (int ik = 0; ik < 2; ik++) {
      gFwoCtrlN[ij][ik] = 0;
      for (int il = 0; il < 4; il++) gFwoCtrl[ij][ik][il] = 0;
      for (int im = 0; im < 2; im++) {
         gFwoEdge[ij][ik][im] = 0; gFwoEdgeI[ij][ik][im] = 0;
         for (int in = 0; in < 4; in++) { gFwoCoin[ij][ik][im][in] = 0; gFwoCoinI[ij][ik][im][in] = 0; }
      }
   }
   for (int io2 = 0; io2 < 12; io2++) { gFwoProfI[io2][0] = 0; gFwoProfI[io2][1] = 0; }
   for (int io5 = 0; io5 < 4; io5++)
      for (int io6 = 0; io6 < 12; io6++) { gFwoProfIR[io5][io6][0] = 0; gFwoProfIR[io5][io6][1] = 0; }
   for (int ip3 = 0; ip3 < 4; ip3++) for (int iq3 = 0; iq3 < 12; iq3++) { gFwoProfR[ip3][iq3][0] = 0; gFwoProfR[ip3][iq3][1] = 0; }

   long nEvSeen = 0, nEvWithFst = 0, nHitTot = 0, nHitOuter = 0, nRawTot = 0, nIOPair = 0;
   long nStop = (nEvents > 0 && nEvents < nAvail) ? nEvents : nAvail;

   // ---- event loop ----
   for (long iev = 0; iev < nStop; iev++) {
      chain->Clear();
      if (kStOK != chain->Make()) { printf("Make() non-OK at event %ld, stopping\n", iev); break; }
      nEvSeen++;

      StMuFstCollection *fc = StMuDst::muFstCollection();
      if (!fc) continue;
      int nh = fc->numberOfHits();
      nRawTot += fc->numberOfRawHits();
      hNH->Fill(nh);
      if (nh <= 0) continue;
      nEvWithFst++;
      nHitTot += nh;

      for (int iz = 0; iz < 4608; iz++) { gFwoOcc[iz] = 0; gFwoOccI[iz] = 0; }
      for (int iz2 = 0; iz2 < 36; iz2++) { gFwoInnN[iz2] = 0; gFwoOutN[iz2] = 0; gFwoInnS[iz2] = -1; gFwoOutS[iz2] = -1; }

      // pass 1: occupancy of the outer sensors
      for (int ip = 0; ip < nh; ip++) {
         StMuFstHit *hp = fc->getHit(ip);
         if (!hp) continue;
         int dp = hp->getDisk();
         int wp = hp->getWedge();
         int sp = hp->getSensor();
         float fp = hp->getMeanPhiStrip();
         hPS->Fill(fp, sp);
         hRS->Fill(hp->getMeanRStrip(), sp);
         hAPV->Fill(fp, hp->getApv());
         {   // cluster shape against distance from the wedge edge, both ends folded
            int cc = (int)(fp + 0.5);
            int dd = (cc <= 11) ? cc : ((cc >= 116) ? 127 - cc : -1);
            if (dd >= 0) {
               if (sp == 0) { hClPhiI->Fill(dd, hp->getNRawHitsPhi()); hClQI->Fill(dd, hp->getCharge()); }
               else         { hClPhi ->Fill(dd, hp->getNRawHitsPhi()); hClQ ->Fill(dd, hp->getCharge()); }
            }
         }
         if (sp == 0 && dp >= 1 && dp <= 3 && wp >= 1 && wp <= 36) {
            int mi = 1 + (wp - 1) % 12;
            int ci = (int)(fp + 0.5);
            if (ci >= 0 && ci <= 127) {
               int kk = (dp - 1) * 12 + (mi - 1);
               gFwoInnN[kk]++;
               gFwoInnS[kk] = ci;
            }
         }
         if (dp < 1 || dp > 3 || wp < 1 || wp > 36) continue;
         if (sp == 0) {                                  // inner occupancy
            int mi2 = 1 + (wp - 1) % 12;
            int ci2 = (int)(fp + 0.5);
            if (ci2 >= 0 && ci2 <= 127)
               gFwoOccI[((dp - 1) * 12 + (mi2 - 1)) * 128 + ci2] = 1;
         }
         if (sp < 1 || sp > 2) continue;                 // outer only below
         int mp = 1 + (wp - 1) % 12;
         int cp = (int)(fp + 0.5);
         if (cp < 0 || cp > 127) continue;
         gFwoOcc[((dp - 1) * 12 + (mp - 1)) * 128 + cp] = 1;
         gFwoOutN[(dp - 1) * 12 + (mp - 1)]++;
         gFwoOutS[(dp - 1) * 12 + (mp - 1)] = cp;
         nHitOuter++;
      }

      // pass 2: for each outer hit, which neighbour edge fired
      for (int iq = 0; iq < nh; iq++) {
         StMuFstHit *hq = fc->getHit(iq);
         if (!hq) continue;
         int dq = hq->getDisk();
         int wq = hq->getWedge();
         int sq = hq->getSensor();
         if (sq < 1 || sq > 2) continue;
         if (dq < 1 || dq > 3 || wq < 1 || wq > 36) continue;
         int mq = 1 + (wq - 1) % 12;
         int cq = (int)(hq->getMeanPhiStrip() + 0.5);
         if (cq < 0 || cq > 127) continue;

         int mPrev = 1 + (mq - 1 + 11) % 12;   // W-1 : HIGHER phi
         int mNext = 1 + (mq - 1 +  1) % 12;   // W+1 : LOWER  phi

         int yy[4];
         yy[0] = fwoAny(dq, mPrev,   0,   3);
         yy[1] = fwoAny(dq, mPrev, 124, 127);
         yy[2] = fwoAny(dq, mNext,   0,   3);
         yy[3] = fwoAny(dq, mNext, 124, 127);

         int par = (mq % 2) ? 0 : 1;           // 0 = odd m, 1 = even m
         {  // r-binned width, outer only
            int rsq = (int)(hq->getMeanRStrip() + 0.5) - 4;
            int dq3 = (cq <= 11) ? cq : ((cq >= 116) ? 127 - cq : -1);
            int en3 = (cq <= 11) ? 0 : 1;
            int fl3 = ((par == 1) ? 1 : 0) ^ ((dq == 2) ? 1 : 0);
            int pr3 = (en3 == 0) ? (fl3 ? 2 : 0) : (fl3 ? 1 : 3);
            int mr3 = (en3 == 0) ? (fl3 ? 0 : 2) : (fl3 ? 3 : 1);
            if (rsq >= 0 && rsq < 4 && dq3 >= 0 && dq3 < 12) {
               if (yy[pr3]) gFwoProfR[rsq][dq3][0]++;
               if (yy[mr3]) gFwoProfR[rsq][dq3][1]++;
            }
         }
         int isEdge = (cq <= 3 || cq >= 124);
         int myEnd  = (cq <= 3) ? 0 : 1;        // 0 = strips 0-3, 1 = strips 124-127
         int isCtrl = (cq >= 60 && cq <= 67);
         if (isEdge) gFwoEdge[dq - 1][par][myEnd]++;
         if (isCtrl) gFwoCtrlN[dq - 1][par]++;

         for (int ir = 0; ir < 4; ir++) {
            if (!yy[ir]) continue;
            hw[dq - 1][mq - 1]->Fill(cq, ir);
            hAll->Fill(cq, ir);
            hDisk[dq - 1]->Fill(cq, ir);
            if (par == 0) hOdd->Fill(cq, ir); else hEven->Fill(cq, ir);
            int fb = ((par == 1) ? 1 : 0) ^ ((dq == 2) ? 1 : 0);   // 0 = front, 1 = back
            if (fb == 0) hFront->Fill(cq, ir); else hBack->Fill(cq, ir);
            if (isEdge) gFwoCoin[dq - 1][par][myEnd][ir]++;
            if (isCtrl) gFwoCtrl[dq - 1][par][ir]++;
         }

         // the wide pair map, same partner convention as hFine
         {
            int dMe = (cq <= 11) ? cq : ((cq >= 116) ? 127 - cq : -1);
            if (dMe >= 0) {
               int nbrP = (cq <= 11) ? mPrev : mNext;
               for (int dNb = 0; dNb < 12; dNb++) {
                  int nstrip = (cq <= 11) ? dNb : (127 - dNb);
                  if (gFwoOcc[((dq - 1) * 12 + (nbrP - 1)) * 128 + nstrip])
                     hPair->Fill(dMe, dNb);
               }
            }
         }
         // fine strip-to-strip, only for my edge strips and the neighbour edge
         // that our mapping says should be the partner
         if (isEdge) {
            int myslot = (cq <= 3) ? cq : (cq - 120);         // 0-3 or 4-7
            for (int is = 0; is < 8; is++) {
               int nstrip = (is <= 3) ? is : (is + 120);
               int nbr = (cq <= 3) ? mPrev : mNext;           // partner under our mapping
               if (gFwoOcc[((dq - 1) * 12 + (nbr - 1)) * 128 + nstrip])
                  hFine->Fill(myslot, is);
            }
         }
      }

      // INNER edge coincidences. The inner sensor spans exactly 128 pitches with
      // no gap, so on nominal geometry it has NO overlap at all and these must be
      // accidental. If the strip PITCH is larger than nominal the inner silicon
      // spills past 30 deg and does overlap; if instead the 1 deg GAP is wider,
      // the inner is untouched. That is what separates the two, and the outer
      // measurement cannot: both move the same combination 127+(gap-30)/pitch.
      //   The inner runs the opposite way to the outer inside a wedge, so its
      // predicted cell is the outer's with the two ends swapped.
      for (int iq2 = 0; iq2 < nh; iq2++) {
         StMuFstHit *hq2 = fc->getHit(iq2);
         if (!hq2) continue;
         if (hq2->getSensor() != 0) continue;
         int dq2 = hq2->getDisk(), wq2 = hq2->getWedge();
         if (dq2 < 1 || dq2 > 3 || wq2 < 1 || wq2 > 36) continue;
         int mq2 = 1 + (wq2 - 1) % 12;
         int cq2 = (int)(hq2->getMeanPhiStrip() + 0.5);
         if (cq2 < 0 || cq2 > 127) continue;
         int mP2 = 1 + (mq2 - 1 + 11) % 12;
         int mN2 = 1 + (mq2 - 1 +  1) % 12;
         int yy2[4];
         yy2[0] = fwoAnyI(dq2, mP2,   0,   1);
         yy2[1] = fwoAnyI(dq2, mP2, 126, 127);
         yy2[2] = fwoAnyI(dq2, mN2,   0,   1);
         yy2[3] = fwoAnyI(dq2, mN2, 126, 127);
         int par2  = (mq2 % 2) ? 0 : 1;
         int flip2 = ((par2 == 1) ? 1 : 0) ^ ((dq2 == 2) ? 1 : 0);
         int dist  = (cq2 <= 11) ? cq2 : ((cq2 >= 116) ? 127 - cq2 : -1);
         int end2  = (cq2 <= 11) ? 0 : 1;
         // inner s0 end pairs with the neighbour the OUTER s124-127 end pairs with
         int pred2 = (end2 == 0) ? (flip2 ? 0 : 2) : (flip2 ? 3 : 1);
         int mirr2 = (end2 == 0) ? (flip2 ? 2 : 0) : (flip2 ? 1 : 3);
         if (dist >= 0 && dist < 12) {
            if (yy2[pred2]) gFwoProfI[dist][0]++;
            if (yy2[mirr2]) gFwoProfI[dist][1]++;
            // same quantity split by radius. rStrip runs 0-3 across the inner
            // sensor (r 5.0-16.5 cm), so no offset here, unlike the outer's -4.
            int rsi = (int)(hq2->getMeanRStrip() + 0.5);
            if (rsi >= 0 && rsi < 4) {
               if (yy2[pred2]) gFwoProfIR[rsi][dist][0]++;
               if (yy2[mirr2]) gFwoProfIR[rsi][dist][1]++;
            }
         }
         {
            int dMeI = (cq2 <= 11) ? cq2 : ((cq2 >= 116) ? 127 - cq2 : -1);
            if (dMeI >= 0) {
               int nbrPI = (cq2 <= 11) ? mP2 : mN2;
               for (int dNbI = 0; dNbI < 12; dNbI++) {
                  int nstripI = (cq2 <= 11) ? dNbI : (127 - dNbI);
                  if (gFwoOccI[((dq2 - 1) * 12 + (nbrPI - 1)) * 128 + nstripI])
                     hPairI->Fill(dMeI, dNbI);
               }
            }
         }
         // fine strip-to-strip for the inner, same construction as hFine
         if (cq2 <= 3 || cq2 >= 124) {
            int myslotI = (cq2 <= 3) ? cq2 : (cq2 - 120);
            for (int isI = 0; isI < 8; isI++) {
               int nstripI = (isI <= 3) ? isI : (isI + 120);
               int nbrI = (cq2 <= 3) ? mP2 : mN2;
               if (gFwoOccI[((dq2 - 1) * 12 + (nbrI - 1)) * 128 + nstripI])
                  hFineI->Fill(myslotI, isI);
            }
         }
         if (cq2 <= 1 || cq2 >= 126) {
            gFwoEdgeI[dq2 - 1][par2][end2]++;
            for (int ir2 = 0; ir2 < 4; ir2++)
               if (yy2[ir2]) gFwoCoinI[dq2 - 1][par2][end2][ir2]++;
         }
      }

      // inner vs outer in the same wedge: require exactly one of each, so no
      // combinatorial smearing
      for (int iw2 = 0; iw2 < 36; iw2++) {
         if (gFwoInnN[iw2] != 1 || gFwoOutN[iw2] != 1) continue;
         hIO->Fill(gFwoInnS[iw2], gFwoOutS[iw2]);
         hIOs->Fill(gFwoInnS[iw2] + gFwoOutS[iw2]);
         hIOd->Fill(gFwoInnS[iw2] - gFwoOutS[iw2]);
         nIOPair++;
      }

      if (nEvSeen % 2000 == 0) printf("  ... %ld events, %ld FST hits\n", nEvSeen, nHitTot);
   }

   // ---- report ----
   printf("\nevents read %ld, with FST hits %ld, FST hits %ld (outer %ld), raw hits %ld\n",
          nEvSeen, nEvWithFst, nHitTot, nHitOuter, nRawTot);
   if (nHitTot > 0)
      printf("mean FST hits/event %.2f\n", (double)nHitTot / (double)(nEvSeen ? nEvSeen : 1));

   printf("\nphiStrip range per sensor (must be 0-127 / 0-63 / 64-127 or the\n"
          "strip numbering assumed here is wrong):\n");
   for (int it = 0; it < 3; it++) {
      int lo = -1, hi = -1;
      for (int iu = 1; iu <= hPS->GetNbinsX(); iu++) {
         if (hPS->GetBinContent(iu, hPS->GetYaxis()->FindBin(it)) > 0) {
            if (lo < 0) lo = (int)hPS->GetXaxis()->GetBinCenter(iu);
            hi = (int)hPS->GetXaxis()->GetBinCenter(iu);
         }
      }
      printf("   sensor %d : phiStrip %d to %d\n", it, lo, hi);
   }

   printf("\nEDGE HITS: coincidences by cell, split by which end the hit is on.\n");
   printf("P marks the cell our mapping predicts, M the mirrored alternative.\n");
   printf("                                  y0:W-1 s0-3  y1:W-1 s124-7  y2:W+1 s0-3  y3:W+1 s124-7   edge hits\n");
   long sumP = 0, sumM = 0;
   for (int iv = 0; iv < 3; iv++) {
      for (int iw = 0; iw < 2; iw++) {
         for (int ix2 = 0; ix2 < 2; ix2++) {
            // our mapping: disks 1,3 odd m -> s0-3 meets W-1 s0-3 (y0) and
            // s124-127 meets W+1 s124-127 (y3); even m and disk 2 each reverse it
            int flip = ((iw == 1) ? 1 : 0) ^ ((iv == 1) ? 1 : 0);
            int pred = (ix2 == 0) ? (flip ? 2 : 0) : (flip ? 1 : 3);
            int mirr = (ix2 == 0) ? (flip ? 0 : 2) : (flip ? 3 : 1);
            printf("   disk %d %-6s s%-8s", iv + 1, (iw == 0) ? "odd m" : "even m",
                   (ix2 == 0) ? "0-3" : "124-127");
            for (int iy2 = 0; iy2 < 4; iy2++)
               printf(" %10ld%s", gFwoCoin[iv][iw][ix2][iy2],
                      (iy2 == pred) ? "P" : ((iy2 == mirr) ? "M" : " "));
            printf("  %9ld\n", gFwoEdge[iv][iw][ix2]);
            sumP += gFwoCoin[iv][iw][ix2][pred];
            sumM += gFwoCoin[iv][iw][ix2][mirr];
         }
      }
   }
   printf("\n   TOTAL in the predicted cells  P = %ld\n", sumP);
   printf("   TOTAL in the mirrored cells   M = %ld\n", sumM);
   printf("\nCONTROL, strips 60-67 (no overlap there; this is the accidental level)\n");
   for (int ix = 0; ix < 3; ix++) {
      for (int iy = 0; iy < 2; iy++) {
         printf("   disk %d %-7s      %10ld   %10ld   %10ld   %10ld   %9ld\n",
                ix + 1, (iy == 0) ? "odd m" : "even m",
                gFwoCtrl[ix][iy][0], gFwoCtrl[ix][iy][1], gFwoCtrl[ix][iy][2], gFwoCtrl[ix][iy][3],
                gFwoCtrlN[ix][iy]);
      }
   }
   printf("\nOUTER WIDTH BY rStrip -- radial shift or pitch?\n");
   printf("A strip fan converging d off the beam axis scales the apparent pitch by\n");
   printf("(1-d/R), so the width would fall with r: 5.51, 5.34, 5.21, 5.11 strips for\n");
   printf("d = -0.176 cm. A pitch error gives 5.25 in all four rows.\n");
   printf("   rStrip   <r>  :");
   for (int ir4 = 0; ir4 < 10; ir4++) printf(" %6d", ir4);
   printf("\n");
   double rmid[4] = {17.94, 20.81, 23.69, 26.56};
   for (int is4 = 0; is4 < 4; is4++) {
      printf("     %d    %5.2f  :", is4 + 4, rmid[is4]);
      for (int it4 = 0; it4 < 10; it4++) printf(" %6ld", gFwoProfR[is4][it4][0]);
      printf("\n              mirr :");
      for (int iu4 = 0; iu4 < 10; iu4++) printf(" %6ld", gFwoProfR[is4][iu4][1]);
      printf("\n");
   }

   printf("\nINNER EDGE COINCIDENCES -- does the inner sensor overlap too?\n");
   printf("Nominal geometry says NO: the inner is 128 pitches = exactly 30 deg.\n");
   printf("A pitch larger than nominal would make it overlap (x1.006 -> 0.76 strips);\n");
   printf("a wider 1 deg gap would not touch it. Edge window s0-1 and s126-127.\n");
   long sumPI = 0, sumMI = 0, sumAI = 0, nEdgeI = 0;
   for (int iv2 = 0; iv2 < 3; iv2++)
      for (int iw2 = 0; iw2 < 2; iw2++)
         for (int ix3 = 0; ix3 < 2; ix3++) {
            int fl3 = ((iw2 == 1) ? 1 : 0) ^ ((iv2 == 1) ? 1 : 0);
            int pr3 = (ix3 == 0) ? (fl3 ? 0 : 2) : (fl3 ? 3 : 1);
            int mr3 = (ix3 == 0) ? (fl3 ? 2 : 0) : (fl3 ? 1 : 3);
            sumPI += gFwoCoinI[iv2][iw2][ix3][pr3];
            sumMI += gFwoCoinI[iv2][iw2][ix3][mr3];
            for (int iy3 = 0; iy3 < 4; iy3++)
               if (iy3 != pr3 && iy3 != mr3) sumAI += gFwoCoinI[iv2][iw2][ix3][iy3];
            nEdgeI += gFwoEdgeI[iv2][iw2][ix3];
         }
   printf("   inner edge hits %ld\n", nEdgeI);
   printf("   predicted cells P = %ld\n", sumPI);
   printf("   mirrored  cells M = %ld\n", sumMI);
   printf("   neither (accidental, 2 cells) = %ld, i.e. %.1f per cell\n", sumAI, 0.5 * sumAI);
   printf("\n   inner coincidences vs distance from the wedge edge:\n");
   printf("   dist :");
   for (int iz3 = 0; iz3 < 12; iz3++) printf(" %6d", iz3);
   printf("\n   pred :");
   for (int ia3 = 0; ia3 < 12; ia3++) printf(" %6ld", gFwoProfI[ia3][0]);
   printf("\n   mirr :");
   for (int ib3 = 0; ib3 < 12; ib3++) printf(" %6ld", gFwoProfI[ib3][1]);
   printf("\n");

   printf("\nINNER WIDTH BY rStrip -- the sharper separation of pitch from radial shift.\n");
   printf("The inner has no 1 deg gap, so its width IS 128*eps - 128*d/r with no\n");
   printf("4.27-strip pedestal. Predictions at rStrip 0,1,2,3 (<r> 6.44, 9.31, 12.19,\n");
   printf("15.06 cm):  pure pitch (eps=+0.77%%)  0.99 0.99 0.99 0.99\n");
   printf("            outer-fit mixture        1.92 1.48 1.24 1.10\n");
   printf("            pure shift (d=-0.176cm)  3.50 2.42 1.85 1.50\n");
   printf("A factor 3.5 apart at rStrip 0, against 5%% across the whole outer range --\n");
   printf("so a coarse width here still separates them, PROVIDED the statistics hold\n");
   printf("up after a 4-way split. Check the pred row sums before believing a slope.\n");
   double rmidI[4] = {6.44, 9.31, 12.19, 15.06};
   printf("   rStrip   <r>  :");
   for (int ic3 = 0; ic3 < 10; ic3++) printf(" %6d", ic3);
   printf("     sum\n");
   for (int id3 = 0; id3 < 4; id3++) {
      long rowSum = 0;
      printf("     %d    %5.2f  :", id3, rmidI[id3]);
      for (int ie3 = 0; ie3 < 10; ie3++) { printf(" %6ld", gFwoProfIR[id3][ie3][0]); rowSum += gFwoProfIR[id3][ie3][0]; }
      printf("  %6ld\n", rowSum);
      printf("              mirr :");
      for (int if3 = 0; if3 < 10; if3++) printf(" %6ld", gFwoProfIR[id3][if3][1]);
      printf("\n");
   }

   printf("\nAPV vs phiStrip, the DB-map check. The channel-map note says outer\n"
          "APV5/6/7/8 carry phi 0-31/32-63/64-95/96-127 and inner APV1/2/3/4 the same\n"
          "ranges. If the data disagrees, the DB map is not the note's indexing and\n"
          "every statement below about strip 0 has to be re-read.\n");
   for (int ia2 = 0; ia2 < 17; ia2++) {
      int lo2 = -1, hi2 = -1;
      long na2 = 0;
      for (int ib2 = 1; ib2 <= hAPV->GetNbinsX(); ib2++) {
         double cc = hAPV->GetBinContent(ib2, hAPV->GetYaxis()->FindBin(ia2));
         if (cc > 0) { if (lo2 < 0) lo2 = (int)hAPV->GetXaxis()->GetBinCenter(ib2);
                       hi2 = (int)hAPV->GetXaxis()->GetBinCenter(ib2); na2 += (long)cc; }
      }
      if (na2 > 0) printf("   apv %2d : phiStrip %3d to %3d   (%ld hits)\n", ia2, lo2, hi2, na2);
   }

   printf("\nINNER vs OUTER in the same wedge, one of each, %ld pairs.\n", nIOPair);
   printf("EXPECTED FLAT, and not a mapping failure: inner covers r 5.0-16.5 and\n"
          "outer 16.5-28.0, so a track gives one or the other, never both, and\n"
          "these pairs are different particles. See the comment at the booking.\n");
   if (nIOPair > 0) {
      printf("   sum  : mean %7.2f  rms %6.2f  peak at %3.0f with %.0f pairs\n",
             hIOs->GetMean(), hIOs->GetRMS(),
             hIOs->GetBinCenter(hIOs->GetMaximumBin()), hIOs->GetMaximum());
      printf("   diff : mean %7.2f  rms %6.2f  peak at %3.0f with %.0f pairs\n",
             hIOd->GetMean(), hIOd->GetRMS(),
             hIOd->GetBinCenter(hIOd->GetMaximumBin()), hIOd->GetMaximum());
      double fs = 0, fd = 0;
      for (int ic2 = 1; ic2 <= hIOs->GetNbinsX(); ic2++) {
         double xc = hIOs->GetBinCenter(ic2);
         if (xc > 118 && xc < 136) fs += hIOs->GetBinContent(ic2);
      }
      for (int id2 = 1; id2 <= hIOd->GetNbinsX(); id2++) {
         double xd = hIOd->GetBinCenter(id2);
         if (xd > -9 && xd < 9) fd += hIOd->GetBinContent(id2);
      }
      printf("   pairs with sum in 119-135 : %.0f (%.1f%%)\n", fs, 100.0 * fs / nIOPair);
      printf("   pairs with |diff| < 9     : %.0f (%.1f%%)\n", fd, 100.0 * fd / nIOPair);
   }

   printf("\nour mapping predicts (stripPhi, the convention production uses):\n"
          "   disks 1,3  odd m : s0-3 -> y0, s124-127 -> y3 | even m : s0-3 -> y2, s124-127 -> y1\n"
          "   disk 2     odd m : s0-3 -> y2, s124-127 -> y1 | even m : s0-3 -> y0, s124-127 -> y3\n"
          "a physical mirror relative to that is exactly y0<->y2 and y1<->y3.\n");

   TFile *fout = new TFile(outFile, "recreate");
   for (int i1 = 0; i1 < 3; i1++) for (int i2 = 0; i2 < 12; i2++) hw[i1][i2]->Write();
   hAll->Write(); hOdd->Write(); hEven->Write(); hFront->Write(); hBack->Write();
   hIO->Write(); hIOs->Write(); hIOd->Write(); hAPV->Write();
   for (int i3 = 0; i3 < 3; i3++) hDisk[i3]->Write();
   hFine->Write(); hFineI->Write(); hPair->Write(); hPairI->Write(); hClPhi->Write(); hClPhiI->Write(); hClQ->Write(); hClQI->Write(); hPS->Write(); hRS->Write(); hNH->Write();

   // ---- the counters, as histograms, so hadd can sum them ----
   // Everything above this point is a TH2F and merges for free. The coincidence
   // counters are plain C arrays that only ever reached stdout, so a condor run
   // merged with hadd would keep the occupancy plots and silently lose the
   // numbers the run is FOR. Flatten them into TH1D, one bin per element, in
   // exactly the loop order written here -- script/fwoRead.C unflattens with the
   // same order and must be changed with this.
   TH1D *cProfR = new TH1D("cProfR", "gFwoProfR[r 0-3][dist 0-11][pred,mirr]", 96, 0, 96);
   int kbR = 0;
   for (int cw1 = 0; cw1 < 4; cw1++)
      for (int cw2 = 0; cw2 < 12; cw2++)
         for (int cw3 = 0; cw3 < 2; cw3++) cProfR->SetBinContent(++kbR, (double) gFwoProfR[cw1][cw2][cw3]);

   TH1D *cProfI = new TH1D("cProfI", "gFwoProfI[dist 0-11][pred,mirr]", 24, 0, 24);
   int kbI = 0;
   for (int cw4 = 0; cw4 < 12; cw4++)
      for (int cw5 = 0; cw5 < 2; cw5++) cProfI->SetBinContent(++kbI, (double) gFwoProfI[cw4][cw5]);

   TH1D *cCoin  = new TH1D("cCoin",  "gFwoCoin[disk][par][end][y]",  48, 0, 48);
   TH1D *cCoinI = new TH1D("cCoinI", "gFwoCoinI[disk][par][end][y]", 48, 0, 48);
   int kbC = 0;
   for (int cw6 = 0; cw6 < 3; cw6++)
      for (int cw7 = 0; cw7 < 2; cw7++)
         for (int cw8 = 0; cw8 < 2; cw8++)
            for (int cw9 = 0; cw9 < 4; cw9++) {
               ++kbC;
               cCoin ->SetBinContent(kbC, (double) gFwoCoin [cw6][cw7][cw8][cw9]);
               cCoinI->SetBinContent(kbC, (double) gFwoCoinI[cw6][cw7][cw8][cw9]);
            }

   TH1D *cEdge  = new TH1D("cEdge",  "gFwoEdge[disk][par][end]",  12, 0, 12);
   TH1D *cEdgeI = new TH1D("cEdgeI", "gFwoEdgeI[disk][par][end]", 12, 0, 12);
   int kbE = 0;
   for (int cx1 = 0; cx1 < 3; cx1++)
      for (int cx2 = 0; cx2 < 2; cx2++)
         for (int cx3 = 0; cx3 < 2; cx3++) {
            ++kbE;
            cEdge ->SetBinContent(kbE, (double) gFwoEdge [cx1][cx2][cx3]);
            cEdgeI->SetBinContent(kbE, (double) gFwoEdgeI[cx1][cx2][cx3]);
         }

   TH1D *cCtrl  = new TH1D("cCtrl",  "gFwoCtrl[disk][par][y]", 24, 0, 24);
   int kbT = 0;
   for (int cx4 = 0; cx4 < 3; cx4++)
      for (int cx5 = 0; cx5 < 2; cx5++)
         for (int cx6 = 0; cx6 < 4; cx6++) cCtrl->SetBinContent(++kbT, (double) gFwoCtrl[cx4][cx5][cx6]);

   TH1D *cCtrlN = new TH1D("cCtrlN", "gFwoCtrlN[disk][par]", 6, 0, 6);
   int kbN = 0;
   for (int cx7 = 0; cx7 < 3; cx7++)
      for (int cx8 = 0; cx8 < 2; cx8++) cCtrlN->SetBinContent(++kbN, (double) gFwoCtrlN[cx7][cx8]);

   // Event/hit totals, so a merged file still knows how much data it came from.
   TH1D *cTot = new TH1D("cTot", "0 nEvSeen 1 nEvWithFst 2 nHitTot 3 nHitOuter 4 nRawTot 5 nIOPair", 6, 0, 6);
   cTot->SetBinContent(1, (double) nEvSeen);
   cTot->SetBinContent(2, (double) nEvWithFst);
   cTot->SetBinContent(3, (double) nHitTot);
   cTot->SetBinContent(4, (double) nHitOuter);
   cTot->SetBinContent(5, (double) nRawTot);
   cTot->SetBinContent(6, (double) nIOPair);

   // gFwoProfIR must be persisted too or a condor run cannot be merged: the
   // printouts above go to stdout, and only what is packed into a TH1D here
   // survives hadd. Same packing order as cProfR.
   TH1D *cProfIR = new TH1D("cProfIR", "gFwoProfIR[r 0-3][dist 0-11][pred,mirr]", 96, 0, 96);
   int kbIR = 0;
   for (int cw6 = 0; cw6 < 4; cw6++)
      for (int cw7 = 0; cw7 < 12; cw7++)
         for (int cw8 = 0; cw8 < 2; cw8++) cProfIR->SetBinContent(++kbIR, (double) gFwoProfIR[cw6][cw7][cw8]);

   cProfR->Write(); cProfI->Write(); cProfIR->Write(); cCoin->Write(); cCoinI->Write();
   cEdge->Write();  cEdgeI->Write(); cCtrl->Write(); cCtrlN->Write(); cTot->Write();

   fout->Close();
   printf("\nwrote %s\n", outFile);
}

void loadLibsFWO()
{
   gSystem->Load("libStarClassLibrary.so");
   gSystem->Load("libStarRoot.so");
   gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
   loadSharedLibraries();
   gSystem->Load("St_base.so");
   gSystem->Load("StUtilities.so");
   gSystem->Load("libPhysics.so");
   gSystem->Load("StStrangeMuDstMaker");
   gSystem->Load("StMuDSTMaker");
   gSystem->Load("StStarLogger.so");
}
