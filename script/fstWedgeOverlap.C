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
int   gFwoOcc[4608];   // [disk-1][m-1][strip] flattened: (d*12 + m)*128 + s
long  gFwoCoin[3][2][2][4];   // edge hits : [disk-1][parity of m][my edge 0=s0-3 1=s124-127][y]
long  gFwoCtrl[3][2][4];   // control     : strips 60-67, same bookkeeping
long  gFwoEdge[3][2][2];      // edge hits seen, same split
long  gFwoCtrlN[3][2];     // control hits seen

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
   // the inner sensor covers r 5.0-13.625 and the outer 13.625-28.0, so one
   // track crosses a wedge at ONE radius and gives an inner hit OR an outer
   // hit, never both. Same-wedge inner+outer pairs are therefore different
   // particles and carry no phi relation. Kept because the flatness is worth
   // seeing rather than rediscovering: measured 15.7% in the sum window against
   // 11.5% in the diff window, i.e. nothing.
   //   To actually test the inner-vs-outer sense you need two DISKS: for a
   // vertex-pointing track r scales with z (r3/r1 = 178.781/151.75 = 1.178), so
   // a track with r1 in 11.57-13.625 is INNER on disk 1 and OUTER on disk 3 at
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
   }

   for (int ij = 0; ij < 3; ij++) for (int ik = 0; ik < 2; ik++) {
      gFwoCtrlN[ij][ik] = 0;
      for (int il = 0; il < 4; il++) gFwoCtrl[ij][ik][il] = 0;
      for (int im = 0; im < 2; im++) {
         gFwoEdge[ij][ik][im] = 0;
         for (int in = 0; in < 4; in++) gFwoCoin[ij][ik][im][in] = 0;
      }
   }

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

      for (int iz = 0; iz < 4608; iz++) gFwoOcc[iz] = 0;
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
         if (sp == 0 && dp >= 1 && dp <= 3 && wp >= 1 && wp <= 36) {
            int mi = 1 + (wp - 1) % 12;
            int ci = (int)(fp + 0.5);
            if (ci >= 0 && ci <= 127) {
               int kk = (dp - 1) * 12 + (mi - 1);
               gFwoInnN[kk]++;
               gFwoInnS[kk] = ci;
            }
         }
         if (sp < 1 || sp > 2) continue;                 // outer only
         if (dp < 1 || dp > 3 || wp < 1 || wp > 36) continue;
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
   printf("EXPECTED FLAT, and not a mapping failure: inner covers r 5.0-13.625 and\n"
          "outer 13.625-28.0, so a track gives one or the other, never both, and\n"
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
   hFine->Write(); hPS->Write(); hRS->Write(); hNH->Write();
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
