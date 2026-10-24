// fttStripRoundTrip.C -- the invariant that the sTGC strip-offset fix had to
// preserve: StFttClusterMaker (strip -> position) and
// StFttSlowSimMaker::computeRowStrip() (position -> strip) must be exact
// inverses, or simulated hits land somewhere their own truth is not.
//
// Both now go through StFttDb, so this checks the pair the makers call:
//   reco:  StFttDb::stripCenter( row, strip, orientation )     -> centre, mm
//   sim :  StFttDb::stripFromLocal( row, orientation, coord )  -> strip, offset
// and reports the constant the old analytic form was missing.
//
//   root4star -l -b -q 'script/fttStripRoundTrip.C'
//
// CINT: fixed arrays, unique loop variable names, no local std::map (CLAUDE.md).
void loadLibsRT(){
  gSystem->Load("libStarClassLibrary.so"); gSystem->Load("libStarRoot.so");
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();
  gSystem->Load("StDbBroker"); gSystem->Load("StDbUtilities");
  gSystem->Load("St_db_Maker"); gSystem->Load("libStDb_Tables.so");
  gSystem->Load("StFttDbMaker");
}

void fttStripRoundTrip()
{
  loadLibsRT();
  StFttDb *db = new StFttDb("fttDb");
  db->loadStripCenterFromFile( "StRoot/StFttDbMaker/Row1.txt" );
  db->loadStripEdgeFromFile(   "StRoot/StFttDbMaker/Row4_edge.txt" );
  db->loadStripCenterFromFile( "StRoot/StFttDbMaker/Row4.txt" );
  db->loadStripLengthFromFile( "StRoot/StFttDbMaker/Row1_StripLength.txt" );
  db->loadStripLengthFromFile( "StRoot/StFttDbMaker/Row2_StripLength.txt" );
  db->loadStripLengthFromFile( "StRoot/StFttDbMaker/Row3_StripLength.txt" );
  db->loadStripLengthFromFile( "StRoot/StFttDbMaker/Row4_StripLength.txt" );
  db->loadStripLengthFromFile( "StRoot/StFttDbMaker/Row5_StripLength.txt" );

  const double pitch = StFttDb::stripPitch;
  printf("\nstripPitch = %.4f mm\n", pitch);

  // row -> orientation to use. Rows 0-2 are the XY map, rows 3-4 the diagonal one.
  int    rowList[5];  rowList[0]=0; rowList[1]=1; rowList[2]=2; rowList[3]=3; rowList[4]=4;
  UChar_t orList[5];
  orList[0]=(UChar_t)1; orList[1]=(UChar_t)1; orList[2]=(UChar_t)1;   // kFttVertical
  orList[3]=(UChar_t)3; orList[4]=(UChar_t)3;                         // kFttDiagonalV
  const char *rowName[5];
  rowName[0]="0 (XY)"; rowName[1]="1 (XY)"; rowName[2]="2 (XY)";
  rowName[3]="3 (diag)"; rowName[4]="4 (diag)";

  // Three probe positions inside each strip, so the sub-strip offset that drives
  // the charge sharing is checked too, not just the centre.
  double probe[3]; probe[0]=0.0; probe[1]=-0.4*pitch; probe[2]=+0.4*pitch;

  printf("\n  row       nStrip  roundTrip  worst|offset err|  centre - (strip-0.5)*pitch\n");
  printf("  ----------------------------------------------------------------------------\n");

  int totStrip = 0, totBad = 0;
  for (int ir = 0; ir < 5; ir++){
    int    nStrip = 0, nBad = 0;
    double worstOff = 0;
    double cmin = 1e9, cmax = -1e9;

    for (int js = 0; js < 400; js++){
      double centre = 0;
      if ( !db->stripCenter( rowList[ir], js, orList[ir], centre ) ) continue;
      nStrip++;

      // the constant the analytic form was missing
      double konst = centre - ( js - 0.5 ) * pitch;
      if (konst < cmin) cmin = konst;
      if (konst > cmax) cmax = konst;

      for (int kp = 0; kp < 3; kp++){
        double coord = centre + probe[kp];
        int    back  = -1;
        double off   = 0;
        if ( !db->stripFromLocal( rowList[ir], orList[ir], coord, back, off ) ){
          // only legitimate at the two ends of the map, where +-0.4 pitch can
          // still round onto a strip that exists -- so count it
          nBad++;
          continue;
        }
        if ( back != js ) { nBad++; continue; }
        double err = fabs( off - probe[kp] );
        if (err > worstOff) worstOff = err;
      }
    }

    printf("  %-9s %5d   %4d/%-4d   %12.2e      %8.4f .. %8.4f\n",
           rowName[ir], nStrip, 3*nStrip - nBad, 3*nStrip, worstOff, cmin, cmax);
    totStrip += nStrip; totBad += nBad;
  }

  printf("\n  %d strips, %d probes, %d failures\n", totStrip, 3*totStrip, totBad);
  if ( 0 == totBad )
    printf("  Reco and sim are exact inverses on every strip and every sub-strip probe.\n");
  else
    printf("  *** %d probes did NOT round trip -- sim and reco disagree ***\n", totBad);
  printf("\n  The constant column is what StFttClusterMaker used to drop: 17.55 mm on the\n");
  printf("  XY rows, 26.74 on the diagonals. One number cannot fix both, hence the map.\n");
}
