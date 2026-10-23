// fttStripCenterCheck.C -- the measured strip centre (scMapXY/scMapDiag, what the
// reconstruction puts on a raw hit) against the analytic (strip-0.5)*pitch that
// StFttSlowSimMaker::computeRowStrip() inverts. One quantity, two owners.
//   root4star -l -b -q 'fttStripCenterCheck.C'
void loadLibsSC(){
  gSystem->Load("libStarClassLibrary.so"); gSystem->Load("libStarRoot.so");
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();
  gSystem->Load("StDbBroker"); gSystem->Load("StDbUtilities");
  gSystem->Load("St_db_Maker"); gSystem->Load("libStDb_Tables.so");
  gSystem->Load("StFttDbMaker");
}
void fttStripCenterCheck(){
  loadLibsSC();
  StFttDb *db = new StFttDb("fttDb");
  db->loadStripCenterFromFile("StRoot/StFttDbMaker/Row1.txt");
  db->loadStripLengthFromFile("StRoot/StFttDbMaker/Row1_StripLength.txt");
  double pitch = StFttDb::stripPitch;
  printf("\nstripPitch = %.4f mm\n", pitch);
  printf("\n strip   measured (stripSegment)   analytic (strip-0.5)*pitch   difference\n");
  double dmin = 1e9, dmax = -1e9;
  int nok = 0;
  for (int isc = 0; isc < 170; isc++){
    double xa,ya,xb,yb;
    if ( !db->stripSegment( 0, isc, (UChar_t)1 /*kFttVertical*/, xa, ya, xb, yb ) ) continue;
    double meas = xa;                       // vertical strip: x is the precision coord
    double ana  = ( isc - 0.5 ) * pitch;
    double d    = meas - ana;
    if (d < dmin) dmin = d;
    if (d > dmax) dmax = d;
    nok++;
    if (isc < 4 || isc > 164)
      printf("  %4d %18.4f %24.4f %15.4f\n", isc, meas, ana, d);
  }
  printf("\n  %d strips compared\n", nok);
  printf("  difference: min %.4f  max %.4f  spread %.4g mm\n", dmin, dmax, dmax-dmin);
  printf("  => %s offset of %.4f mm = %.3f strip pitches\n",
         (dmax-dmin) < 1e-3 ? "CONSTANT" : "NON-CONSTANT", dmin, dmin/pitch);
}
