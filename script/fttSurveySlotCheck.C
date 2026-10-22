// fttSurveySlotCheck.C -- does the DB survey slot reproduce StFttDb's hardcoded numbers?
//
//   root4star -l -b -q 'fttSurveySlotCheck.C(20211012)'   // the survey slot (-11..-14)
//   root4star -l -b -q 'fttSurveySlotCheck.C(20211010)'    // the ideal slot, for contrast
//
// The hardcoded X_shift_Quad*/Y_shift_Quad*/idealPlaneZLocations_Quad* are Run 22
// survey of the built detector. They are now also a DB entry at 2021-10-11 00:00:00,
// written as (hardcoded origin - AGML placeholder) so that the DB path composes back
// to the same numbers. This macro checks that, end to end, through StFttDb itself:
// the hardcoded values are read out of the object first, then the DB tables are
// loaded over them, and the two are differenced. Nothing is inlined here, so the
// test cannot drift from the code it checks.
//
// Reads ./StarDb only, so it tests the macro-directory files before any push.
// CINT: fixed-size arrays, unique loop variables (see CLAUDE.md).

void loadLibsFttSlot()
{
  gSystem->Load("libStarClassLibrary.so");
  gSystem->Load("libStarRoot.so");
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();
  gSystem->Load("StarMagField");
  gSystem->Load("StMagF");
  gSystem->Load("StDetectorDbMaker");
  gSystem->Load("StTpcDb");
  gSystem->Load("StDbBroker");
  gSystem->Load("StDbUtilities");
  gSystem->Load("St_db_Maker");
  gSystem->Load("libStDb_Tables.so");
  gSystem->Load("StFttDbMaker");
}

void fttSurveySlotCheck(int date = 20211012, int time = 0)
{
  loadLibsFttSlot();

  // files only, and only ours: $STAR/StarDb ships Geometry/stgc too
  St_db_Maker *dbMk = new St_db_Maker("db", "StarDb");
  dbMk->SetFlavor("ofl");
  dbMk->SetDateTime(date, time);
  dbMk->Init();
  dbMk->Make();

  StFttDb *fttDb = new StFttDb("fttDb");

  // the hardcoded numbers, straight out of the object
  double hcX[4][4], hcY[4][4], hcZ[4][4];
  fttDb->resetGeometryToHardcoded();
  for (int iqh = 0; iqh < 4; iqh++) {
    for (int iph = 0; iph < 4; iph++) {
      hcX[iqh][iph] = fttDb->xShift(iqh, iph);
      hcY[iqh][iph] = fttDb->yShift(iqh, iph);
      hcZ[iqh][iph] = fttDb->zLoc (iqh, iph);
    }
  }

  TDataSet *DB = dbMk->GetInputDB("Geometry/stgc");
  if (!DB) { printf("\nfttSurveySlotCheck: no Geometry/stgc at %d %06d\n", date, time); return; }

  const char *nm[3];
  nm[0] = "stgcOnTpc"; nm[1] = "stationOnStgc"; nm[2] = "pentOnStation";
  St_Survey *tb[3];
  printf("\nfttSurveySlotCheck: Geometry/stgc from ./StarDb at %d %06d\n", date, time);
  for (int itb = 0; itb < 3; itb++) {
    tb[itb] = (St_Survey *) DB->Find(nm[itb]);
    if (!tb[itb]) { printf("  %-15s NOT FOUND\n", nm[itb]); continue; }
    TDatime val[2];
    dbMk->GetValidity((TTable *)tb[itb], val);
    printf("  %-15s %2d rows   begin %s\n",
           nm[itb], (int)tb[itb]->GetNRows(), val[0].AsSQLString());
  }
  if (!tb[2]) { printf("  no pentOnStation -- nothing to compare\n"); return; }

  fttDb->loadGeometryFromDb(tb[0], tb[1], tb[2]);
  if (!fttDb->geometryFromDb()) { printf("  DB load did not take effect\n"); return; }

  const char *qn[4]; qn[0] = "A"; qn[1] = "B"; qn[2] = "C"; qn[3] = "D";
  double wx = 0, wy = 0, wz = 0;
  printf("\n  quad st    DB x,y (mm)        hardcoded x,y      dx     dy  |"
         "   DB z    hc z     dz (cm)\n");
  for (int iqc = 0; iqc < 4; iqc++) {
    for (int ipc = 0; ipc < 4; ipc++) {
      double dx = fttDb->xShift(iqc, ipc) - hcX[iqc][ipc];
      double dy = fttDb->yShift(iqc, ipc) - hcY[iqc][ipc];
      double dz = fttDb->zLoc (iqc, ipc) - hcZ[iqc][ipc];
      if (fabs(dx) > wx) wx = fabs(dx);
      if (fabs(dy) > wy) wy = fabs(dy);
      if (fabs(dz) > wz) wz = fabs(dz);
      printf("   %s  %d  %9.3f %8.3f   %9.3f %8.3f  %6.3f %6.3f  | %8.3f %8.3f %8.4f\n",
             qn[iqc], ipc + 1,
             fttDb->xShift(iqc, ipc), fttDb->yShift(iqc, ipc),
             hcX[iqc][ipc], hcY[iqc][ipc], dx, dy,
             fttDb->zLoc(iqc, ipc), hcZ[iqc][ipc], dz);
    }
  }
  printf("\n  max |dx| %.4f mm   max |dy| %.4f mm   max |dz| %.5f cm\n", wx, wy, wz);
  if (wx < 1e-4 && wy < 1e-4 && wz < 1e-6)
    printf("  PASS: the DB slot reproduces the hardcoded survey exactly.\n\n");
  else
    printf("  DIFFERENT geometry (expected for the ideal slot, not for the survey slot).\n\n");
}
