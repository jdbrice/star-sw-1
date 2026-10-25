// fstSensorPhiZ.C -- for every FST silicon node, the azimuth of the CENTRE of the
// segment and its z. The node translation is useless for azimuth: FTUS is a TUBS
// placed at x=y=0 with the azimuth carried by the rotation, so atan2(ty,tx) is 0.
// Take the shape's mid-phi / mid-r in the LOCAL frame and transform that instead.
void fstSensorPhiZ(const char* fin, const char* fout)
{
  TGeoManager *gm = TGeoManager::Import(fin);
  if (!gm) { printf("fstSensorPhiZ: cannot import %s\n", fin); return; }
  FILE *fp = fopen(fout, "w");
  if (!fp) { printf("fstSensorPhiZ: cannot write %s\n", fout); return; }

  TGeoIterator it(gm->GetTopVolume());
  TGeoNode *nd = 0;
  int nFound = 0;
  while ((nd = it.Next())) {
    TString vn = nd->GetVolume()->GetName();
    if (!vn.BeginsWith("FTUS")) continue;
    TGeoTubeSeg *ts = (TGeoTubeSeg *) nd->GetVolume()->GetShape();
    if (!ts) continue;
    double rmid  = 0.5 * (ts->GetRmin() + ts->GetRmax());
    double p1    = ts->GetPhi1();
    double p2    = ts->GetPhi2();
    double pmid  = 0.5 * (p1 + p2);
    double lp[3], gp[3];
    lp[0] = rmid * TMath::Cos(pmid * TMath::DegToRad());
    lp[1] = rmid * TMath::Sin(pmid * TMath::DegToRad());
    lp[2] = 0.0;
    const TGeoHMatrix *mx = it.GetCurrentMatrix();
    mx->LocalToMaster(lp, gp);
    double phig = TMath::ATan2(gp[1], gp[0]) * TMath::RadToDeg();
    if (phig < 0) phig += 360.0;
    fprintf(fp, "%s %.4f %.4f %.4f %.4f\n", nd->GetName(), phig, gp[2], rmid, ts->GetDz());
    nFound++;
  }
  fclose(fp);
  printf("fstSensorPhiZ: %s -> %s : %d FTUS nodes\n", fin, fout, nFound);
}
