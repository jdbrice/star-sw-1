// Are the two OUTER sensors distinguishable by their placement matrix alone,
// or only by their shape?
void fstTwoOuter(){
  TGeoManager::Import("/star/u/akio/fcstrk11/star-sw-fwd/fGeom.root");
  for (int s=1; s<=3; s++){
    TString path = Form("/HALL_1/CAVE_1/FSTM_1/FSTD_4/FSTW_2/FTUS_%d", s);
    if(!gGeoManager->cd(path)){ printf(">>> cd failed %s\n",path.Data()); continue; }
    TGeoHMatrix* m = gGeoManager->GetCurrentMatrix();
    const double* t = m->GetTranslation();
    const double* r = m->GetRotationMatrix();
    TGeoTubeSeg* ts = (TGeoTubeSeg*)gGeoManager->GetCurrentNode()->GetVolume()->GetShape();
    printf(">>> FTUS_%d  T=(%.4f,%.4f,%9.4f)  R=[%+.3f %+.3f %+.3f | %+.3f %+.3f %+.3f]  shapePhi=[%7.2f,%7.2f] rmin=%5.2f\n",
           s, t[0],t[1],t[2], r[0],r[1],r[2], r[3],r[4],r[5], ts->GetPhi1(), ts->GetPhi2(), ts->GetRmin());
  }
}
