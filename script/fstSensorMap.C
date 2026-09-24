// fstSensorMap.C -- does global sensor index -> FTUS path (as getFstSensorOrigin builds it)
// land on the sensor the HIT thinks it is?  Hit sensor 0 = inner (r 5-16.5), 1,2 = outer.
void fstSensorMap(){
  TGeoManager::Import("/star/u/akio/fcstrk11/star-sw-fwd/fGeom.root");
  const int kElecToGeantWedge[12]      = {2, 7, 1, 12, 6, 11, 5, 10, 4, 9, 3, 8};
  const int kElecToGeantWedgeDisk2[12] = {7, 1, 12, 6, 11, 5, 10, 4, 9, 3, 8, 2};
  printf(">>> globalIdx = disk*36 + wedge*3 + sensor  (sensor 0=inner per the HIT)\n");
  printf(">>> idx  disk wedge sens   FTUS path                              rmin   shapePhi        z\n");
  for (int d=0; d<1; d++) for (int w=0; w<2; w++) for (int sn=0; sn<3; sn++){
    int index = d*36 + w*3 + sn;
    int sensorIndex = (index % 3) + 1;
    int electronicWedge = (index / 3) % 12;
    int planeIndex = (index / 36) + 4;
    const int* wmap = (planeIndex==5) ? kElecToGeantWedgeDisk2 : kElecToGeantWedge;
    int wedgeIndex = wmap[electronicWedge];
    TString path = Form("/HALL_1/CAVE_1/FSTM_1/FSTD_%d/FSTW_%d/FTUS_%d", planeIndex, wedgeIndex, sensorIndex);
    if (!gGeoManager->cd(path)) { printf(">>> %3d  cd FAILED %s\n", index, path.Data()); continue; }
    TGeoNode* nd = gGeoManager->GetCurrentNode();
    TGeoShape* sh = nd->GetVolume()->GetShape();
    double rmin=-1, p1=-999, p2=-999;
    if (sh->InheritsFrom("TGeoTubeSeg")){
      rmin=((TGeoTubeSeg*)sh)->GetRmin(); p1=((TGeoTubeSeg*)sh)->GetPhi1(); p2=((TGeoTubeSeg*)sh)->GetPhi2();
    }
    double z = gGeoManager->GetCurrentMatrix()->GetTranslation()[2];
    printf(">>> %3d   %d    %2d    %d   FTUS_%d  rmin=%5.2f  [%7.2f,%7.2f]  z=%9.4f  %s\n",
           index, d, w, sn, sensorIndex, rmin, p1, p2, z,
           (sn==0) ? ((rmin<10)?"inner ok":"<-- HIT says inner, GEOM is OUTER")
                   : ((rmin>10)?"outer ok":"<-- HIT says outer, GEOM is INNER"));
  }
}
