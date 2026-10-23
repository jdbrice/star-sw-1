// dumpGeomNodes.C -- write the global position of every sTGC pentagon and FST
// sensor node in a geometry cache, one per line, for differencing two caches.
//
//   root4star -l -b -q 'dumpGeomNodes.C("fGeom_dev2026_ideal.root","ideal.txt")'
//
// Run once per file in its OWN process: gGeoManager is a singleton and importing
// a second geometry into the same session is not reliable.
// CINT: fixed arrays, unique loop variable names (see CLAUDE.md).
void dumpGeomNodes(const char* fin = "fGeom.root", const char* fout = "nodes.txt")
{
  TGeoManager *gm = TGeoManager::Import(fin);
  if (!gm) { printf("dumpGeomNodes: cannot import %s\n", fin); return; }

  FILE *fp = fopen(fout, "w");
  if (!fp) { printf("dumpGeomNodes: cannot write %s\n", fout); return; }

  TGeoIterator it(gm->GetTopVolume());
  TGeoNode *nd = 0;
  int nStfm = 0, nFtus = 0;
  while ((nd = it.Next())) {
    TString vn = nd->GetVolume()->GetName();
    bool isPent = vn.BeginsWith("STFM");
    bool isSens = vn.BeginsWith("FTUS");
    if (!isPent && !isSens) continue;
    const TGeoHMatrix *mtx = it.GetCurrentMatrix();
    const Double_t *tr = mtx->GetTranslation();
    // Node names (FTUS_1..FTUS_108, STFM_1..STFM_16) are the key. NOT
    // TGeoIterator::GetPath(TString&) -- in CINT that call kills the loop
    // silently and the file comes out empty.
    fprintf(fp, "%s %.6f %.6f %.6f\n", nd->GetName(), tr[0], tr[1], tr[2]);
    if (isPent) nStfm++;
    if (isSens) nFtus++;
  }
  fclose(fp);
  printf("dumpGeomNodes: %s -> %s : %d STFM nodes, %d FTUS nodes\n",
         fin, fout, nStfm, nFtus);
}
