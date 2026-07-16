// Classifies StFwdResidualMaker output files by whether they contain a given
// top-level key (default: ProjXYByCharge, added 2026-07-15) -- the only
// reliable way to tell old-library vs new-library output, since a job that
// started before a rebuild but finished writing after it would have a "new"
// mtime with old-library content (already dlopen'd into that process's
// memory at job start -- rebuilding the .so on disk doesn't affect it).
//
// Usage: root4star -b -q 'classifyResidual.C("filelist.txt","newlist.txt","oldlist.txt")'
//    or: root4star -b -q 'classifyResidual.C("filelist.txt","newlist.txt","oldlist.txt","SomeOtherKey")'
void classifyResidual(const char* inList, const char* newList, const char* oldList,
                       const char* checkKey = "ProjXYByCharge") {
    std::ifstream in(inList);
    std::ofstream outNew(newList), outOld(oldList);
    std::string line;
    int nNew = 0, nOld = 0, nBad = 0;
    gErrorIgnoreLevel = kFatal; // silence routine TFile::Open messages
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        TFile* f = TFile::Open(line.c_str());
        if (!f || f->IsZombie()) { nBad++; if (f) delete f; continue; }
        bool hasKey = (f->GetListOfKeys()->Contains(checkKey) != 0);
        f->Close();
        delete f;
        if (hasKey) { outNew << line << "\n"; nNew++; }
        else        { outOld << line << "\n"; nOld++; }
    }
    printf("=== %s: has-%s=%d  missing=%d  unreadable=%d ===\n", inList, checkKey, nNew, nOld, nBad);
}
