//usr/bin/env root4star -l -b -q  $0; exit $?
// that is a valid shebang to run script as executable

// build_fgeom.C -- THE builder for the StFwdTrackMaker geometry cache.
//
// Supersedes build_geom_misalign.C and build_geom_run24.C, which were two
// near-copies that differed in TWO ways neither filename recorded:
//
//   build_geom_misalign.C   misalign configurable, FCS calorimeter material
//                           LEFT IN PLACE, St_Survey_Table.h pre-parsed
//   build_geom_run24.C      misalign + sdt20240622 hardcoded, FCS calorimeter
//                           material REPLACED WITH AIR, no pre-parse
//
// So every Run22 cache was built with FCS material present and every Run24
// cache without it -- measured with TGeoManager::Import: the Run22 cache keeps
// 7 WcalGeo + 8 HcalGeo volumes, the Run24 cache has all 15 as AIR_FCS. That
// difference sat inside the Run22-vs-Run24 comparison with nothing naming it.
// Here it is the fcsAir argument, and it goes into the output filename.
//
// Arguments
//   geomtag   AGML tag passed to StarGeometry::Construct.
//               dev2022sm  Run22        y2024  Run24
//   sdt       DB timestamp for the bfc "sdt" option, YYYYMMDD[.HHMMSS]. This
//             selects the survey/misalign tables: the latest entry in
//             StarDb/Geometry/{fst,stgc} with beginTime <= sdt wins. It does
//             NOT identify which iteration you get -- sdt20240622 picks
//             pentOnStation.20231210.000001 today and picked .000000 before
//             that table was added. Use the output name to record the
//             iteration, or you cannot tell two caches apart.
//   misalign  1 = add the bfc "misalign" option so AgPosition reads the survey
//             tables. 0 = ideal geometry, tables never read.
//   fcsAir    1 = replace the FCS ECAL/HCAL materials with air so Genfit can
//             extrapolate to the shower max at z=725 cm without crossing
//             calorimeter material. The field is already zero at |z|>450 cm
//             (STARField.h), so zero-material + zero-field is straight-line
//             propagation, which Genfit handles trivially.
//             0 = leave the calorimeter material in the geometry.
//   output    "" (default) auto-names the file from everything above, so two
//             caches that differ in content cannot share a name:
//                 fGeom.<geomtag>.misalign<N>.<air|mat>.sdt<sdt>.root
//             Pass it explicitly to reproduce a legacy name exactly.
//
// To reproduce the caches already in use, pass the legacy name:
//
//   Run22, as every Run22 production used it (NO air):
//     build_fgeom.C("dev2022sm","20211112.000003",1,0,
//                   "fGeom.dev2022sm.misalign1.sdt20211112.000003.root")
//
//   Run24 iteration 1, as the run24align/run24pickup campaigns used it (air):
//     build_fgeom.C("y2024","20240622",1,1,
//                   "fGeom.y2024.iter1.sdt20240622.root")
//
// An existing output file is NEVER overwritten -- the macro reports it and
// exits. Delete it first, or pass a new name.

void build_fgeom( TString     geomtag  = "dev2022sm",
                  const char* sdt      = "20211112.000003",
                  int         misalign = 1,
                  int         fcsAir   = 1,
                  TString     output   = "" )
{
    gSystem->Load( "libStarRoot.so" );

    // The "misalign" option makes AgPosition ask St_db_Maker for the survey
    // tables, which live as CINT macros under StarDb/Geometry/{fst,stgc}/.
    // Each one does #include "tables/St_Survey_Table.h", and CINT cannot parse
    // that header's ClassDefTable(St_Survey,Survey_st) unless the dictionary is
    // already loaded -- without it the include dies with "unrecognized language
    // construct", LoadTable returns an error and St_db_Maker.cxx:932 asserts.
    // Must be loaded BEFORE bfc/StarGeometry::Construct.
    gSystem->Load( "libSt_base.so" );
    gSystem->Load( "libStDb_Tables.so" );
    if ( misalign != 0 && !TClass::GetClass("St_Survey") ) {
        cout << "FATAL: St_Survey dictionary not loaded -- survey tables would fail to parse" << endl;
        return;
    }

    if ( output == "" ) {
        output = Form( "fGeom.%s.misalign%1d.%s.sdt%s.root",
                       geomtag.Data(), misalign, (fcsAir ? "air" : "mat"), sdt );
    }

    // Check this BEFORE constructing anything -- the build takes minutes.
    TFile *exists = TFile::Open( output.Data() );
    if ( 0 != exists ) {
        cout << "=== build_fgeom: [" << output.Data() << "] ALREADY EXISTS ===" << endl;
        cout << "=== nothing was built. delete it or pass a different output name. ===" << endl;
        delete exists;
        return;
    }

    // No hardcoded macro path: the old copies named /star-sw/... and
    // /afs/rhic.bnl.gov/... directories that do not exist on this machine, and
    // SetMacroPath REPLACES the default rather than adding to it. Load bfc.C by
    // absolute path so it cannot resolve to something else.
    gROOT->SetMacroPath( Form(".:%s/StRoot/macros:%s",
                              gSystem->Getenv("STAR"), gROOT->GetMacroPath()) );
    gROOT->LoadMacro( Form("%s/StRoot/macros/bfc.C", gSystem->Getenv("STAR")) );

    if ( misalign == 0 ) bfc( 0, Form("fzin agml sdt%s",          sdt), "" );
    else                 bfc( 0, Form("fzin agml sdt%s misalign", sdt), "" );

    gSystem->Load( "libStarClassLibrary.so" );
    gSystem->Load( "libStEvent.so" );

    // Parse the survey table header NOW, while the interpreter is idle.
    // St_db_Maker reads StarDb/Geometry/{fst,stgc}/*.C with a nested ".L" from
    // inside AgML module execution, and in that state CINT cannot parse
    // St_Survey_Table.h:20 (ClassDefTable). Parsing it here sets the include
    // guard so the nested .L never has to parse it again. Verified: without
    // this, construction aborts on the first table; with it, all 7 load.
    // build_geom_run24.C omitted this and happened to get away with it.
    if ( misalign != 0 ) {
        TInterpreter::EErrorCode incErr;
        gInterpreter->ProcessLine( "#include \"tables/St_Survey_Table.h\"", &incErr );
        if ( incErr != TInterpreter::kNoError ) {
            cout << "FATAL: could not pre-parse tables/St_Survey_Table.h" << endl;
            return;
        }
    }

    AgModule::SetStacker( new StarTGeoStacker() );
    AgPosition::SetDebug(2);
    cout << "Building geometry for tag [" << geomtag.Data() << "] sdt" << sdt
         << " misalign=" << misalign << " fcsAir=" << fcsAir << endl;
    StarGeometry::Construct( geomtag.Data() );

    int nReplaced = 0;
    if ( fcsAir != 0 ) {
        TGeoMaterial *airMat = new TGeoMaterial("AIR_FCS", 14.61, 7.3, 1.205e-3);
        TGeoMedium   *airMed = new TGeoMedium("AIR_FCS",
                                   gGeoManager->GetListOfMedia()->GetEntries()+1, airMat);
        TObjArray *vols = gGeoManager->GetListOfVolumes();
        for ( int iv = 0; iv < vols->GetEntries(); iv++ ) {
            TGeoVolume *v = (TGeoVolume*)vols->At(iv);
            TString mname = v->GetMaterial()->GetName();
            if ( mname.BeginsWith("WcalGeo") || mname.BeginsWith("HcalGeo") ) {
                v->SetMedium( airMed );
                nReplaced++;
            }
        }
    }

    // Genfit requires the geometry is cached in a ROOT file
    gGeoManager->Export( output.Data() );

    cout << endl;
    cout << "=== build_fgeom done ===" << endl;
    cout << "  output     " << output.Data() << endl;
    cout << "  geomtag    " << geomtag.Data() << endl;
    cout << "  sdt        " << sdt << endl;
    cout << "  misalign   " << misalign << (misalign ? "  (survey tables read)" : "  (ideal geometry)") << endl;
    cout << "  fcsAir     " << fcsAir   << "  (" << nReplaced
         << " FCS calorimeter volumes replaced with air)" << endl;
    if ( fcsAir != 0 && nReplaced == 0 ) {
        cout << "  WARNING: fcsAir was requested but NOTHING was replaced --" << endl;
        cout << "  the geomtag probably has no FCS. Check before using this cache." << endl;
    }
}
