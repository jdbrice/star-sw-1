//usr/bin/env root4star -l -b -q  $0; exit $?
// that is a valid shebang to run script as executable

// build_geom_misalign.C -- THE builder for the StFwdTrackMaker geometry cache.
//
// Not to be confused with StRoot/StFwdTrackMaker/macro/build_geom.C, which
// reproduces none of our caches: it hardcodes sdt20211016, has no misalign, and
// takes no arguments.
//
// This supersedes the former build_geom_run24.C, which was a near-copy of this
// file differing in TWO ways neither filename recorded:
//
//   this file (before)   misalign configurable, FCS calorimeter material LEFT
//                        IN PLACE, St_Survey_Table.h pre-parsed
//   build_geom_run24.C   misalign and sdt20240622 hardcoded, FCS calorimeter
//                        material REPLACED WITH AIR, no pre-parse
//
// So every Run22 cache was built with FCS material present and every Run24
// cache without it -- measured with TGeoManager::Import: the Run22 cache keeps
// 7 WcalGeo + 8 HcalGeo volumes, the Run24 cache has all 15 as AIR_FCS. That
// difference sat inside the Run22-vs-Run24 comparison with nothing naming it.
// Here it is the fcsAir argument, and it goes into the output filename.
//
// ARGUMENTS
//   geomtag   AGML tag passed to StarGeometry::Construct.
//               dev2022sm  Run22        y2024  Run24
//   sdt       DB timestamp for the bfc "sdt" option, YYYYMMDD[.HHMMSS]. This
//             selects the survey/misalign tables: the latest entry in
//             StarDb/Geometry/{fst,stgc} with beginTime <= sdt wins.
//   misalign  1 = add the bfc "misalign" option so AgPosition reads the survey
//             tables. 0 = ideal geometry, tables never read.
//   fcsAir    1 = replace the FCS ECAL/HCAL materials with air, so Genfit can
//             extrapolate to the shower max at z=725 cm without crossing
//             calorimeter material. The field is already zero at |z|>450 cm
//             (STARField.h), so zero-material + zero-field is straight-line
//             propagation, which Genfit handles trivially.
//             0 = leave the calorimeter material in the geometry.
//   label     Optional extra token for the filename ONLY. Use it for anything
//             the other arguments cannot express -- above all WHICH ITERATION
//             of the DB tables you got, because sdt does not say. sdt20240622
//             picks pentOnStation.20231210.000001 today and picked .000000
//             before that row was added, so two caches built from the same
//             arguments on different days are different geometries. Pass
//             "iter1", "seeded", etc.
//
// OUTPUT NAME is always constructed from the arguments, never passed in:
//
//     fGeom.<geomtag>.misalign<N>.fcs<F>.sdt<sdt>[.<label>].root
//
// e.g. Run22 as every Run22 production used it (NO air):
//        build_geom_misalign.C("dev2022sm","20211112.000003",1,0)
//        -> fGeom.dev2022sm.misalign1.fcs0.sdt20211112.000003.root
//
//      Run24 iteration 1, as the run24align/run24pickup campaigns used it:
//        build_geom_misalign.C("y2024","20240622",1,1,"iter1")
//        -> fGeom.y2024.misalign1.fcs1.sdt20240622.iter1.root
//
// The caches already on disk predate this naming and keep their old names
// (fGeom.dev2022sm.misalign1.sdt20211112.000003.root,
//  fGeom.y2024.iter1.sdt20240622.root). Content is the same; only the name
// differs, so there is no need to rebuild them.
//
// An existing output file is NEVER overwritten -- the macro reports it and
// exits. Delete it first, or change an argument.

void build_geom_misalign( TString     geomtag  = "dev2022sm",
                          const char* sdt      = "20211112.000003",
                          int         misalign = 1,
                          int         fcsAir   = 1,
                          TString     label    = "" )
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

    TString output = Form( "fGeom.%s.misalign%1d.fcs%1d.sdt%s",
                           geomtag.Data(), misalign, fcsAir, sdt );
    if ( label != "" ) { output += "."; output += label; }
    output += ".root";

    // Check this BEFORE constructing anything -- the build takes minutes.
    TFile *exists = TFile::Open( output.Data() );
    if ( 0 != exists ) {
        cout << "=== [" << output.Data() << "] ALREADY EXISTS ===" << endl;
        cout << "=== nothing was built. delete it, or change an argument. ===" << endl;
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
    cout << "=== build_geom_misalign done ===" << endl;
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
