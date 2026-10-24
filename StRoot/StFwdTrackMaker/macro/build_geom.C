//usr/bin/env root4star -l -b -q  $0; exit $?
// that is a valid shebang to run script as executable


// misalign = true adds the BFC "misalign" option, which is the ONLY thing that
// makes StBFChain call AgPosition::SetReal() (StBFChain.cxx:1011) and so the only
// way the <Misalign> tags in the XML reach the built geometry. Without it
// AgMLDbFunctor's default operator() returns the unit matrix and the geometry is
// nominal however many Misalign tags the module carries.
// sdt picks WHICH tables: the misalignment is read at that timestamp.
//
// sdt is a date AND a time, "sdtYYYYMMDD[.HHMMSS]", and the trailing field of a
// StarDb file name is that time, not a version counter: pentOnStation.20211112
// .000004.C begins at 2021-11-12 00:00:04. So plain "sdt20211112" is midnight and
// sorts BEFORE it -- it silently selects the newest table at or before midnight
// (pentOnStation.20211110.000001.C, i.e. the nominal sTGC) and you get a cache
// with no alignment in it and no warning. To pick up the .000004 lane you must
// pass sdt20211112.000004. Check the log for the "LoadTable: .L .../StarDb/..."
// lines and confirm the file names are the ones you meant.
void build_geom( TString geomtag = "y2023", TString output="fGeom.root",
                 bool misalign = false, TString sdt = "sdt20211016" ) {

    gSystem->Load( "libStarRoot.so" );
    
    //gROOT->SetMacroPath("/star-sw/StRoot/macros/");
    // gROOT->SetMacroPath(".:./StRoot/macros:./StRoot/macros/graphics:./StRoot/macros/analysis:./StRoot/macros/test:./StRoot/macros/examples:./StRoot/macros/html:./StRoot/macros/qa:./StRoot/macros/calib:./StRoot/macros/mudst:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/graphics:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/analysis:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/test:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/examples:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/html:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/qa:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/calib:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/mudst:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/macros:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/tutorials");
    // The hardcoded /afs/rhic.bnl.gov path below does not exist on every node and
    // overrode the one the STAR environment already sets, so bfc.C was not found.
    // Upstream removed it the same way in PR #33; keeping the string commented out
    // for anyone running outside the standard environment.
    // gROOT->SetMacroPath(".:/star-sw/StRoot/macros/:./StRoot/macros:./StRoot/macros/graphics:./StRoot/macros/analysis:./StRoot/macros/test:./StRoot/macros/examples:./StRoot/macros/html:./StRoot/macros/qa:./StRoot/macros/calib:./StRoot/macros/mudst:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/graphics:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/analysis:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/test:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/examples:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/html:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/qa:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/calib:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/mudst:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/macros:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/tutorials");
    gROOT->LoadMacro("bfc.C");
    TString bfcOpt = "fzin agml " + sdt;
    if ( misalign ) bfcOpt += " misalign";
    cout << "build_geom: bfc options [" << bfcOpt.Data() << "]" << endl;
    bfc(0, bfcOpt.Data(), "" );

    gSystem->Load("libStarClassLibrary.so");
    gSystem->Load("libStEvent.so" );

    // Force build of the geometry
    TFile *geom = TFile::Open( output.Data() );

    if ( 0 == geom ) {
        // Pre-warm CINT with the Survey table wrapper. The <Misalign> tags make
        // AgPosition ask St_db_Maker for Geometry/{fst,stgc}/* during Construct()
        // below, and St_db_Maker reads the StarDb/*.C table macros through CINT.
        // By that point CINT can fail to parse tables/St_Survey_Table.h
        // ("Error: unrecognized language construct ... St_Survey_Table.h:20",
        // the ClassDefTable line), and the failed load then trips
        // St_db_Maker::LoadTable's assert(!ee) and aborts the whole build. The
        // same header parses fine right after bfc(), so parsing it once here --
        // while the interpreter is still clean -- avoids it. Verified harmless:
        // with and without this line, dev2022sm/sdt20211112.000004 produces
        // caches whose 16 STFM and 108 FTUS node positions agree exactly.
        gROOT->ProcessLine("#include \"tables/St_Survey_Table.h\"");

        AgModule::SetStacker( new StarTGeoStacker() );
        AgPosition::SetDebug(2);
        cout << "Building geometry for tag [" << geomtag.Data() << "]" << endl;
        StarGeometry::Construct( geomtag.Data() );

        // Replace FCS ECAL (WcalGeo1_*) and FCS HCAL (HcalGeo1_*) materials with air
        // so that Genfit can extrapolate to z=725 cm without hitting calorimeter material.
        // Field is already zero at |z|>450 cm (STARField.h), so zero-material + zero-field
        // means straight-line propagation, which Genfit handles trivially.
        TGeoMaterial *airMat = new TGeoMaterial("AIR_FCS", 14.61, 7.3, 1.205e-3);
        TGeoMedium   *airMed = new TGeoMedium("AIR_FCS", gGeoManager->GetListOfMedia()->GetEntries()+1, airMat);
        TObjArray *vols = gGeoManager->GetListOfVolumes();
        int nReplaced = 0;
        for (int iv = 0; iv < vols->GetEntries(); iv++) {
            TGeoVolume *v = (TGeoVolume*)vols->At(iv);
            TString mname = v->GetMaterial()->GetName();
            if (mname.BeginsWith("WcalGeo") || mname.BeginsWith("HcalGeo")) {
                v->SetMedium(airMed);
                nReplaced++;
            }
        }
        cout << "Replaced " << nReplaced << " FCS calorimeter volumes with air for tracking geometry." << endl;

        // Genfit requires the geometry is cached in a ROOT file
        gGeoManager->Export( output.Data() );
        cout << "Writing output to geometry file [" << output.Data() << "]" << endl;
    }
    else {
        cout << "WARNING:  Geometry file [" << output.Data() << "] already exists." << endl;
        cout << "Existting without doing anything!" << endl;
        delete geom;
    }

}
