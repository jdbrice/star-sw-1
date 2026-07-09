//usr/bin/env root4star -l -b -q  $0; exit $?
// that is a valid shebang to run script as executable


void build_geom( TString geomtag = "y2023", TString output="fGeom.root" ) {

    gSystem->Load( "libStarRoot.so" );
    
    //gROOT->SetMacroPath("/star-sw/StRoot/macros/");
    // gROOT->SetMacroPath(".:./StRoot/macros:./StRoot/macros/graphics:./StRoot/macros/analysis:./StRoot/macros/test:./StRoot/macros/examples:./StRoot/macros/html:./StRoot/macros/qa:./StRoot/macros/calib:./StRoot/macros/mudst:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/graphics:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/analysis:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/test:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/examples:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/html:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/qa:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/calib:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/mudst:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/macros:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/tutorials");
    gROOT->SetMacroPath(".:/star-sw/StRoot/macros/:./StRoot/macros:./StRoot/macros/graphics:./StRoot/macros/analysis:./StRoot/macros/test:./StRoot/macros/examples:./StRoot/macros/html:./StRoot/macros/qa:./StRoot/macros/calib:./StRoot/macros/mudst:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/graphics:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/analysis:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/test:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/examples:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/html:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/qa:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/calib:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/mudst:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/macros:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/tutorials");
    gROOT->LoadMacro("bfc.C");
    bfc(0, "fzin agml sdt20211016", "" );

    gSystem->Load("libStarClassLibrary.so");
    gSystem->Load("libStEvent.so" );

    // Force build of the geometry
    TFile *geom = TFile::Open( output.Data() );

    if ( 0 == geom ) {
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
