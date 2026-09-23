//usr/bin/env root4star -l -b -q  $0; exit $?
// that is a valid shebang to run script as executable

void build_geom_misalign( TString geomtag = "dev2022sm", TString output="fGeom.root", int misalign=1, const char* sdt="20211110.000001") {

    gSystem->Load( "libStarRoot.so" );

    // The "misalign" option makes AgPosition ask St_db_Maker for the survey
    // tables, which live as CINT macros under $STAR/StarDb/Geometry/{fst,stgc}/.
    // Each one does #include "tables/St_Survey_Table.h", and CINT cannot parse
    // that header's ClassDefTable(St_Survey,Survey_st) unless the dictionary is
    // already loaded -- without it the include dies with "unrecognized language
    // construct", LoadTable returns an error and St_db_Maker.cxx:932 asserts.
    // build_geom.C never hits this because with no misalign the tables are
    // never read. Must be loaded BEFORE bfc/StarGeometry::Construct.
    gSystem->Load( "libSt_base.so" );
    gSystem->Load( "libStDb_Tables.so" );
    if ( !TClass::GetClass("St_Survey") ) {
        cout << "FATAL: St_Survey dictionary not loaded -- survey tables would fail to parse" << endl;
        return;
    }
    

    // No SetMacroPath needed: root4star's default macro path already contains
    // $STAR/StRoot/macros, which is where bfc.C lives (verified: LoadMacro
    // resolves it with no path set). Setting it REPLACES the default, and the
    // path that used to be here also named /star-sw/... and /afs/rhic.bnl.gov/...
    // directories that do not exist on this machine.
    gROOT->LoadMacro("bfc.C");

    if(misalign==0){
      bfc(0, Form("fzin agml sdt%s",sdt), "" );
    }else{
      bfc(0, Form("fzin agml sdt%s misalign",sdt), "" );
    }
    output.ReplaceAll(".root",Form(".%s.misalign%1d.sdt%s.root",geomtag.Data(),misalign,sdt));
    
    gSystem->Load("libStarClassLibrary.so");
    gSystem->Load("libStEvent.so" );

    // Parse the survey table header NOW, while the interpreter is idle.
    // St_db_Maker reads StarDb/Geometry/{fst,stgc}/*.C with a nested ".L" from
    // inside AgML module execution, and in that state CINT cannot parse
    // St_Survey_Table.h:20 (ClassDefTable) -- it dies with "unrecognized
    // language construct", LoadTable returns an error and St_db_Maker.cxx:932
    // asserts. Parsing it here sets the include guard so the nested .L never
    // has to parse it again. Verified: without this, construction aborts on the
    // first table; with it, all 7 survey tables load.
    TInterpreter::EErrorCode incErr;
    gInterpreter->ProcessLine("#include \"tables/St_Survey_Table.h\"", &incErr);
    if ( incErr != TInterpreter::kNoError ) {
        cout << "FATAL: could not pre-parse tables/St_Survey_Table.h" << endl;
        return;
    }

    // Force build of the geometry
    TFile *geom = TFile::Open( output.Data() );

    if ( 0 == geom ) {
        AgModule::SetStacker( new StarTGeoStacker() );
        AgPosition::SetDebug(2);
        cout << "Building geometry for tag [" << geomtag.Data() << "]" << endl;
        StarGeometry::Construct( geomtag.Data() );

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
