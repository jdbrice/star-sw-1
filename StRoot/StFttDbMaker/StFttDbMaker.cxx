/***************************************************************************
 * StFttDbMaker.cxx
 * jdb & Zhen
 ***************************************************************************
 * Description: This maker is the interface between FTT and the STAR database
 ***************************************************************************/

#include "StFttDbMaker.h"
#include "StFttDb.h"
#include "St_db_Maker/St_db_Maker.h"
#include "StMessMgr.h"
#include "TSystem.h"
#include "TDatime.h"
#include "tables/St_fttDataWindowsB_Table.h"
#include "tables/St_Survey_Table.h"

ClassImp(StFttDbMaker)

StFttDbMaker::StFttDbMaker(const char *name) : StMaker(name){
  LOG_INFO << "******** StFttDbMaker::StFttDbMaker = "<<name<<endm;
  mFttDb = new StFttDb("fttDb");
  AddData(mFttDb,".const");
}; 

StFttDbMaker::~StFttDbMaker() {
  delete mFttDb;
}

int StFttDbMaker::Init(){
  mFttDb->Init();
  return StMaker::Init();
}

void StFttDbMaker::Clear(Option_t *option){
  StMaker::Clear(option);
}

int StFttDbMaker::Make(){
  return StMaker::Make();
}

int StFttDbMaker::InitRun(int runNumber) {
  LOG_INFO << "StFttDbMaker::InitRun - run = " << runNumber << endm;

    // mFttDb->loadHardwareMapFromFile( "StRoot/StFttDbMaker/vmm_map.dat" );
    mFttDb->loadStripCenterFromFile( "StRoot/StFttDbMaker/Row1.txt" );
    mFttDb->loadStripEdgeFromFile(   "StRoot/StFttDbMaker/Row4_edge.txt" );
    mFttDb->loadStripCenterFromFile( "StRoot/StFttDbMaker/Row4.txt" );
    mFttDb->loadStripLengthFromFile( "StRoot/StFttDbMaker/Row1_StripLength.txt" );
    mFttDb->loadStripLengthFromFile( "StRoot/StFttDbMaker/Row2_StripLength.txt" );
    mFttDb->loadStripLengthFromFile( "StRoot/StFttDbMaker/Row3_StripLength.txt" );
    mFttDb->loadStripLengthFromFile( "StRoot/StFttDbMaker/Row4_StripLength.txt" );
    mFttDb->loadStripLengthFromFile( "StRoot/StFttDbMaker/Row5_StripLength.txt" );

    std::ifstream file("vmm_map.dat");
    if(file.is_open()){ // debugging / calibration only
        file.close();
        LOG_INFO << "Loading Hardware Map from FILE!!" << endm;
        LOG_INFO << "Remove / rename file to load from DB" << endm;
        mFttDb->loadHardwareMapFromFile( "vmm_map.dat" );
        // mFttDb->loadHardwareMapFromFile( "/star/u/wangzhen/sTGC/Commissioning/ClusterFinder/PointMaker_building_test_0616/star-sw-1/StRoot/StFwdTrackMaker/macro/vmm_map.dat" );
    } else { // default

        TDataSet *mDbDataSet = GetDataBase("Geometry/ftt/fttHardwareMap");
        if (mDbDataSet){
          St_fttHardwareMap *dataset = (St_fttHardwareMap*) mDbDataSet->Find("fttHardwareMap");
          mFttDb->loadHardwareMapFromDb( dataset );
        } else {
          LOG_WARN << "Cannot access Geometry/ftt/fttHardwareMap and no local map given" << endm;
        }
    }

    loadGeometry();

    loadDataWindows( runNumber );

  return kStOK;
}

// Per-quadrant sTGC offsets from the survey tables, replacing the hardcoded numbers
// in StFttDb unless setUseDbGeometry(false) was called. Missing tables are not an
// error -- StFttDb keeps the hardcoded values and says so.
//
// Namespace is Geometry/stgc for now. Dmitry: stgc is not a DB namespace and this
// moves to Geometry/ftt together with the StgmGeo1.xml <Misalign> paths, agreed for
// October after the GitHub/Gitea migration.
void StFttDbMaker::setUseDbGeometry( bool v ){ if ( mFttDb ) mFttDb->setUseDbGeometry( v ); }

void StFttDbMaker::loadGeometry(){
  const char *nm[3] = { "Geometry/stgc/stgcOnTpc",
                        "Geometry/stgc/stationOnStgc",
                        "Geometry/stgc/pentOnStation" };
  St_Survey *t[3] = { 0, 0, 0 };
  for ( int i = 0; i < 3; i++ ) {
    TDataSet *ds = GetDataBase( nm[i] );
    if ( !ds ) { LOG_WARN << "Cannot access " << nm[i] << endm; continue; }
    const char *leaf = strrchr( nm[i], '/' ) + 1;
    t[i] = (St_Survey *) ds->Find( leaf );
    if ( !t[i] ) LOG_WARN << "No " << leaf << " table in " << nm[i] << endm;
  }
  mFttDb->loadGeometryFromDb( t[0], t[1], t[2] );
}

void StFttDbMaker::loadDataWindows( int runNumber ){
  mFttDb->clearDataWindows();
  // A local per-run file overrides the DB -- for calibration and testing, same
  // convention as vmm_map.dat above. Condor jobs run in a scratch directory without
  // it, so production is unaffected.
  TString localFile = Form( "fttDataWindow/fttDataWindow.%d.txt", runNumber );
  if ( !gSystem->AccessPathName( localFile.Data() ) ) {
    LOG_WARN << "StFttDbMaker::loadDataWindows - LOADING DATA WINDOWS FROM LOCAL FILE " << localFile
             << " (remove/rename it to use the DB)" << endm;
    mFttDb->loadDataWindowsFromFile( localFile.Data() );
  } else {
    TDataSet *mDbDataSetDW = GetDataBase("Calibrations/ftt/fttDataWindowsB");
    if ( mDbDataSetDW ) {
      St_fttDataWindowsB *dataset = (St_fttDataWindowsB*) mDbDataSetDW->Find("fttDataWindowsB");
      mFttDb->loadDataWindowsFromDb( dataset );
      if ( dataset ) {
        TDatime val[2];
        int iv = St_db_Maker::GetValidity( dataset, val );
        // The DB query time is the St_db_Maker's own clock (SetDateTime), NOT this
        // maker's StMaker::GetDateTime(), which reads the event header and is never
        // set in a hand-built MuDst chain.
        StMaker *dbmk = GetMakerInheritsFrom( "St_db_Maker" );
        TString qt = dbmk ? TString( dbmk->GetDateTime().AsSQLString() ) : TString( "(no St_db_Maker)" );
        LOG_INFO << "StFttDbMaker::loadDataWindows - DB entry validity rc=" << iv
                 << " begin " << val[0].AsSQLString() << " end " << val[1].AsSQLString()
                 << " ; DB query time " << qt << endm;
      }
    } else {
      LOG_WARN << "Cannot access Calibrations/ftt/fttDataWindowsB" << endm;
    }
  }

  // Run-stamp guard. VMM anchors are reset at random at every run start, so an entry
  // measured for another run must never supply anchors. St_db_Maker correctly returns
  // the latest entry at or before this run's start -- which for a run without its own
  // entry is the previous run's. The windows (min/max) are still applied either way.
  int stamp = mFttDb->dataWindowStampRun();
  bool valid = ( stamp == runNumber );
  mFttDb->setDataWindowAnchorsValid( valid );
  if ( valid ) {
    LOG_INFO << "StFttDbMaker::loadDataWindows - per-VMM time anchors ENABLED: entry stamped for run "
             << stamp << endm;
  } else if ( stamp < 0 ) {
    LOG_INFO << "StFttDbMaker::loadDataWindows - unstamped (legacy) data-window entry: windows used, "
             << "anchors ignored; StFttHitCalibMaker derives anchors on the fly" << endm;
  } else {
    LOG_WARN << "StFttDbMaker::loadDataWindows - data-window entry is stamped for run " << stamp
             << " but this is run " << runNumber << ": its anchors are IGNORED (they are not valid "
             << "across runs); StFttHitCalibMaker derives anchors on the fly" << endm;
  }
}