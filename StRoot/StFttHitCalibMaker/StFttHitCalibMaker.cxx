#include "StFttHitCalibMaker.h"

// #include "StFttRawHitMaker/StFttRawHitMaker.h"

#include "StEvent/StFttRawHit.h"
#include "StEvent/StEvent.h"
#include "StEvent/StFttCollection.h"

#include "StEvent/StFttCluster.h"

#include "StFttDbMaker/StFttDb.h"

#include "TFile.h"
#include "TCanvas.h"

#include <set>
//_____________________________________________________________                                                       
StFttHitCalibMaker::StFttHitCalibMaker(const char *name):
StMaker("fttHitCalib",name), mHelper(nullptr)
{                                            
    LOG_DEBUG << "StFttHitCalibMaker::ctor"  << endm;
    this->mCalibMode = StFttHitCalibMaker::CalibMode::Production;
    mHelper = new HitCalibHelper();
}
//_____________________________________________________________                                                       
StFttHitCalibMaker::~StFttHitCalibMaker()
{
  if(mHelper) delete mHelper;
  mHelper = nullptr;
}

//_____________________________________________________________                                                       
Int_t StFttHitCalibMaker::Init()
{
    LOG_INFO << "StFttHitCalibMaker::Init" << endm;
    return kStOk;
}
//_____________________________________________________________                                                       
Int_t StFttHitCalibMaker::InitRun(Int_t runnumber)
{ 
    mHelper->clear();
    mNHitsDbAnchor = 0; mNHitsOnTheFly = 0; mNHitsNotReady = 0; mBookkeepingPrinted = false;
    return kStOk;
}

//_____________________________________________________________                                                       
// How each hit's time was obtained. Printed from FinishRun AND Finish, once: in a
// hand-built MuDst chain the chain run number is never set, so StMaker::Finish()
// does not call FinishRun, and fwd_afterburner_db.C does not call chain->Finish()
// at all. There the reliable indicator is StFttDbMaker's InitRun line
// ("per-VMM time anchors ENABLED" / "IGNORED").
void StFttHitCalibMaker::printTimeBookkeeping( const char *where )
{
    if ( mBookkeepingPrinted ) return;
    if ( mNHitsDbAnchor + mNHitsOnTheFly + mNHitsNotReady == 0 ) return;
    LOG_INFO << "StFttHitCalibMaker::" << where << " hit times: from DB anchor " << mNHitsDbAnchor
             << ", on-the-fly anchor " << mNHitsOnTheFly
             << ", not calibrated (time=-4097) " << mNHitsNotReady << endm;
    mBookkeepingPrinted = true;
}

Int_t StFttHitCalibMaker::FinishRun(Int_t runnumber)
{ 
    printTimeBookkeeping( Form( "FinishRun(%d)", runnumber ) );
    mHelper->clear();
    return kStOk;
}

//-------------------------------------------------------------                                                       
Int_t StFttHitCalibMaker::Finish()
{ 
    LOG_INFO << "StFttHitCalibMaker::Finish()" << endm;
    printTimeBookkeeping( "Finish" );

    if (this->mCalibMode == StFttHitCalibMaker::CalibMode::Calibration) {
        LOG_INFO << "Writing StFttHitCalib parameters to plaintext: " << endm;
        WriteCalibrationToPlainText();
    }

    return kStOk;
}


void StFttHitCalibMaker::WriteCalibrationToPlainText() {

    ofstream outf( "fttRawHitTime.dat" );
    for ( int uuid = 0; uuid <= 400; uuid++ ){
        Short_t anchor = mHelper->anchor( uuid );
        auto hist = mHelper->histFor( uuid );
        size_t counts = hist.size(); 
        size_t samples = mHelper->samples( uuid );
        outf << TString::Format( "%d\t%d\t%lu\t%lu", uuid, (int) anchor, counts, samples ) << endl;
    }
    outf.close();

}

//_____________________________________________________________                                                       
Int_t StFttHitCalibMaker::Make()
{ 
    LOG_INFO << "StFttHitCalibMaker::Make()" << endm;

    mEvent = (StEvent*)GetInputDS("StEvent");
    if(mEvent) {
        LOG_DEBUG<<"Found StEvent"<<endm;
    } else {
        return kStOk;
    }
    mFttCollection=mEvent->fttCollection();
    if(!mFttCollection) {
        return kStOk;
    } else {
        LOG_DEBUG <<"Found StFttCollection"<<endm;
    }

    mFttDb = static_cast<StFttDb*>(GetDataSet("fttDb"));


    for ( auto rawHit : mFttCollection->rawHits() ) {

        // Per-run anchor from Calibrations/ftt/fttDataWindowsB, when StFttDbMaker has
        // confirmed the entry is stamped for THIS run. Removes the per-file warm-up,
        // which on sparse data (zero-field physics stream) never converges.
        Short_t dbAnchor = 0;
        if ( mFttDb->getAnchor( rawHit, dbAnchor ) ) {
            // same definition as HitCalibHelper::time(): shortest signed distance on
            // the 12-bit circular dbcid counter
            Short_t diff = rawHit->dbcid() - dbAnchor;
            if ( diff > 2048 ) diff -= 4096;
            if ( diff < -2048 ) diff += 4096;
            rawHit->setTime( diff );
            mNHitsDbAnchor++;
            continue;
        }

        // Otherwise derive the anchor on the fly (legacy behaviour).
        UShort_t fob = (UShort_t)mFttDb->fob( rawHit );
        UShort_t uuid = rawHit->vmm() + ( StFttDb::nVMMPerFob * fob );

        mHelper->fill( uuid, rawHit->dbcid() );

        if ( mHelper->ready( uuid ) ){
            rawHit->setTime( mHelper->time( uuid, rawHit->dbcid()) );
            mNHitsOnTheFly++;
        } else {
            rawHit->setTime( -4097 );
            mNHitsNotReady++;
        }
    }

    return kStOk;
}