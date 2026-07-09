#ifndef ST_FWD_MOCK_MAKER_H
#define ST_FWD_MOCK_MAKER_H

#include "TClonesArray.h"
#ifndef __CINT__
#include "GenFit/Track.h"
#include "StFwdTrackMaker/include/Tracker/FwdHit.h"
#include "StMuDSTMaker/COMMON/StMuFwdTrack.h"
#endif

#include "StChain/StMaker.h"
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "StEvent/StEnumerations.h"
#include "StThreeVectorD.hh"
#include "StThreeVectorF.hh"
#include "StPhysicalHelixD.hh"

#include <map>


class StMuDstMaker;
class StMuDst;
class StMuFwdTrackCollection;
class StMuFcsCollection;
class StFwdTrackMaker;
class StEvent;
class StFwdTrack;

class StFwdMockMaker : public StMaker {

    ClassDef(StFwdMockMaker, 1);

  public:
    StFwdMockMaker() {}
    ~StFwdMockMaker(){/* nada */};

    int Init() {
        LOG_INFO << "StFwdMockMaker::Init()" << endm;
    
        LOG_INFO << "StFwdMockMaker::Init() done" << endm;
        return kStOK;
    }
    int Finish(){
        return kStOk;
    }
    int Make();
    void FillEvent();
    StFwdTrack * makeStFwdTrack( size_t indexTrack );


    void Clear(const Option_t *opts = ""){

    }

    vector<TVector3> mFwdHitsFtt; // FTT hits

};


#endif
