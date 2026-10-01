#ifndef STFTTHITCALIBMAKER_HELPER_H
#define STFTTHITCALIBMAKER_HELPER_H

#include <map>
#include <vector>

// Minimum number of HITS a VMM must have seen before its on-the-fly anchor is
// trusted.
#define MIN_BCID_SAMPLES 200

class HitCalibHelper {
public:
    HitCalibHelper(){

    }

    // FIX (2026-09-14): this used to test dbcidHist[uuid].size() > 200, i.e. the
    // number of DISTINCT dbcid values seen, not the number of hits. A well-timed
    // VMM concentrates its hits on a handful of dbcid values and so could stay
    // "not ready" indefinitely, while a noisy VMM spreading over many values
    // became ready fast -- the criterion rewarded noise. Harmless while production
    // ran with an AcceptAll time cut; decisive once the calibrated cut was on.
    // Measured on zero-field run 23063028 after 300 events: 4 of 373 VMMs ready
    // (horizontal strips 0/128), ~21 distinct values per VMM against ~212 hits.
    // Field-on run 23081049 reached the old threshold in <100 events only because
    // it has 31x more FTT hits per event.
    bool ready( UShort_t uuid ){
        return ( nSamples.count( uuid ) > 0 && nSamples[uuid] > MIN_BCID_SAMPLES );
    }

    void fill( UShort_t uuid, Short_t dbcid ){
        // dbcid is a 12-bit counter, so a VMM's histogram has at most 4096 keys
        // and needs no cap. (The previous cap -- stop adding NEW keys once there
        // were >200 distinct ones -- could freeze a histogram before its real
        // peak had been seen, if noise filled those keys first.)
        auto& hist = dbcidHist[ uuid ];
        Int_t c = ++hist[ dbcid ];
        nSamples[ uuid ]++;

        // Anchor = the most populated dbcid, same definition as before, but
        // maintained incrementally. It used to be a full std::max_element scan of
        // the map on every hit once past threshold, which is O(keys) per hit.
        Int_t& best = anchorCount[ uuid ];
        if ( c > best ) {
            best = c;
            dbcidAnchor[ uuid ] = dbcid;
        }
    }

    Short_t time( UShort_t uuid, Short_t dbcid ){
        // dbcid is a 12-bit (4096) circular counter -- a channel whose true
        // anchor sits near the 0/4095 edge would otherwise get hits on the
        // far side of the wrap reported as a difference near +-4096 instead
        // of their true small offset. Wrap to the shortest signed distance.
        Short_t diff = dbcid - dbcidAnchor[ uuid ];
        if ( diff > 2048 ) diff -= 4096;
        if ( diff < -2048 ) diff += 4096;
        return diff;
    }

    Short_t anchor( UShort_t uuid ) {
        if ( dbcidAnchor.count( uuid ) == 0 )
            return -1;
        return dbcidAnchor[ uuid ];
    }

    
    size_t samples( UShort_t uuid ) {
        if ( nSamples.count(uuid) > 0 ) return nSamples[uuid];
        return 0;
    }

     const map<Short_t, Int_t>& histFor( UShort_t uuid ) {
        return dbcidHist[ uuid ];
     }

    map< UShort_t, Short_t> &getAnchorMap(){
        return dbcidAnchor;
    }

    void clear(){
        dbcidHist.clear();
        dbcidAnchor.clear();
        nSamples.clear();
        anchorCount.clear();
    }


protected:

    // key - unique VMM Id (0,386]
    // value - histogram (map) with key: deltaBCID, value: counts;
    map< UShort_t, map<Short_t, Int_t> > dbcidHist;

    // key - unique VMM Id (0, 386]
    // value - dbcid reference
    map< UShort_t, Short_t> dbcidAnchor;

    // key - unique VMM Id; value - total hits seen (what ready() tests)
    map< UShort_t, Int_t> nSamples;

    // key - unique VMM Id; value - count in the current anchor bin
    map< UShort_t, Int_t> anchorCount;

};

#endif