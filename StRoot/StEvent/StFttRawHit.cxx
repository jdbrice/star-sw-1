/***************************************************************************
 *
 * $Id: StFttRawHit.cxx,v 2.1 2021/11/18 14:53:48 jdb Exp $
 *
 * Author: jdb, Nov 2021
 ***************************************************************************
 *
 * Description: 
 *
 ***************************************************************************/ 
#include "StFttRawHit.h"
#include <cmath>


StFttRawHit::StFttRawHit() 
: mSector(0),
mRDO(0),
mFEB(0),
mVMM(0),
mChannel(0),
mADC(0),
mBCID(0),
mTB(-32000),
mBCIDDelta(-32000),
mTime(-32000),
mPlane(255),
mQuadrant(kFttUnknownQuadrant),
mRow(255),
mStrip(255),
mOrientation(kFttUnknownOrientation),
mIdTruth(0),
mQaTruth(0),
// These three were in NO constructor, so a raw hit built outside the DAQ path
// -- StFttSlowSimMaker builds its own -- carried whatever was on the stack
// until something called setStripEdges(). -1 is what StFttDb uses for "no
// value", e.g. the edges of an XY strip, which has no edge table.
mStripCenter(-1),
mStripLeftEdge(-1),
mStripRightEdge(-1),
mStripLength(-1)
{ /*noop*/ }

StFttRawHit::StFttRawHit(   UChar_t mSector, UChar_t mRDO, UChar_t mFEB, 
                            UChar_t mVMM, UChar_t mChannel, UShort_t mADC, 
                            UShort_t mBCID, Short_t mTB, Short_t mBCIDDelta){
    setRaw( mSector, mRDO, mFEB, mVMM, mChannel, mADC, mBCID, mTB, mBCIDDelta);
    // This ctor does not delegate, so the strip geometry has to be defaulted
    // here as well -- see the note in the default ctor.
    mPlane       = 255;
    mQuadrant    = kFttUnknownQuadrant;
    mRow         = 255;
    mStrip       = 255;
    mOrientation = kFttUnknownOrientation;
    mIdTruth     = 0;
    mQaTruth     = 0;
    mTime        = -32000;
    mStripCenter = -1;
    mStripLeftEdge  = -1;
    mStripRightEdge = -1;
    mStripLength    = -1;
}

void StFttRawHit::setRaw(   UChar_t mSector, UChar_t mRDO, UChar_t mFEB, 
                            UChar_t mVMM, UChar_t mChannel, UShort_t mADC, 
                            UShort_t mBCID, Short_t mTB, Short_t mBCIDDelta){
    this->mSector    = mSector;
    this->mRDO       = mRDO;
    this->mFEB       = mFEB;
    this->mVMM       = mVMM;
    this->mChannel   = mChannel;
    this->mADC       = mADC;
    this->mBCID      = mBCID;
    this->mTB        = mTB;
    this->mBCIDDelta = mBCIDDelta;
}

void StFttRawHit::setMapping(   UChar_t mPlane, UChar_t mQuadrant, 
                                UChar_t mRow, UChar_t mStrip, UChar_t mOrientation ){
    this->mPlane       = mPlane;
    this->mQuadrant    = mQuadrant;
    this->mRow         = mRow;
    this->mStrip       = mStrip;
    this->mOrientation = mOrientation;
}

void StFttRawHit::setStripEdges( Float_t mStripCenter, Float_t mStripLeftEdge, Float_t mStripRightEdge ) {
    this->mStripCenter      = mStripCenter;
    this->mStripLeftEdge    = mStripLeftEdge;
    this->mStripRightEdge   = mStripRightEdge;
}

void StFttRawHit::setStripLength( Float_t mStripLength ) {
    this->mStripLength = mStripLength;
}

ostream&
operator<<( ostream &os, const StFttRawHit& rh )
{
    os << " StFttRawHit( "    << endl;
    os << "\tmSector = "      << (int)rh.sector()      << endl;
    os << "\tmRDO = "         << (int)rh.rdo()         << endl;
    os << "\tmFEB = "         << (int)rh.feb()         << endl;
    os << "\tmVMM = "         << (int)rh.vmm()         << endl;
    os << "\tmChannel = "     << (int)rh.channel()     << endl;
    os << "\tmADC = "         << (int)rh.adc()         << endl;
    os << "\tmBCID = "        << (int)rh.bcid()        << endl;
    os << "\tmBCIDDelta = "   << (int)rh.dbcid()       << endl;
    os << "\tmTB = "          << (int)rh.tb()          << endl;
    os << "\tmTime = "        << (int)rh.time()        << endl;
    os << "\tmPlane = "       << (int)rh.plane()       << endl;
    os << "\tmQuadrant = "    << (int)rh.quadrant()    << endl;
    os << "\tmRow = "         << (int)rh.row()         << endl;
    os << "\tmStrip = "       << (int)rh.strip()       << endl;
    os << "\tmOrientation = " << (int)rh.orientation() << endl;
    os << "\tmStripCenter = " << rh.stripCenter()      << endl;
    os << "\tmStripLength = " << rh.stripLength()      << endl;
    os << "\tmStripLeftEdge = " << rh.stripLeftEdge()  << endl;
    os << "\tmStripRightEdge = " << rh.stripRightEdge()<< endl;
    os << "\tidTruth = "      << (int)rh.idTruth()     << endl;
    os << "\tqaTruth = "      << (int)rh.qaTruth()     << endl;
    os << " ) " << endl;
    return os;
}