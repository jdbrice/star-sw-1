/***************************************************************************
 * StFttDb.h
 * jdb Feb, 2022
 ***************************************************************************
 *
 * Description: FTT DB Utility
 ***************************************************************************/

#ifndef STFTTDB_H
#define STFTTDB_H

#ifndef StMaker_H
#include "StMaker.h"
#endif
#include "StEvent/StEnumerations.h"
#include <stdint.h>
#include "TGeoMatrix.h"
#include <map>
#include <vector>


class StFttRawHit;
class StFttCluster;
class StFttPoint;

struct FttDataWindow {
    Short_t uuid;
    UChar_t mode;
    Short_t min;
    Short_t max;
    Short_t anchor;
};


class St_fttHardwareMap;
class St_fttDataWindowsB;

class St_Survey;

class StFttDb : public TDataSet {

public: 
  StFttDb(const char *name="fttDb");
  virtual       ~StFttDb();
  int  Init();
  int  InitRun(int runNumber);
 
  void SetDebug(int v=1) {mDebug=v;}  //! debug level
  void setDebug(int v=1) {mDebug=v;}  //! debug level
  void setDbAccess(int v=1);  //! enable(1) or disable(0) offline DB access
  void setRun(int run);       //! set run# 

  static size_t uuid( StFttRawHit * h, bool includeStrip = false ) ;
  static size_t uuid( StFttCluster * c ) ;
  static size_t vmmId( StFttRawHit * h ) ;

  // HARDWARE Mapping StFttRawHits
    uint16_t packKey( int feb, int vmm, int ch ) const;
    void unpackKey( int key, int &feb, int &vmm, int &ch ) const;
    uint16_t packVal( int row, int strip ) const;
    void unpackVal( int val, int &row, int &strip ) const;
    void loadHardwareMapFromFile( std::string fn );
    bool loadStripCenterFromFile( std::string fn );
    bool loadStripEdgeFromFile( std::string fn );
    bool loadStripLengthFromFile( std::string fn );
    void loadHardwareMapFromDb( St_fttHardwareMap * );
    void loadDataWindowsFromFile( std::string fn );
    void loadDataWindowsFromDb( St_fttDataWindowsB * );


    UChar_t plane( StFttRawHit * hit );
    UChar_t quadrant( StFttRawHit * hit );
    UChar_t fob( StFttRawHit * hit );
    UChar_t rob( StFttRawHit * hit );
    UChar_t rob( StFttCluster * clu );

    static double stripPitch; // mm
    static double gapPitch; // mm
    static double stripWidth; // mm
    static double rowLength; // mm
    static double HVStripShift; // mm
    static double DiagStripShift; // mm
    static double lowerQuadOffsetX; //mm
    static double idealPlaneZLocations[4];

    static const size_t nPlane        = 4;
    static const size_t nQuadPerPlane = 4;
    static const size_t nFobPerQuad   = 6;
    static const size_t nVMMPerFob    = 4;
    static const size_t nChPerVMM     = 64;
    static const size_t nStripGroupEdge = 8;

    //name for the cluster direction
    static TString Direction_name[nQuadPerPlane+1];
    static double FirstStripEdge[2]; //mm

    //for idealPlaneZLocations_QuadX, now using the cm as unit because that in the old version this using the cm as unit
    static double LocalStripZLocations[nPlane];
    static double idealPlaneZLocations_QuadA[nPlane];//cm
    static double idealPlaneZLocations_QuadB[nPlane];
    static double idealPlaneZLocations_QuadC[nPlane];
    static double idealPlaneZLocations_QuadD[nPlane];
    static double X_shift_QuadA[nQuadPerPlane];//mm , x shift from pin hole to (0,0)
    static double X_shift_QuadB[nQuadPerPlane];//mm , x shift from pin hole to (0,0)
    static double X_shift_QuadC[nQuadPerPlane];//mm , x shift from pin hole to (0,0)
    static double X_shift_QuadD[nQuadPerPlane];//mm , x shift from pin hole to (0,0)
    static double Y_shift_QuadA[nQuadPerPlane];//mm , y shift from pin hole to (0,0)
    static double Y_shift_QuadB[nQuadPerPlane];//mm , y shift from pin hole to (0,0)
    static double Y_shift_QuadC[nQuadPerPlane];//mm , y shift from pin hole to (0,0)
    static double Y_shift_QuadD[nQuadPerPlane];//mm , y shift from pin hole to (0,0)
    static double YX_StripGroupEdge[nStripGroupEdge]; // mm , from left to right, reference : https://drupal.star.bnl.gov/STAR/system/files/StFttPointMaker_ZhenWang_20221013_fwdsoftmeeting.pdf Slide 14
    static double D_StripGroupEdge[nStripGroupEdge]; // mm , from left to right, reference : https://drupal.star.bnl.gov/STAR/system/files/StFttPointMaker_ZhenWang_20221013_fwdsoftmeeting.pdf Slide 14
    static double X_StripGroupEdge[nStripGroupEdge]; // mm , from left to right, reference : https://drupal.star.bnl.gov/STAR/system/files/
    static double Y_StripGroupEdge[nStripGroupEdge]; // mm , from left to right, reference : https://drupal.star.bnl.gov/STAR/system/files/StFttPointMaker_ZhenWang_20221013_fwdsoftmeeting.pdf Slide 14

    static const size_t nQuad = nQuadPerPlane * nPlane; // 4 * 4 = 16 Total number of Quadrants
    static const size_t nRob = nQuad;

    static const size_t nFobPerPlane = nFobPerQuad * nQuadPerPlane; // 4 * 6 = 24 Total number of FOB
    static const size_t nFob         = nFobPerPlane * nPlane; // 24 * 4 = 96 total FOB;

    static const size_t nVMMPerQuad  = nFobPerQuad * nVMMPerFob; // 6 * 4 = 24 VMM per Quad
    static const size_t nVMMPerPlane = nVMMPerQuad * nQuadPerPlane; // 24 * 4 = 96 VMM per Plane
    static const size_t nVMM         = nVMMPerPlane * nPlane; // 96 * 4 = 384 VMM total

    static const size_t nChPerFob    = nVMMPerFob * nChPerVMM; // 4 * 64 = 256 VMM per Plane
    static const size_t nChPerQuad   = nChPerFob * nFobPerQuad; // 256 * 6 = 1536 VMM per Plane
    static const size_t nChPerPlane  = nChPerQuad * nQuadPerPlane; // 1536 * 4 = 6144 Ch per Plane
    static const size_t nCh          = nChPerPlane * nPlane; // 6144 * 4 = 24576 Ch per Plane

    static const size_t maxADC       = 1024 + 1; // 10 bits;
    static const size_t maxBCID      = 4096 + 1; // 12 bits;
    static const int minTb           = -32768 - 1000; // get the under and overflow
    static const int maxTb           = 32768 + 1000;
    
    static const size_t nRowsPerQuad       = 5; // 0 - 2 (Horizontal or Vertical) + 3, 4 for Diagonal
    static const size_t nStripOrientations = 5; // real orientations (4) + 1 for unkown                       = 5 total
    static const size_t maxStripPerRow     = 166; // max number of strips in a row (only row 0 has this many)
    
    static vector<string> orientationLabels;



    UChar_t orientation( StFttRawHit * hit );
    UChar_t getOrientation( int rob, int feb, int vmm, int row ) const;
    bool hardwareMap( int rob, int feb, int vmm, int ch, int &row, int &strip, UChar_t &orientation ) const;
    bool hardwareMap( StFttRawHit * rawHit ) const;

    // ---- geometry offsets: from the DB (default) or the original hardcoded numbers ----
    // The X_shift_Quad* / Y_shift_Quad* / idealPlaneZLocations_Quad* statics above are
    // the ORIGINAL hardcoded numbers. They stay as the placeholder, and are exactly what
    // setUseDbGeometry(false) restores, so the old behaviour is always one call away.
    //
    // With the DB on (the default) the hardcoded origins above are NOT used at all.
    // The effective offset is built from the AGML nominal instead:
    //     AGML placeholder + stgcOnTpc + stationOnStgc[station] + pentOnStation[4*station+pent]
    // The tables were written as (StFttDb origin - AGML placeholder) + alignment, so this
    // reproduces the true position exactly while being the SAME quantity AGML computes.
    // Geometry and hit positions then come from one source rather than two kept in step
    // by hand, and an alignment iteration only has to touch the table.
    // with pent = (4 - quad) % 4: AGML orders the four pentagons A,D,C,B by rotation
    // while StFttDb numbers its quadrants A,B,C,D. That map was measured from the built
    // geometry, see FwdGeomUtils::fttChamberIndex and script/fttQuadConvention.C.
    //
    // Units: Survey translations are cm. X/Y_shift are mm, so t0/t1 are scaled by 10;
    // t2 adds straight to the z locations, which are already cm.
    //
    // Rows are taken in table order, which is what AGML's <Misalign row="N"/> indexes.
    // ---- survey-table index conventions ---------------------------------
    // AGML orders the four pentagons of a station A,D,C,B by rotation about z,
    // StFttDb numbers its quadrants A,B,C,D. Measured from the built geometry,
    // see script/fttQuadConvention.C. One copy of that permutation lives in
    // StFttDb.cxx and these are the only ways to reach it.
    //
    // FTT survey rows are POSITIONAL -- row = 4*station + pent, which is what
    // AGML's <Misalign row="N"/> indexes, and what loadGeometryFromDb() reads.
    // (The FST is the other way round: its rows are addressed by the Id column,
    // see StFstDb::sensorSurveyId. A writer that confuses the two produces a
    // shuffled table that still loads.)
    //
    // WRITING a table needs surveyRow(): a correction for (station, quad) goes
    // in that row. The reverse pair is for reading a table back and saying
    // which quadrant a row belongs to.
    static int  quadToPent( int quad );     // StFttDb A,B,C,D -> AGML pent index
    static int  pentToQuad( int pent );     // and back
    static int  surveyRow( int station, int quad );            // -> 0..15, or -1
    static bool surveyRowToStationQuad( int row, int &station, int &quad );

    void setUseDbGeometry( bool v ) { mUseDbGeometry = v; }
    bool useDbGeometry() const      { return mUseDbGeometry; }
    bool geometryFromDb() const     { return mGeoFromDb; }   // a DB load actually happened
    void resetGeometryToHardcoded();
    int  loadGeometryFromDb( St_Survey *stgcOnTpc, St_Survey *stationOnStgc, St_Survey *pentOnStation );
    // The values getGloablOffset*() actually uses. Prefer these over the statics, which
    // are only the placeholder once the DB is in play.
    double xShift( int quad, int plane ) const { return mXShift[quad][plane]; }  // mm
    double yShift( int quad, int plane ) const { return mYShift[quad][plane]; }  // mm
    double zLoc  ( int quad, int plane ) const { return mZLoc  [quad][plane]; }  // cm

    void getGloablOffset( UChar_t plane, UChar_t quad, float &dx, float &sx, float &dy, float &sy, float &dz, float &sz );
    void getGloablOffset_ClusterPoint( UChar_t plane, UChar_t quad, float &dx, float &sx, float &dy, float &sy, float &dz, float &sz );

    // ---- quadrant local <-> global, as a real transform ------------------------
    // getGloablOffset_ClusterPoint() returns (d,s) pairs with s = +-1, i.e. a
    // translation and axis REFLECTIONS. AGML rotates the four pentagons instead:
    // with quad->pent {0,3,2,1}, quadrant A gets the identity and C gets 180
    // degrees -- those two happen to agree with the reflections -- while B gets
    // -90 and D gets +90, which a single-axis reflection cannot reproduce. That
    // is why StFttDb and the built geometry agree exactly in quadrants A and C
    // and differ by ~140-200 mm in B and D.
    //
    // These matrices carry the rotation properly and compose the survey levels as
    // matrix PRODUCTS rather than by adding translations, so a survey row with a
    // real rotation, or a new level such as TpcPosition, is one more factor
    // instead of a re-derivation. Units are cm, as TGeo uses.
    const TGeoHMatrix& quadToGlobal( int quad, int plane ) const;
    bool localToGlobal( UChar_t plane, UChar_t quad, const double l[3], double g[3] ) const;
    bool globalToLocal( UChar_t plane, UChar_t quad, const double g[3], double l[3] ) const;
    // Rebuild the matrices from the current offsets. Called by the loaders; public
    // so a test can force it after poking the shifts.
    void buildQuadMatrices();

    // ---- AGML quadrant-local <-> StFttDb quadrant-local ------------------------
    // The simulation reports a hit in AGML's quadrant-local frame (g2t_fts_hit's
    // xlocal, which is STMG-local: STMG sits at x=0,y=0 inside STFM, so its x,y
    // ARE the pentagon's). StFttDb's local frame is a different one -- different
    // origin AND different handedness, because AGML rotates the four pentagons
    // while StFttDb reflects them. Neither is wrong; they are two conventions,
    // and this is the bridge between them.
    //
    // Composing StFttDb's reflection S = diag(sx,sy) with AGML's pentagon
    // rotation R gives T = S^-1 R, which works out to
    //     quad A  S=diag( 1, 1)  R=identity   ->  T = identity
    //     quad B  S=diag( 1,-1)  R=rot(-90)   ->  T = (x,y) -> (y,x)
    //     quad C  S=diag(-1,-1)  R=rot(180)   ->  T = identity
    //     quad D  S=diag(-1, 1)  R=rot(+90)   ->  T = (x,y) -> (y,x)
    // i.e. a pure frame rotation with NO translation -- measured as exactly 0 in
    // all 16 (quadrant, station) at sdt20211112.000001, and that is only true
    // when both sides were built from the same survey tables. In hardcoded mode
    // the two origins are 8-48 mm apart, so this REQUIRES geometryFromDb() and
    // refuses otherwise rather than returning a quietly wrong point.
    //
    // The map is its own inverse, so the same call converts either way; the
    // second name exists only so the calling code reads correctly.
    // Units are whatever the caller uses -- this rotates, it does not rescale,
    // so do not mix the cm of quadToGlobal() with the mm of getGloablOffset().
    bool agmlQuadLocalToDbLocal( UChar_t plane, UChar_t quad, const double a[3], double d[3] ) const;
    bool dbLocalToAgmlQuadLocal( UChar_t plane, UChar_t quad, const double d[3], double a[3] ) const
        { return agmlQuadLocalToDbLocal( plane, quad, d, a ); }
    // ---- strip geometry: where a strip physically is, in quadrant-local mm ----
    // The four orientations are two strip layers on each of two foils, running at
    // 45 degrees to each other, so an X strip and a U strip either cross inside the
    // sensitive area or they cannot have come from the same particle. That makes a
    // geometric gate possible on any pairing of clusters, which is what
    // StFttClusterPointMaker has to decide and what a charge-correlation study has
    // to control for.
    //
    // stripSegment() returns the strip's two endpoints. The precision coordinate is
    // (strip-0.5)*stripPitch, matching StFttSlowSimMaker::computeRowStrip(); the
    // along-strip extent comes from the Row*_StripLength tables already loaded here.
    // Rows 0-2 are Vertical/Horizontal, rows 3-4 Diagonal, in the rotated frame
    //     x' = (x + y - 2a)/sqrt2,  y' = (y - x)/sqrt2,  a = D_StripGroupEdge[0]
    // with the sign of y' set by row AND orientation -- kFttDiagonalV runs positive
    // on row 4, kFttDiagonalH on row 3 (see StFttClusterPointMaker::MakeLocalPoints).
    // The precision coordinate of one strip, in quadrant-local mm: the single
    // place that answers "where is strip N". Everything else should ask here
    // rather than recompute, which is the whole point -- the formula
    // (strip-0.5)*stripPitch that StFttClusterMaker.cxx:363 and
    // StFttSlowSimMaker::computeRowStrip() both use is this value minus a
    // CONSTANT 17.5500 mm for XY strips and 26.7400 mm for diagonals
    // (script/fttStripCenterCheck.C measures both over every strip, spread
    // ~1e-13, so the maps carry no per-strip information beyond the offset).
    //
    // Which is right is settled empirically, not geometrically: assuming the
    // FST->sTGC fit absorbed 17.55*(sx,sy), the residual against
    // pentOnStation.20211015 is mean -0.07 mm, rms 1.97, over all 16
    // quadrant/stations. The offset is real and our "alignment" is largely it.
    // There is NO first-principles derivation of 17.55 in the repository; it
    // would need the chamber drawing.
    //
    // Note the two offsets DIFFER, so a single added constant cannot fix both
    // XY and diagonal strips -- which is why this is a lookup and not a number.
    bool stripCenter( int row, int strip, UChar_t orientation, double &center ) const;

    // Everything hardwareMap() puts on a raw hit about where its strip is, in
    // one call: centre, the two edges (diagonals only -- XY strips have no edge
    // table and keep -1, as they always did) and the length. Quadrant-local mm.
    //
    // hardwareMap() itself now goes through this, so a maker that builds raw
    // hits outside the DAQ path -- StFttSlowSimMaker does -- can fill the same
    // fields the same way instead of leaving them unset. They are NOT
    // initialised by StFttRawHit's constructors, so "unset" means whatever was
    // on the stack.
    bool stripEdges( int row, int strip, UChar_t orientation,
                     float &center, float &leftEdge, float &rightEdge, float &length ) const;

    // The inverse: which strip does a quadrant-local precision coordinate fall
    // on? Uses the same maps, so it round trips with stripCenter() by
    // construction rather than by two formulas agreeing.
    bool stripFromLocal( int row, UChar_t orientation, double precisionCoord,
                         int &strip, double &offsetInStrip ) const;

    bool stripSegment( int row, int strip, UChar_t orientation,
                       double &x1, double &y1, double &x2, double &y2 ) const;

    // True if the two strips cross within both of their lengths. tolerance (mm)
    // widens each segment at both ends, for the strip width and for the fact that
    // a cluster's position is not exactly its max-ADC strip centre.
    // crossX/crossY get the crossing point when there is one.
    bool stripsCross( int rowA, int stripA, UChar_t orientA,
                      int rowB, int stripB, UChar_t orientB,
                      double &crossX, double &crossY, double tolerance = 0.0 ) const;

    bool reverseHardwareMap(int &rob, int &feb, int &vmm, int &ch, int plane, int quad, int row, int strip, UChar_t &orientation) const;
    bool reverseHardwareMap( int &feb, int &vmm, int &ch, int row, int strip ) const;

    enum TimeCutMode {
      CalibratedBunchCrossingMode = 0,
      TimebinMode = 1
    };

    void setTimeCut(TimeCutMode mode, int low=-70, int high=30)
    {
        mUserDefinedTimeCut = true;
        mTimeCutMode=mode; 
        mTimeCutLow=low; 
        mTimeCutHigh=high;
    }

    void getTimeCut( StFttRawHit * hit, int &mode, int &l, int &h );

    // ---- per-run VMM time anchors (Calibrations/ftt/fttDataWindowsB) ------------
    // VMM dbcid anchors are reset at random at every run start, so an anchor is only
    // meaningful for the run it was measured in. STAR DB entries carry a begin time
    // only, so St_db_Maker hands a run with no entry of its own the PREVIOUS run's
    // entry. Hence the RUN STAMP written by script/fttDataWindow_db.C into the spare
    // slot 384 (uuid -1, min = run/10000, max = run%10000), and anchors are used only
    // when StFttDbMaker has confirmed stamp == run being reconstructed.
    int  dataWindowStampRun() const            { return mDwStampRun; }
    void setDataWindowAnchorsValid( bool v )   { mDwAnchorsValid = v; }
    bool dataWindowAnchorsValid() const        { return mDwAnchorsValid; }
    // true, with the anchor, only if anchors are validated for this run AND this VMM
    // has one (anchor -1 in the table = no usable calibration for that VMM)
    bool getAnchor( StFttRawHit * hit, Short_t &anchor );
    // Forget any previously loaded windows/anchors/stamp. Called before every load so a
    // run whose fetch fails cannot silently keep the previous run's entry.
    void clearDataWindows() { dwMap.clear(); mDwStampRun = -1; mDwAnchorsValid = false; mDwFallbackWarned = false; }

 private:
  bool   mUseDbGeometry = true;          //! true = placeholder + DB survey tables (default)
  bool   mGeoFromDb     = false;         //! set once a DB load has actually been applied
  double mXShift[nQuadPerPlane][nPlane]; //! mm, effective per quadrant/plane
  double mYShift[nQuadPerPlane][nPlane]; //! mm
  double mZLoc  [nQuadPerPlane][nPlane];
    TGeoHMatrix mQuadToGlobal[nQuadPerPlane][nPlane];   //! quadrant local -> global //! cm

  int   mDbAccess=1;                     //! enable(1) or disabe(0) DB access
  int   mRun=0;                          //! run#
  // int   mDebug=1;                     //! >0 dump tables to text files    
  int   mDebug=0;                        //! >0 dump tables to text files    

  bool mUserDefinedTimeCut = false;
  // Returned by getTimeCut() for a VMM with no data-window entry. These used to have
  // no initializer, so such a VMM got an undefined mode and window -- and since
  // StFttClusterMaker defaults to kTimeCutModeDB, that reaches every chain that does
  // not call SetTimeCut() (BFC production, simulation). A VMM is "missing" when no
  // fttDataWindowsB entry was loaded at all (query time before 2021-10-01, no DB
  // access), when an entry's uuid column skips ids, or when the hit's vmmId is out of
  // range. Default = calibrated time, [-65, 100]: exactly the global DB entry in use
  // since 2021-10-01, so a missing entry behaves like the DB default.
  TimeCutMode mTimeCutMode = CalibratedBunchCrossingMode;
  int mTimeCutLow  = -65;
  int mTimeCutHigh = 100;
  bool mDwFallbackWarned = false;  //! warn once per load, not once per hit

  std :: map< uint16_t , uint16_t > mMap;
  std :: map< uint16_t , std::vector<uint16_t> > rMap; // reverse map -- a given (row,strip) has one channel per orientation (H and V both exist at the same row/strip in the real hardware map), so this must hold all of them, not just the last one loaded
  //  data windows map
  std :: map< uint16_t , FttDataWindow > dwMap;
  int  mDwStampRun = -1;          //! run stamped into the loaded data-window entry, -1 = unstamped
  bool mDwAnchorsValid = false;   //! set by StFttDbMaker only when the stamp matches the run
  std :: map< int , Float_t > scMapXY; // strip center map 
  std :: map< int , Float_t > scMapDiag; // strip center map 
  std :: map< int , Float_t > slMapRow1; // strip length map Row1
  std :: map< int , Float_t > slMapRow2; // strip length map Row2
  std :: map< int , Float_t > slMapRow3; // strip length map Row3
  std :: map< int , Float_t > slMapRow4; // strip length map Row4
  std :: map< int , Float_t > slMapRow5; // strip length map Row5
  std :: map< int , Float_t > seMapDiagLeft; // strip left edge map, just for the diagonal strips 
  std :: map< int , Float_t > seMapDiagRight; // strip right edge map, just for the diagonal strips

  ClassDef(StFttDb,0)   //StAF chain virtual base class for Makers        
};

#endif

