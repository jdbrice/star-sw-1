/***************************************************************************
 * StFttDb.cxx
 * jdb & Zhen Feb, 2022
 ***************************************************************************
 * Description: This interface between FTT and the STAR database
 ***************************************************************************/

#include "StFttDb.h"
#include "StMaker.h"
#include "StMessMgr.h"
#include "StEvent/StFttRawHit.h"
#include "StEvent/StFttCluster.h"
#include "StEvent/StFttPoint.h"
#include <math.h>

#include "tables/St_fttHardwareMap_Table.h"
#include "tables/St_fttDataWindowsB_Table.h"
#include "tables/St_Survey_Table.h"


ClassImp(StFttDb)


TString StFttDb::Direction_name[] = {"kFttHorizontal","kFttVertical","kFttDiagonalH","kFttDiagonalV","kFttUnknownOrientation"};

double StFttDb::stripPitch = 3.2; // mm
double StFttDb::gapPitch = 0.5; // mm
double StFttDb::stripWidth = 2.7; // mm
double StFttDb::rowLength = 180; // mm
double StFttDb::lowerQuadOffsetX = 101.6; // mm
// double StFttDb::idealPlaneZLocations[] = { 281.082,304.062,325.058,348.068 };//ideal position
double StFttDb::idealPlaneZLocations[] = { 312.342,329.953,347.637,365.422 };//suvery data, cm or mm? now just use quad A's data
double StFttDb::LocalStripZLocations[] = { 1.1    ,1.53   ,2.18   ,2.61    };// from yingying's measurement, cm ,
double StFttDb::idealPlaneZLocations_QuadA[] = {312.342,329.953,347.637,365.422};//cm
double StFttDb::idealPlaneZLocations_QuadB[] = {312.165,330.094,347.543,365.409};//cm
double StFttDb::idealPlaneZLocations_QuadC[] = {312.317,329.723,347.346,365.311};//cm
double StFttDb::idealPlaneZLocations_QuadD[] = {312.098,329.731,347.444,365.455};//cm
 
double StFttDb::HVStripShift = 15.95;//mm
double StFttDb::DiagStripShift = 19.42;//mm
double StFttDb::FirstStripEdge[] = {14.6, 19.42};
vector<string> StFttDb::orientationLabels = { "Horizontal", "Vertical", "DiagonalH", "DiagonalV", "Unknown" };
double StFttDb::X_shift_QuadA[] = {8.09, 8.34, 6.62,7.54 };//mm 
double StFttDb::X_shift_QuadB[] = {112.74, 112.14, 113.30, 113.49};//mm 
double StFttDb::X_shift_QuadC[] = {-107.51, -108.22, -109.87, -108.91};//mm 
double StFttDb::X_shift_QuadD[] = {-3.69, -4.96, -4.36, -3.75};//mm 
double StFttDb::Y_shift_QuadA[] = {95.34, 94.33, 96.03, 95.01};//mm 
double StFttDb::Y_shift_QuadB[] = {84.24, 83.37, 83.61, 83.81};//mm
double StFttDb::Y_shift_QuadC[] = {83.60, 83.42, 84.37, 82.81};//mm
double StFttDb::Y_shift_QuadD[] = {95.70, 94.40, 95.55, 94.16};//mm
double StFttDb::YX_StripGroupEdge[] = {11.49, 172.29, 360.09};//mm 
double StFttDb::D_StripGroupEdge[] = {11.49};//mm 
double StFttDb::X_StripGroupEdge[] = {14.60, 172.29, 216.89, 315.4, 360.09, 410.9, 504.2, 548.7};//mm 
double StFttDb::Y_StripGroupEdge[] = {14.60, 172.29, 216.89, 315.4, 360.09, 410.9, 504.2, 548.7};//mm 


StFttDb::StFttDb(const char *name) : TDataSet(name) { resetGeometryToHardcoded(); }; 

StFttDb::~StFttDb() {}

int StFttDb::Init(){

  return kStOK;
}

void StFttDb::setDbAccess(int v) {mDbAccess =  v;}
void StFttDb::setRun(int run) {mRun = run;}

int StFttDb::InitRun(int runNumber) {
    mRun=runNumber;
    return kStOK;
}


size_t StFttDb::uuid( StFttRawHit * h, bool includeStrip ) {
    // this UUID is not really universally unique
    // it is unique up to the hardware location 
    // at the give precision
    // NOT including strip level is useful for cluster 
    // calculations that combine all strips from given 
    // plane, quad, row, orientation

    
    size_t _uuid = (size_t)h->orientation() + (nStripOrientations) * ( h->row() + nRowsPerQuad * ( h->quadrant() + nQuadPerPlane * h->plane() ) );
    
    if ( includeStrip ){
        _uuid = (size_t) h->strip() * maxStripPerRow *( h->orientation() + (nStripOrientations) * ( h->row() + nRowsPerQuad * ( h->quadrant() + nQuadPerPlane * h->plane() ) ) );
    } 

    return _uuid;
}

size_t StFttDb::uuid( StFttCluster * c ) {
    // this UUID is not really universally unique
    // it is unique up to the hardware location

    size_t _uuid = (size_t)c->orientation() + (nStripOrientations) * ( c->row() + nRowsPerQuad * ( c->quadrant() + nQuadPerPlane * c->plane() ) );
    return _uuid;
}

size_t StFttDb::vmmId( StFttRawHit * h ) {
    // Calculate VMM hardware ID based on electronic readout structure
    // VMM_ID = vmm + nVMMPerFob * (feb + nFobPerQuad * (quadrant + nQuadPerPlane * plane))
    // Where: plane [0-3], quadrant [0-3], feb [0-5], vmm [0-3]
    // Valid range: 0-383 (total of 384 VMMs)

    u_char iPlane = h->sector() - 1;     // sector is 1-based
    u_char iQuad  = h->rdo() - 1;        // rdo is 1-based
    u_char iFeb   = h->feb();            // feb is 0-based
    u_char iVmm   = h->vmm();            // vmm is 0-based

    size_t vmm_id = iVmm + nVMMPerFob * ( iFeb + nFobPerQuad * ( iQuad + nQuadPerPlane * iPlane ) );

    return vmm_id;
}


void StFttDb::getTimeCut( StFttRawHit * hit, int &mode, int &l, int &h ){
        mode = mTimeCutMode;
        l = mTimeCutLow;
        h = mTimeCutHigh;
        if (mUserDefinedTimeCut)
            return;

        // load calibrated data windows from DB
        // NOTE: dwMap is indexed by VMM hardware ID, not geometric UUID
        size_t hit_vmmid = vmmId( hit );

        // Validate VMM ID is in expected range
        if ( hit_vmmid >= nVMM ) {
            LOG_WARN << "StFttDb::getTimeCut - VMM ID out of range: " << hit_vmmid
                     << " (max=" << (nVMM-1) << ")" << endm;
            LOG_WARN << "  Hit: plane=" << (int)plane(hit)
                     << " quad=" << (int)quadrant(hit)
                     << " feb=" << (int)hit->feb()
                     << " vmm=" << (int)hit->vmm() << endm;
            return;
        }

        if ( dwMap.count( hit_vmmid ) ){
            mode = dwMap[ hit_vmmid ].mode;
            l = dwMap[ hit_vmmid ].min;
            h = dwMap[ hit_vmmid ].max;
        } else if ( !mDwFallbackWarned ) {
            LOG_WARN << "StFttDb::getTimeCut - no data-window entry for VMM " << hit_vmmid
                     << (dwMap.empty() ? " (no fttDataWindowsB loaded at all)" : "")
                     << "; using default mode " << (int)mode << " window [" << l << ", " << h
                     << "]. Warned once; other VMMs may also fall back." << endm;
            mDwFallbackWarned = true;
        }

    }


uint16_t StFttDb::packKey( int feb, int vmm, int ch ) const{
    // feb = [1 - 6] = 3 bits
    // vmm = [1 - 4] = 3 bits
    // ch  = [1 - 64] = 7 bits
    return feb + (vmm << 3) + (ch << 6);
}
void StFttDb::unpackKey( int key, int &feb, int &vmm, int &ch ) const{
    feb = key & 0b111;
    vmm = (key >> 3) & 0b111;
    ch  = (key >> 6) & 0b1111111;
    return;
}
uint16_t StFttDb::packVal( int row, int strip ) const{
    // row = [1 - 4] = 3 bits
    // strip = [1 - 152] = 8 bits
    return row + ( strip << 3 );
}
void StFttDb::unpackVal( int val, int &row, int &strip ) const{
    row = val & 0b111; // 3 bits
    strip = (val >> 3) & 0b11111111; // 8 bit
    return;
}

void StFttDb::loadDataWindowsFromDb( St_fttDataWindowsB * dataset ) {
    if (dataset) {
        Int_t rows = dataset->GetNRows();

        if ( !rows ) return;

        dwMap.clear();
        mDwStampRun = -1;
        mDwAnchorsValid = false;

        fttDataWindowsB_st *table = dataset->GetTable();
        const int nEntries = sizeof(table[0].uuid) / sizeof(table[0].uuid[0]);
        for (Int_t i = 0; i < rows; i++) {
            for ( int j = 0; j < nEntries; j++ ) {
                // Run stamp (see StFttDb.h): not a VMM, never goes into dwMap.
                if ( table[i].uuid[j] == -1 ) {
                    mDwStampRun = table[i].min[j] * 10000 + table[i].max[j];
                    continue;
                }
                if ( table[i].uuid[j] < 0 ) continue;
                // printf( "[feb=%d, vmm=%d, ch=%d] ==> [row=%d, strip%d]\n", table[i].feb[j], table[i].vmm[j], table[i].vmm_ch[j], table[i].row[j], table[i].strip[j] );


                // uint16_t key = packKey( table[i].feb[j], table[i].vmm[j], table[i].vmm_ch[j] );
                // uint16_t val = packVal( table[i].row[j], table[i].strip[j] );
                // mMap[ key ] = val;
                // rMap[ val ] = key;
                FttDataWindow fdw;
                fdw.uuid   = table[i].uuid[j];
                fdw.mode   = table[i].mode[j];
                fdw.min    = table[i].min[j];
                fdw.max    = table[i].max[j];
                fdw.anchor = table[i].anchor[j];
                dwMap[ fdw.uuid ] = fdw;

                // std::cout << (int)table[i].feb[j] << std::endl;
            }
            // sample output of first member variable
        }
    } else {
        LOG_ERROR << "dataset does not contain requested table" << endm;
    }
}

// Reads fttDataWindow.<run>.txt as written by script/fttDataWindow.C: '#' lines are
// comments (the "# run N" header supplies the run stamp); data lines start with
// "uuid mode min max anchor" and any further QA columns are ignored. Same content
// and conventions as the DB table, so the same run-stamp check applies.
void StFttDb::loadDataWindowsFromFile( std::string fn ) {
    std::ifstream inf( fn.c_str() );
    if ( !inf.good() ) {
        LOG_ERROR << "StFttDb::loadDataWindowsFromFile - cannot open " << fn << endm;
        return;
    }
    dwMap.clear();
    mDwStampRun = -1;
    mDwAnchorsValid = false;

    std::string line;
    int nRead = 0;
    while ( std::getline( inf, line ) ) {
        if ( line.empty() ) continue;
        if ( line[0] == '#' ) {
            int r = 0;
            if ( sscanf( line.c_str(), "# run %d", &r ) == 1 ) mDwStampRun = r;
            continue;
        }
        int u, m, lo, hi, an;
        if ( sscanf( line.c_str(), "%d %d %d %d %d", &u, &m, &lo, &hi, &an ) != 5 ) continue;
        if ( u < 0 || u >= (int)nVMM ) continue;
        FttDataWindow fdw;
        fdw.uuid   = u;
        fdw.mode   = m;
        fdw.min    = lo;
        fdw.max    = hi;
        fdw.anchor = an;
        dwMap[ fdw.uuid ] = fdw;
        nRead++;
    }
    LOG_INFO << "StFttDb::loadDataWindowsFromFile - " << nRead << " VMM entries from " << fn
             << ", run stamp " << mDwStampRun << endm;
}

bool StFttDb::getAnchor( StFttRawHit * hit, Short_t &anchor ) {
    if ( !mDwAnchorsValid ) return false;
    size_t id = vmmId( hit );
    if ( id >= nVMM ) return false;
    std::map< uint16_t, FttDataWindow >::const_iterator it = dwMap.find( (uint16_t)id );
    if ( it == dwMap.end() ) return false;
    if ( it->second.anchor < 0 ) return false;
    anchor = it->second.anchor;
    return true;
}


void StFttDb::loadHardwareMapFromDb( St_fttHardwareMap * dataset ) {
    if (dataset) {
        Int_t rows = dataset->GetNRows();
        if (rows > 1) {
            std::cout << "INFO: found INDEXED table with " << rows << " rows" << std::endl;
        }

        mMap.clear();
        rMap.clear();

        fttHardwareMap_st *table = dataset->GetTable();
        const int nEntries = sizeof(table[0].feb) / sizeof(table[0].feb[0]);
        for (Int_t i = 0; i < rows; i++) {
            for ( int j = 0; j < nEntries; j++ ) {
                uint16_t key = packKey( table[i].feb[j], table[i].vmm[j], table[i].vmm_ch[j] );
                uint16_t val = packVal( table[i].row[j], table[i].strip[j] );
                mMap[ key ] = val;
                rMap[ val ].push_back( key ); // a given (row,strip) has one channel per orientation -- keep all of them
            }
            // sample output of first member variable
        }
    } else {
        std::cout << "ERROR: dataset does not contain requested table" << std::endl;
    }
}


void StFttDb::loadHardwareMapFromFile( std::string fn ){
    std::ifstream inf;
    inf.open( fn.c_str() );

    mMap.clear();
    if ( !inf ) {
        LOG_WARN << "sTGC Hardware map file not found" << endm;
        return;
    }

    rMap.clear();
    string hs0, hs1, hs2, hs3, hs4;
    // HEADER:
    // Row_num    FEB_num    VMM_num    VMM_ch         strip_ch
    inf >> hs0 >> hs1 >> hs2 >> hs3 >> hs4;
    
    if ( mDebug ){
        printf( "Map Header: %s, %s, %s, %s, %s", hs0.c_str(), hs1.c_str(), hs2.c_str(), hs3.c_str(), hs4.c_str() );
        puts("");
    }
    
    uint16_t row, feb, vmm, ch, strip;
    while( inf >> row >> feb >> vmm >> ch >> strip ){
        // pack the key (feb, vmm, ch)
        uint16_t key = packKey( feb, vmm, ch );
        uint16_t val = packVal( row, strip );
        mMap[ key ] = val;
        rMap[ val ].push_back( key ); // a given (row,strip) has one channel per orientation -- keep all of them
        if ( mDebug ){
            printf( "key=%d", key );
            printf( "in=(feb=%d, vmm=%d, ch=%d)\n", feb, vmm, ch );
            int ufeb, uvmm, uch;
            unpackKey( key, ufeb, uvmm, uch );
            printf( "key=(feb=%d, vmm=%d, ch=%d)\n", ufeb, uvmm, uch ); puts("");
            assert( feb == ufeb && vmm == uvmm && ch == uch );
            int urow, ustrip;
            printf( "val=%d", val );
            unpackVal( val, urow, ustrip );
            assert( row == urow && strip == ustrip );
            printf( "(row=%d, strip=%d)\n", row, strip );
        }
    }
    inf.close();
    LOG_INFO << "sTGC Hardware map loaded from File: " << fn << endm;
}

//read the strip center file to get the position information
//for the local coordinate, zero point will be the pin hole, that may changed depends on survey result
bool StFttDb::loadStripCenterFromFile( std::string fn ){
    std::ifstream inf;
    if (mDebug)
        {
            printf( "Opening file: %s \n", fn.c_str());
        }
    inf.open( fn.c_str() );
    if ( !inf ) {
        LOG_WARN << "Strip center file not found" << endm;
        return kFALSE;
    }
    std::string st1 = "Row1";
    std::string st2 = "Row4";

    //check the input file, input file should be Row1 or Row4
    //Row1 for H&V; Row4 for Diagonal
    size_t idx1 = fn.find(st1);
    size_t idx2 = fn.find(st2);
    if( idx1 == string::npos && idx2 == string::npos )
    {
        cout<< "Wrong Input Strip Center File !!!!!!!!!!!" << endl;
        return kFALSE;
    }

    if ( idx1 == string::npos ) // load the Row4(Diagonal)
    {
        scMapDiag.clear();

        //File Header
        std::string nStrip;
        std::string StripCenter;
        inf >> nStrip >> StripCenter;
        if (mDebug)
        {
            printf( "File Header: %s, %s", nStrip.c_str(), StripCenter.c_str());
        }

        //read strip Center info
        int n_strip; Float_t pos_strip_center;
        while (inf >> nStrip >> StripCenter)
        {
            n_strip = atoi(nStrip.c_str());
            pos_strip_center = atof(StripCenter.c_str());

            scMapDiag[n_strip] = pos_strip_center;
        }

    }

    if ( idx2 == string::npos ) // load the Row1(H&V)
    {
        scMapXY.clear();

        //File Header
        std::string nStrip, StripCenter;
        inf >> nStrip >> StripCenter;
        if (mDebug)
        {
            printf( "File Header: %s, %s", nStrip.c_str(), StripCenter.c_str());
        }

        //read strip Center info
        int n_strip; Float_t pos_strip_center;
        while (inf >> nStrip >> StripCenter)
        {
            n_strip = atoi(nStrip.c_str());
            pos_strip_center = atof(StripCenter.c_str());

            scMapXY[n_strip] = pos_strip_center;
        }

    }

    inf.close();
    return kTRUE;
}

//read the strip edge file, it will be used for reject ghost hits
//for the local coordinate, zero point will be the pin hole, that may changed depends on survey result
bool StFttDb::loadStripEdgeFromFile( std::string fn ){
    std::ifstream inf;
    if (mDebug)
        {
            printf( "Opening file: %s \n", fn.c_str());
        }
    inf.open( fn.c_str() );
    if ( !inf ) {
        LOG_WARN << "sTGC Strip edge file not found" << endm;
        return kFALSE;
    }
    std::string st1 = "Row4_edge";

    //check the input file, input file should include Row4_edge
    //Row4 for Diagonal
    size_t idx1 = fn.find(st1);
    if( idx1 == string::npos)
    {
        LOG_ERROR << "Wrong Input Strip Edge File !!!!!!!!!!!" << endm;
        return kFALSE;
    }

    if ( idx1 != string::npos ) // load the Row4(Diagonal)
    {
        seMapDiagLeft.clear();
        seMapDiagRight.clear();

        //File Header
        std::string nStrip, StripEdge_L, StripEdge_R;
        inf >> nStrip >> StripEdge_L >> StripEdge_R;
        if (mDebug)
        {
            printf( "File Header: %s, %s, %s", nStrip.c_str(), StripEdge_L.c_str(), StripEdge_R.c_str());
        }

        //read strip Center info
        int n_strip; Float_t pos_strip_edge_L, pos_strip_edge_R;
        while (inf >> nStrip >> StripEdge_L >> StripEdge_R)
        {
            n_strip = atoi(nStrip.c_str());
            pos_strip_edge_L = atof(StripEdge_L.c_str());
            pos_strip_edge_R = atof(StripEdge_R.c_str());

            seMapDiagLeft[n_strip] = pos_strip_edge_L;
            seMapDiagRight[n_strip] = pos_strip_edge_R;
        }

    }

    inf.close();
    LOG_INFO << "sTGC Strip Edges loaded from File: " << fn << endm;
    return kTRUE;
}

//load the strip length information from the files. this will be used to set as the sigma along the strip direction 
bool StFttDb::loadStripLengthFromFile( std::string fn ){
    std::ifstream inf;
    if (mDebug)
        {
            printf( "Opening file: %s \n", fn.c_str());
        }
    inf.open( fn.c_str() );
    if ( !inf ) {
        LOG_WARN << "sTGC Stirp length file not found" << endm;
        return kFALSE;
    }
    std::string st1 = "Row1";
    std::string st2 = "Row2";
    std::string st3 = "Row3";
    std::string st4 = "Row4";
    std::string st5 = "Row5";

    //check the input file, input file should be Row1 or Row4
    //Row1 for H&V; Row4 for Diagonal
    size_t idx1 = fn.find(st1);
    size_t idx2 = fn.find(st2);
    size_t idx3 = fn.find(st3);
    size_t idx4 = fn.find(st4);
    size_t idx5 = fn.find(st5);
    if( idx1 == string::npos && idx2 == string::npos  && idx3 == string::npos  && idx4 == string::npos  && idx5 == string::npos )
    {
        cout<< "Wrong Input Strip Length File !!!!!!!!!!!" << endl;
        return kFALSE;
    }

    if ( idx1 != string::npos ) // load the Row1(H/V with most strips)
    {
        slMapRow1.clear();

        //File Header
        std::string nStrip;
        std::string StripLength;
        inf >> nStrip >> StripLength;
        if (mDebug)
        {
            printf( "File Header: %s, %s", nStrip.c_str(), StripLength.c_str());
        }

        //read strip Center info
        int n_strip; Float_t pos_strip_length;
        while (inf >> nStrip >> StripLength)
        {
            n_strip = atoi(nStrip.c_str());
            pos_strip_length = atof(StripLength.c_str());

            if (mDebug)
            {
                LOG_INFO << "Strip " << n_strip << " with Length" << pos_strip_length << endm;
            }
            slMapRow1[n_strip] = pos_strip_length;
        }
    }

    if ( idx2 != string::npos ) // load the Row2(H/V with second largest strip group)
    {
        slMapRow2.clear();

        //File Header
        std::string nStrip;
        std::string StripLength;
        inf >> nStrip >> StripLength;
        if (mDebug)
        {
            printf( "File Header: %s, %s", nStrip.c_str(), StripLength.c_str());
        }

        //read strip Center info
        int n_strip; Float_t pos_strip_length;
        while (inf >> nStrip >> StripLength)
        {
            n_strip = atoi(nStrip.c_str());
            pos_strip_length = atof(StripLength.c_str());

            if (mDebug)
            {
                LOG_INFO << "Strip " << n_strip << " with Length" << pos_strip_length << endm;
            }
            slMapRow2[n_strip] = pos_strip_length;
        }
    }

    if ( idx3 != string::npos ) // load the Row3(H/V with third largest strip group)
    {
        slMapRow3.clear();

        //File Header
        std::string nStrip;
        std::string StripLength;
        inf >> nStrip >> StripLength;
        if (mDebug)
        {
            printf( "File Header: %s, %s", nStrip.c_str(), StripLength.c_str());
        }

        //read strip Center info
        int n_strip; Float_t pos_strip_length;
        while (inf >> nStrip >> StripLength)
        {
            n_strip = atoi(nStrip.c_str());
            pos_strip_length = atof(StripLength.c_str());

            if (mDebug)
            {
                LOG_INFO << "Strip " << n_strip << " with Length" << pos_strip_length << endm;
            }
            slMapRow3[n_strip] = pos_strip_length;
        }
    }

    if ( idx4 != string::npos ) // load the Row4(digonal with largest strips)
    {
        slMapRow4.clear();

        //File Header
        std::string nStrip;
        std::string StripLength;
        inf >> nStrip >> StripLength;
        if (mDebug)
        {
            printf( "File Header: %s, %s", nStrip.c_str(), StripLength.c_str());
        }

        //read strip Center info
        int n_strip; Float_t pos_strip_length;
        while (inf >> nStrip >> StripLength)
        {
            n_strip = atoi(nStrip.c_str());
            pos_strip_length = atof(StripLength.c_str());

            if (mDebug)
            {
                LOG_INFO << "Strip " << n_strip << " with Length" << pos_strip_length << endm;
            }
            slMapRow4[n_strip] = pos_strip_length;
        }
    }

    if ( idx5 != string::npos ) // load the Row5(digonal with second largest strips)
    {
        slMapRow5.clear();

        //File Header
        std::string nStrip;
        std::string StripLength;
        inf >> nStrip >> StripLength;
        if (mDebug)
        {
            printf( "File Header: %s, %s", nStrip.c_str(), StripLength.c_str());
        }

        //read strip Center info
        int n_strip; Float_t pos_strip_length;
        while (inf >> nStrip >> StripLength)
        {
            n_strip = atoi(nStrip.c_str());
            pos_strip_length = atof(StripLength.c_str());

            if (mDebug)
            {
                LOG_INFO << "Strip " << n_strip << " with Length" << pos_strip_length << endm;
            }
            slMapRow5[n_strip] = pos_strip_length;
        }
    }

    inf.close();
    LOG_INFO << "sTGC Strip Edges loaded from File: " << fn << endm;
    return kTRUE;
}

// same for all planes
// we have quadrants like:
//
// D | A
// ------
// C | B
// Row 3 and 4 are always diagonal
// odd (even) FOB are horizontal (vertical) for A and C (B and D)
// even (odd) FOB are vertical (horizontal) for A and C (B and D)
//
// REVERTED to original (2026-07-22, akio+Claude): a full H<->V swap of this
// formula was tested at big statistics (script/checkAllExcessQuad.C,
// fttgeom/real_data_excess_check.html) and found NOT to produce a clean
// signal relocation -- traced the full chain (StFttClusterMaker's
// hStripsPerRob/vStripsPerRob bucketing -> FindClusters, identical
// algorithm on both streams -> StFttClusterPointMaker::MakeLocalPoints,
// which assigns the SAME orientation-agnostic clu->x() [=strip*3.2-1.6,
// StFttClusterMaker.cxx CalculateClusterInfo] to global X for kFttVertical
// or global Y for kFttHorizontal) -- a clean swap of this formula SHOULD
// have caused the FST-blind diagnostic's V and H populations to trade
// places entirely (same physical hits, relabeled). They didn't -- both
// states show real but somewhat weaker peaks in the SAME slots, not a
// trade. This means either the real bug isn't a clean global parity-based
// mislabeling (this formula may need per-FEB correction, not a blanket
// flip), or something else not yet traced explains the mismatch. See
// script/testXYMirror.C (2026-07-22) for a cleaner, isolated test that
// mirrors X<->Y purely in the FST-blind diagnostic, without touching this
// formula/clustering/real point objects at all.
UChar_t StFttDb::getOrientation( int rob, int feb, int vmm, int row ) const {
    if ( mDebug ) {
        printf( "getOrientation( %d, %d, %d, %d )", rob, feb, vmm, row ); puts("");
    }

    if ( rob % 2 == 0 ){ // even rob
        if ( feb % 2 == 0 ) { // even feb
            // row 3 and 4 are always diagonal
            if ( 3 == row || 4 == row )
                return kFttDiagonalH;
            return kFttHorizontal;
        }
        // row 3 and 4 are always diagonal
        if ( 3 == row || 4 == row )
            return kFttDiagonalV;
        // even
        return kFttVertical;
    } else { // odd rob

        if ( feb % 2 == 0 ) { // even feb
            // row 3 and 4 are always diagonal
            if ( 3 == row || 4 == row )
                return kFttDiagonalV;
            return kFttVertical;
        }
        // row 3 and 4 are always diagonal
        if ( 3 == row || 4 == row )
            return kFttDiagonalH;
        // even
        return kFttHorizontal;
    }
    // should never get here!
    if ( mDebug ) {
        LOG_DEBUG << "kFttUnknownOrientation = " << kFttUnknownOrientation << endm;
    }
    return kFttUnknownOrientation;
}

/* get 
 * returns the mapping for a given input
 * 
 * input:
 *      rob: 1 - 16
 *      feb: 1 - 6
 *      vmm: 1 - 4
 *      ch : 0 - 63
 *
 * output:
 *      row: 0 - 4
 *      strip: 0 - 162
 *      orientation: one of {Horizontal, Vertical, Diagonal, Unknown}
 *
 */
bool StFttDb::hardwareMap( int rob, int feb, int vmm, int ch, int &row, int &strip, UChar_t &orientation ) const{
    uint16_t key = packKey( feb, vmm, ch );
    if ( mMap.count( key ) ){
        uint16_t val = mMap.at( key );
        unpackVal( val, row, strip );
        orientation = getOrientation( rob, feb, vmm, row );
        return true;
    }
    return false;
}

bool StFttDb::hardwareMap( StFttRawHit * hit ) const{
    uint16_t key = packKey( hit->feb()+1, hit->vmm()+1, hit->channel() );
    if ( mMap.count( key ) ){
        uint16_t val = mMap.at( key );
        int row=-1, strip=-1;
        unpackVal( val, row, strip );
        
        u_char iPlane = hit->sector() - 1;
        u_char iQuad = hit->rdo() - 1;
        int rob = iQuad + ( iPlane *nQuadPerPlane ) + 1;

        UChar_t orientation = getOrientation( rob, hit->feb()+1, hit->vmm()+1, row );
        hit->setMapping( iPlane, iQuad, row, strip, orientation );

        // set strip info
        Float_t stripCenter = -1;
        Float_t stripLeftEdge = -1;
        Float_t stripRightEdge = -1;
        Float_t stripLength = -1;
        if (orientation == kFttHorizontal || orientation == kFttVertical){
            if ( scMapXY.count( strip ) > 0 )
                stripCenter = scMapXY.at(strip);
            else {
                LOG_ERROR << "Cannot find StripCenter for " << strip << endm;
            }
            if ( slMapRow1.count( strip ) > 0 && row == 0)// for Strip Length infomation
                stripLength = slMapRow1.at(strip);
            else if (slMapRow2.count( strip ) > 0 && row == 1)
            {
                stripLength = slMapRow2.at(strip);
            } else if (slMapRow3.count( strip ) > 0 && row == 2)
            {
                stripLength = slMapRow3.at(strip);
            } else {
                LOG_ERROR << "Cannot find StripLength for row " << row << " Strip " << strip << endm;
            }
        }
        if (orientation == kFttDiagonalH || orientation == kFttDiagonalV) {
            if ( scMapDiag.count(strip) > 0 )
                stripCenter    = scMapDiag.at(strip);
            else {
                LOG_ERROR << "Cannot find StripCenter for Diag " << strip << endm;
            }
            if (seMapDiagLeft.count(strip) > 0)
                stripLeftEdge  = seMapDiagLeft.at(strip)-gapPitch/2.;
            else {
                LOG_ERROR << "Cannot find StripLeftEdge for Diag " << strip << endm;
            }
            if (seMapDiagRight.count(strip) > 0)
                stripRightEdge = seMapDiagRight.at(strip)+gapPitch/2.;
            else {
                LOG_ERROR << "Cannot find StripRightEdge for " << strip << endm;
            }
            if (slMapRow4.count(strip) > 0 && row == 3)// for Strip Length infomation
                stripLength = slMapRow4.at(strip);
            else if (slMapRow5.count(strip) > 0 && row == 4)
            {
                stripLength = slMapRow5.at(strip);
            } else
            {
                LOG_ERROR << "Cannot find StripLength for row " << row << " Strip " << strip << endm;
            }
        }
        hit->setStripEdges( stripCenter, stripLeftEdge, stripRightEdge );
        hit->setStripLength(stripLength);

        return true;
    }
    return false;
}
//used to reversve the map, from the hardware to electronic map
// plane, quad, row and strip can be calculated from the simulation maker
// the key issue if to figure out which feb, row, and strip is for the selected channel
bool StFttDb::reverseHardwareMap( int &rob, int &feb, int &vmm, int &ch, int plane, int quad, int row, int strip, UChar_t &orientation ) const {
    uint16_t val = packVal( row, strip );
    if ( !rMap.count( val ) ) return false;

    rob = quad + ( plane *nQuadPerPlane ) + 1;// input plane and quad should start from 0;

    // (row,strip) alone is ambiguous -- both an H and a V (or DiagonalH/DiagonalV
    // for row 3/4) channel exist there in the real hardware map. If the caller
    // pre-set orientation to a specific value, return the channel matching it;
    // if left at kFttUnknownOrientation, return whichever comes first (matches
    // the old, pre-multi-channel behavior for callers that don't care).
    UChar_t wantOrientation = orientation;
    for ( uint16_t key : rMap.at( val ) ) {
        int f, v, c;
        unpackKey( key, f, v, c );
        UChar_t o = getOrientation( rob, f, v, row );
        if ( wantOrientation == kFttUnknownOrientation || o == wantOrientation ) {
            feb = f; vmm = v; ch = c;
            orientation = o;
            return true;
        }
    }
    return false;
}
//used to reversve the map, from the hardware to electronic map
// plane, quad, row and strip can be calculated from the simulation maker
// the key issue if to figure out which feb, row, and strip is for the selected channel
bool StFttDb::reverseHardwareMap( int &feb, int &vmm, int &ch, int row, int strip ) const{
    // No rob/orientation available here to disambiguate which of the (row,strip)
    // channels (H vs V, or DiagonalH vs DiagonalV) is wanted -- returns the first
    // one on record. Prefer the 9-arg overload when orientation matters.
    uint16_t val = packVal( row, strip );
    if ( rMap.count( val ) && !rMap.at( val ).empty() ){
        uint16_t key = rMap.at( val ).front();
        unpackKey( key, feb, vmm, ch );//get the feb, vmm and channel information
        return true;
    }
    return false;
}

UChar_t StFttDb::plane( StFttRawHit * hit ){
    if ( hit->plane() < nPlane )
        return hit->plane();
    return hit->sector() - 1;
}

UChar_t StFttDb::quadrant( StFttRawHit * hit ){
    // Bug (found 2026-07-22): was checking against nQuad (=16, total
    // quadrants across all 4 planes) instead of nQuadPerPlane (=4). Since
    // the "unset" sentinel kFttUnknownQuadrant=4 is < 16, the check always
    // passed and the rdo()-1 fallback below never fired for a hit whose
    // StFttRawHit::mQuadrant hadn't been mapped yet (e.g. anything read
    // before StFttClusterMaker::ApplyHardwareMap() runs this event --
    // notably StFttHitCalibMaker::Make(), which calls StFttDb::fob(),
    // which calls this, and runs earlier in the chain than
    // ApplyHardwareMap in production). Matches StFttDb::plane()'s
    // (correct) use of nPlane just above.
    if ( hit->quadrant() < nQuadPerPlane )
        return hit->quadrant();
    return hit->rdo() - 1;
}

UChar_t StFttDb::rob( StFttRawHit * hit ){
    // NOTE: 1-based, range [1,16] -- NOT the same convention as rob(StFttCluster*) below.
    return quadrant(hit) + ( plane(hit) * nQuadPerPlane ) + 1;
}

UChar_t StFttDb::rob( StFttCluster * clu ){
    // NOTE: 0-based, range [0,15] -- NOT the same convention as rob(StFttRawHit*) above.
    return clu->quadrant() + ( clu->plane() * StFttDb::nQuadPerPlane );
}

UChar_t StFttDb::fob( StFttRawHit * hit ){
    return hit->feb() + ( quadrant( hit ) * nFobPerQuad ) + ( plane(hit) * nFobPerPlane ) + 1;
}

UChar_t StFttDb::orientation( StFttRawHit * hit ){
    if ( hit->orientation() < kFttUnknownOrientation ){
        return hit->orientation();
    }
    return kFttUnknownOrientation;
}

// ---------------------------------------------------------------------------
// Geometry offsets. See the block comment in StFttDb.h for the convention.
// ---------------------------------------------------------------------------
void StFttDb::resetGeometryToHardcoded() {
    for ( size_t ip = 0; ip < nPlane; ip++ ) {
        mXShift[0][ip] = X_shift_QuadA[ip];  mYShift[0][ip] = Y_shift_QuadA[ip];
        mXShift[1][ip] = X_shift_QuadB[ip];  mYShift[1][ip] = Y_shift_QuadB[ip];
        mXShift[2][ip] = X_shift_QuadC[ip];  mYShift[2][ip] = Y_shift_QuadC[ip];
        mXShift[3][ip] = X_shift_QuadD[ip];  mYShift[3][ip] = Y_shift_QuadD[ip];
        mZLoc[0][ip]   = idealPlaneZLocations_QuadA[ip];
        mZLoc[1][ip]   = idealPlaneZLocations_QuadB[ip];
        mZLoc[2][ip]   = idealPlaneZLocations_QuadC[ip];
        mZLoc[3][ip]   = idealPlaneZLocations_QuadD[ip];
    }
    mGeoFromDb = false;
}

// Row i of a Survey table, or 0 if it is missing or too short. Rows are taken in
// table order -- the same thing AGML's <Misalign row="N"/> indexes -- and the
// 1-based Id column is only checked, never used to address.
static Survey_st* fttSurveyRow( St_Survey *t, int irow, const char *name ) {
    if ( !t ) return 0;
    if ( irow < 0 || irow >= t->GetNRows() ) {
        LOG_WARN << "StFttDb: " << name << " has " << ( t ? t->GetNRows() : 0 )
                 << " rows, need row " << irow << " -- treating as identity" << endm;
        return 0;
    }
    Survey_st *r = ((Survey_st*) t->GetTable()) + irow;
    if ( r->Id != irow + 1 )
        LOG_WARN << "StFttDb: " << name << " row " << irow << " has Id " << r->Id
                 << ", expected " << irow + 1 << " (using row order regardless)" << endm;
    return r;
}

// The AGML nominal ("placeholder") position of each pentagon, from StgmGeo1.xml.
// STGM sits in CAVE at (0, 5.9, 338.8385) and the four STFM are placed inside it at
// x = 0, 0, -6.5, +6.5 (pentagon order A, D, C, B) and z = zplane, so the nominal
// global position of a pentagon is (phx[pent], 5.9, 338.8385 + zplane[station]).
//
// These are DESIGN values, not survey: they are the symmetric placeholder AGML starts
// from, and everything real is carried by the misalign tables on top.
//
// They duplicate StgmGeo1.xml and must change with it. That duplication is avoidable:
// the built geometry already has the tables applied, so the STFM_n node position IS
// placeholder + tables -- reading it would need neither these constants nor the table
// walk below. Verified equal on 2026-09-28 (quad A station 1: node (2.152, 12.203,
// 312.342) vs this path (2.15244, 12.2029, 312.342)). Not done that way yet because
// StFttDb would then need gGeoManager loaded before first use, which InitRun cannot
// guarantee; this path stays as the fallback for chains that never load fGeom.
static const double kAgmlPentX[4]  = {  0.0,   0.0,  -6.5,   6.5 };   // cm, by pent (A,D,C,B)
static const double kAgmlPentY     =    5.9;                          // cm, from STGM in CAVE
static const double kAgmlStationZ[4] = { 312.342, 329.953, 347.637, 365.422 };  // cm

int StFttDb::loadGeometryFromDb( St_Survey *stgcOnTpc, St_Survey *stationOnStgc,
                                 St_Survey *pentOnStation ) {
    if ( !mUseDbGeometry ) {
        resetGeometryToHardcoded();
        LOG_INFO << "StFttDb: DB geometry DISABLED, using the hardcoded quadrant offsets" << endm;
        return kStOK;
    }
    // pentOnStation is what carries the per-quadrant position. Without it the AGML
    // placeholder alone is symmetric and would be centimetres wrong, so fall back
    // rather than produce something plausible-looking and wrong.
    if ( !pentOnStation ) {
        resetGeometryToHardcoded();
        LOG_WARN << "StFttDb: no Geometry/stgc/pentOnStation table, "
                 << "falling back to the hardcoded quadrant offsets" << endm;
        return kStWarn;
    }

    // In DB mode the hardcoded StFttDb origins are NOT used. The tables were built as
    //     pentOnStation = (StFttDb origin - AGML placeholder) + alignment
    // so placeholder + table reproduces the true position exactly, and it is the same
    // quantity AGML itself computes -- the geometry and the hit positions then come
    // from one source instead of two that have to be kept in step by hand.
    Survey_st *glob = fttSurveyRow( stgcOnTpc, 0, "stgcOnTpc" );
    const int quad2pent[4] = { 0, 3, 2, 1 };   // StFttDb A,B,C,D -> AGML pent order A,D,C,B

    for ( size_t ip = 0; ip < nPlane; ip++ ) {
        Survey_st *st = fttSurveyRow( stationOnStgc, (int)ip, "stationOnStgc" );
        for ( size_t iq = 0; iq < nQuadPerPlane; iq++ ) {
            const int ipent = quad2pent[iq];
            Survey_st *pe = fttSurveyRow( pentOnStation, 4 * (int)ip + ipent, "pentOnStation" );
            double tx = kAgmlPentX[ipent], ty = kAgmlPentY, tz = kAgmlStationZ[ip];   // cm
            if ( glob ) { tx += glob->t0; ty += glob->t1; tz += glob->t2; }
            if ( st   ) { tx += st->t0;   ty += st->t1;   tz += st->t2;   }
            if ( pe   ) { tx += pe->t0;   ty += pe->t1;   tz += pe->t2;   }
            mXShift[iq][ip] = 10.0 * tx;          // shifts are mm, survey is cm
            mYShift[iq][ip] = 10.0 * ty;
            mZLoc  [iq][ip] =        tz;          // z locations are already cm
        }
    }
    mGeoFromDb = true;

    const char *qn[4] = { "A", "B", "C", "D" };
    LOG_INFO << "StFttDb: quadrant offsets built from AGML placeholder + DB survey tables"
             << " (hardcoded origins NOT used)" << endm;
    for ( size_t jq = 0; jq < nQuadPerPlane; jq++ )
        LOG_INFO << "  quad " << qn[jq]
                 << "  dx(mm) " << mXShift[jq][0] << " " << mXShift[jq][1] << " "
                                << mXShift[jq][2] << " " << mXShift[jq][3]
                 << " | dy(mm) " << mYShift[jq][0] << " " << mYShift[jq][1] << " "
                                 << mYShift[jq][2] << " " << mYShift[jq][3]
                 << " | z(cm) "  << mZLoc[jq][0]   << " " << mZLoc[jq][1]   << " "
                                 << mZLoc[jq][2]   << " " << mZLoc[jq][3] << endm;
    return kStOK;
}

void StFttDb::getGloablOffset( UChar_t plane, UChar_t quad, 
                                float &dx, float &sx,
                                float &dy, float &sy, 
                                float &dz, float &sz ){
    // TODO: connect to DB for calibrated positions. 
    // for now we use the ideal positions (from simulated geometry)
    // calibration will come later

    // scale factors
    sx = 1.0;
    sy = 1.0;
    sz = 1.0;

    // shifts
    dx = 0.0;
    dy = 6.0;
    dz = 0.0;

    if ( plane < 4 )
    {
        // upper quadrants are not displace
        // there have a issue, for the xy shift, the unit is mm, but for the z, the unit is cm 
        // for Z location, suppose that z at the center of the chamber
        if ( quad < nQuadPerPlane ) {
            dx = mXShift[quad][plane];
            dy = mYShift[quad][plane];
            dz = mZLoc[quad][plane] - (LocalStripZLocations[2]+LocalStripZLocations[3])/2.;
        }
    }
        else dz = -999;


    // these are the reflections of a pentagon into the symmetric shape for quadrants A, B, C, D
    if ( quad == 1 )
        sy = -1.0;
    else if ( quad == 2 ){
        sx = -1.0;
        sy = -1.0;
    } else if ( quad == 3 )
        sx = -1.0;

}

void StFttDb::getGloablOffset_ClusterPoint( UChar_t plane, UChar_t quad, 
                                float &dx, float &sx,
                                float &dy, float &sy, 
                                float &dz, float &sz ){
    // TODO: connect to DB for calibrated positions. 
    // for now we use the ideal positions (from simulated geometry)
    // calibration will come later

    // scale factors
    sx = 1.0;
    sy = 1.0;
    sz = 1.0;

    // shifts
    dx = 0.0;
    dy = 6.0;
    dz = 0.0;

    if ( plane < 4 )
        dz = StFttDb::idealPlaneZLocations[plane];

    // upper quadrants are not displace
    // there have a issue, for the xy shift, the unit is mm, but for the z, the unit is cm 
    if ( quad < nQuadPerPlane && plane < nPlane ) {
        dx = mXShift[quad][plane];
        dy = mYShift[quad][plane];
        dz = mZLoc[quad][plane];
    }

    // these are the reflections of a pentagon into the symmetric shape for quadrants A, B, C, D
    if ( quad == 1 )
        sy = -1.0;
    else if ( quad == 2 ){
        sx = -1.0;
        sy = -1.0;
    } else if ( quad == 3 )
        sx = -1.0;

}