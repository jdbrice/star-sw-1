#ifndef FWD_GEOM_UTILS_H
#define FWD_GEOM_UTILS_H

#include "TGeoVolume.h"
#include "TGeoNode.h"
#include "TGeoMatrix.h"
#include "TGeoNavigator.h"
#include "TGeoTube.h"   // TGeoTubeSeg: the FST sensor shapes, used for their phi limits
#include "StEvent/StFstConsts.h"   // kFstNumSensors
#include "StMessMgr.h"             // LOG_INFO / LOG_WARN / LOG_ERROR

class FwdGeomUtils {
    public:



        FwdGeomUtils( TGeoManager * gMan ) {
            if ( gMan != nullptr ){
                _navigator = gMan->AddNavigator();
                _gMan = gMan;
            }
        }

        ~FwdGeomUtils(){
            if ( _gMan != nullptr && _navigator != nullptr){
                _gMan->RemoveNavigator( _navigator );
            }
        }

        bool cd( const char* path ){
            // Change to the specified path
            bool ret = _navigator -> cd(path);
            // If successful, set the node, the volume, and the GLOBAL transformation
            // for the requested node.  Otherwise, invalidate these
            if ( ret ) {
                _matrix = _navigator->GetCurrentMatrix();
                _node   = _navigator->GetCurrentNode();
                _volume = _node->GetVolume();
            } else {
                _matrix = 0; _node = 0; _volume = 0;
            }
            return ret;
        }

        vector<double> fttZ( vector<double> defaultZ ) {
            double z0 = fttZ(0);
            if ( z0 > 1.0 ) { // returns 0 on faiure
                vector<double> z = {z0, fttZ(1), fttZ(2), fttZ(3)};
                return z;
            }
            return defaultZ;
        }
        double fttZ( int index ) {

            // This ftt_z_delta is needed to match the z location of hits (midpint of active volume?) to the z location of the mother volume.
            // NOTE: It may be possible to improve this when the higher precision FTT geometry model is added
            const double ftt_z_delta = -0.5825245;
            stringstream spath;
            spath << "/HALL_1/CAVE_1/STGM_1/STFM_" << (index + 1) * 4 << "/";
            bool can = cd( spath.str().c_str() );
            if ( can && _matrix != nullptr ){
                return _matrix->GetTranslation()[2] + ftt_z_delta;
            }
            return 0.0;
        }


        vector<double> fstZ( vector<double> defaultZ ) {
            double z0 = fstZ(0);
            if ( z0 > 1.0 ) { // returns 0 on faiure
                vector<double> z = {z0, fstZ(1), fstZ(2)};
                return z;
            }
            return defaultZ;
        }

        double fstZ( int index ) {
            // starting in FtsmGeom v1.16? or 1.17
            // the index are now 4,5,6
            // hence +4 below
            // also fixed typo, previously was incorrectly FTSD_
            const double z_delta = 1.755;
            stringstream spath;
            spath << "/HALL_1/CAVE_1/FSTM_1/FSTD_" << (index + 4) << "/";
            bool can = cd( spath.str().c_str() );
            if ( can && _matrix != nullptr ){
                return _matrix->GetTranslation()[2] + z_delta;
            }
            return 0.0;
        }

        TVector3 getFttQuadrant( int index, TVector3 &u, TVector3 &v){
            // 0 - 15 is the front face
            // 16 - 31 is the back face

            int iquad = index % 16 + 1; // geometry is 1 - 16
            int iplane = index / 16 + 1; // geometry is 1 - 2

            stringstream spath;
            spath << "/HALL_1/CAVE_1/STGM_1/STFM_" << (iquad) << "/STMG_" << iplane << "/";
            bool can = cd( spath.str().c_str() );
            if ( can && _matrix != nullptr ){
                double x = _matrix->GetTranslation()[0];
                double y = _matrix->GetTranslation()[1];
                double z = _matrix->GetTranslation()[2];
                // Column 0 of R = local x-axis in global space = u
                // Column 1 of R = local y-axis in global space = v
                // GetRotationMatrix() is row-major: element [i*3+j] = R[i][j]
                // Column j is elements [0*3+j], [1*3+j], [2*3+j]
                u.SetXYZ(_matrix->GetRotationMatrix()[0], _matrix->GetRotationMatrix()[3], _matrix->GetRotationMatrix()[6]);
                v.SetXYZ(_matrix->GetRotationMatrix()[1], _matrix->GetRotationMatrix()[4], _matrix->GetRotationMatrix()[7]);
                return TVector3(x, y, z);
            }
            std ::cerr << "Failed to get FTT quadrant origin for index " << index << std::endl;
            return TVector3(0,0,0);
        }

        // ---------------------------------------------------------------------
        // FST sensor lookup, driven by the geometry rather than by a hardcoded
        // volume path.
        //
        // The path differs between FST geometry variants, so a fixed string
        // cannot work for all of them:
        //    FSTMv1  (full, ideal)  HALL/CAVE_1/FSTM_1/FSTD_4/FSTW_1/FTUS_1
        //    FSTMv1sm / v2sm (misaligned)   HALL/CAVE_1/FSTM_1/FTUS_1 .. FTUS_108
        // The misaligned variants drop the FSTD disk mother and lift the sensors
        // out of FSTW up into FSTM, so that each of the 108 can carry its own
        // <Misalign> matrix. The sensor ordering inside a wedge is opposite
        // between the two as well (full: FTUS_1,2 outer and FTUS_3 inner;
        // misaligned: copy 1+3w+36d inner, then the two outer). So neither the
        // path nor the copy number is a safe key.
        //
        // What IS stable is the physics: a sensor is identified by where it is
        // and what shape it has.
        //    disk    from its global z            (~150 / ~165 / ~179)
        //    wedge   from its local +x axis in global, which is the wedge
        //            bisector and sits at 75 - 30*w degrees (unaffected by the
        //            180 deg x-rotation that flips front/back wedges)
        //    sensor  from the shape: rmin 5 -> inner (hit sensor 0), rmin 16.5
        //            with phi1 > 180 -> hit sensor 1, else hit sensor 2
        //            (that assignment was measured from data, script/fstLocalMap.C)
        // This also retires the kElecToGeantWedge translation tables, which only
        // existed because copy numbers under the old hierarchy did not follow phi.
        // ---------------------------------------------------------------------
        void walkFstSensors( TGeoNode* node, TGeoHMatrix mat, int depth ){
            TGeoHMatrix here = mat; here.Multiply( node->GetMatrix() );
            TString vname = node->GetVolume()->GetName();
            if ( vname.BeginsWith("FTUS") ){
                TGeoShape* shape = node->GetVolume()->GetShape();
                if ( !shape || !shape->InheritsFrom("TGeoTubeSeg") ) return;
                TGeoTubeSeg* seg = (TGeoTubeSeg*)shape;
                const double* t = here.GetTranslation();
                const double* r = here.GetRotationMatrix();
                int disk   = ( t[2] < 158.0 ) ? 0 : ( ( t[2] < 172.0 ) ? 1 : 2 );
                double phiW = TMath::ATan2( r[3], r[0] ) * TMath::RadToDeg();
                if ( phiW < 0 ) phiW += 360.0;
                int wedge = (int)TMath::Nint( (75.0 - phiW) / 30.0 );
                wedge = ( (wedge % 12) + 12 ) % 12;
                int sensor = ( seg->GetRmin() < 10.0 ) ? 0
                           : ( ( seg->GetPhi1() > 180.0 ) ? 1 : 2 );
                int idx = disk * 36 + wedge * 3 + sensor;
                if ( idx >= 0 && idx < kFstNumSensors ){
                    if ( _fstSensorOK[idx] ) {
                        LOG_WARN << "FwdGeomUtils: duplicate FST sensor index " << idx
                                 << " (disk " << disk << " wedge " << wedge
                                 << " sensor " << sensor << ")" << endm;
                    }
                    _fstSensorMat[idx]  = here;
                    _fstSensorPhi1[idx] = seg->GetPhi1();
                    _fstSensorPhi2[idx] = seg->GetPhi2();
                    _fstSensorOK[idx]   = true;
                }
                return;
            }
            if ( depth > 8 ) return;
            for (int i = 0; i < node->GetNdaughters(); i++)
                walkFstSensors( node->GetDaughter(i), here, depth + 1 );
        }

        void buildFstSensorMap(){
            if ( _fstSensorMapped ) return;
            _fstSensorMapped = true;
            for (int i = 0; i < kFstNumSensors; i++) _fstSensorOK[i] = false;
            if ( !gGeoManager || !gGeoManager->GetTopNode() ){
                LOG_ERROR << "FwdGeomUtils: no geometry, cannot map FST sensors" << endm;
                return;
            }
            TGeoHMatrix identity;
            walkFstSensors( gGeoManager->GetTopNode(), identity, 0 );
            int n = 0;
            for (int i = 0; i < kFstNumSensors; i++) if ( _fstSensorOK[i] ) n++;
            LOG_INFO << "FwdGeomUtils: mapped " << n << " / " << kFstNumSensors
                     << " FST sensors from the geometry" << endm;
            if ( n != kFstNumSensors )
                LOG_ERROR << "FwdGeomUtils: incomplete FST sensor map (" << n
                          << "/" << kFstNumSensors << ") -- geometry not understood" << endm;
        }

        // phi1/phi2 (optional, degrees) come back as the SENSOR SHAPE's azimuthal
        // limits in its own local frame -- that is where AGML keeps the 1 deg
        // outer-sensor gap (outer shapes are +-0.5..15.5, inner is a clean +-15).
        TVector3 getFstSensorOrigin (int index, TVector3 &u, TVector3 &v,
                                     double *phi1 = 0, double *phi2 = 0) {
            buildFstSensorMap();
            if ( index < 0 || index >= kFstNumSensors || !_fstSensorOK[index] ){
                std::cerr << "Failed to get FST sensor origin for index " << index << std::endl;
                return TVector3(0,0,0);
            }
            const double* t = _fstSensorMat[index].GetTranslation();
            const double* r = _fstSensorMat[index].GetRotationMatrix();
            // Column 0 of R = local x-axis in global space = u
            // Column 1 of R = local y-axis in global space = v
            // V is deliberately NOT normalised to counterclockwise: the front/back
            // wedge flip is exactly what this rotation encodes, and normalising it
            // away would force the strip decode to re-supply the orientation from
            // hardcoded constants.
            u.SetXYZ( r[0], r[3], r[6] );
            v.SetXYZ( r[1], r[4], r[7] );
            if (phi1) *phi1 = _fstSensorPhi1[index];
            if (phi2) *phi2 = _fstSensorPhi2[index];
            if ( _verbose ){
                LOG_INFO << "FST Sensor " << index << " origin: " << t[0] << ", " << t[1]
                         << ", " << t[2] << " shapePhi [" << _fstSensorPhi1[index] << ","
                         << _fstSensorPhi2[index] << "]" << endm;
            }
            return TVector3( t[0], t[1], t[2] );
        }

    protected:
    // geometry-driven FST sensor map, indexed by disk*36 + wedge*3 + sensor
    TGeoHMatrix _fstSensorMat[kFstNumSensors];
    double      _fstSensorPhi1[kFstNumSensors] = {0};
    double      _fstSensorPhi2[kFstNumSensors] = {0};
    bool        _fstSensorOK[kFstNumSensors]   = {false};
    bool        _fstSensorMapped = false;

    TGeoVolume    *_volume    = nullptr;
    TGeoNode      *_node      = nullptr;
    TGeoHMatrix   *_matrix    = nullptr;
    TGeoIterator  *_iter      = nullptr;
    TGeoNavigator *_navigator = nullptr;
    TGeoManager   *_gMan      = nullptr;

    const int _verbose = 1;
};

#endif
