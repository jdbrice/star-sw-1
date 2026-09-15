#ifndef genfit_STARField_h
#define genfit_STARField_h

#include "TVector3.h"
#include "StarMagField/StarMagField.h"
#include "GenFit/AbsBField.h"

//_______________________________________________________________________________________
// Adaptor for STAR magnetic field loaded via StarMagField Maker
//
// B = {0,0,0} if there is no StarMagField instance (the MuDst afterburner macros
// must create one, see macro/mudst/fwd_afterburner_db.C).
// Otherwise the STAR map at the instance's scale, with two options
// (config keys read in TrackFitter.h):
//
//   fstConstBz ("TrackFitter:fieldFstConstBz", default true)
//       Inside the FST box |z| < 250 && r < 50 use a uniform Bz equal to the map's
//       Bz(0,0,0) at the current scale, i.e. with the run's sign and magnitude:
//       -4.9836 kG for the 2022 reversed full field, +4.9798 for MC at scale +1,
//       ~0 for field-off runs. This replaces a hardcoded +4.97979927 kG, which was
//       the wrong sign for reversed-field runs (including BFC production) and wrong
//       for field-off runs. It halves StFwdTrackMaker CPU relative to the map
//       (run 23081049, 40 events: 207-223 s -> 108-110 s) with no change in track
//       count, pT or chi2. false = use the map everywhere.
//
//   zMaxField  ("TrackFitter:fieldZMax", default 450 cm)
//       B = 0 for |z| > zMaxField; <= 0 means no cut. The map beyond 450 cm is
//       the real extended map (reversed full field: -0.4 kG at z=450, -0.04 kG at
//       z=700), not an edge extrapolation. Removing the cut changed FCS projection
//       failures by < 0.1% and cost no CPU with fstConstBz=true.
class StarFieldAdaptor : public genfit::AbsBField {
  public:
    StarFieldAdaptor( bool fstConstBz = true, double zMaxField = 450. )
      : mFstConstBz(fstConstBz), mZMaxField(zMaxField), mCachedFactor(-999.), mBz0(0.) {};

    virtual TVector3 get(const TVector3 &position) const {
        double x[] = {position[0], position[1], position[2]};
        double B[] = {0, 0, 0};

        get( x[0], x[1], x[2], B[0], B[1], B[2] );

        return TVector3(B);
    };

    inline virtual void get(const double &_x, const double &_y, const double &_z, double &Bx, double &By, double &Bz) const {
        double x[] = {_x, _y, _z};
        double B[] = {0, 0, 0};

        StarMagField *mf = StarMagField::Instance();
        if (mf){
            double az = fabs(x[2]);
            if ( mFstConstBz && az < 250. && x[0]*x[0] + x[1]*x[1] < 2500. ) {
                if ( mf->GetFactor() != mCachedFactor ) {     // scale can change at InitRun
                    double o[3] = {0, 0, 0}, b0[3] = {0, 0, 0};
                    mf->Field(o, b0);
                    mBz0 = b0[2];
                    mCachedFactor = mf->GetFactor();
                }
                B[2] = mBz0;
            } else if ( mZMaxField > 0 && az > mZMaxField ) {
                // B stays {0,0,0}
            } else {
                mf->Field(x, B);
            }
        }

        Bx = B[0];
        By = B[1];
        Bz = B[2];
        return;
    };

  private:
    bool           mFstConstBz;
    double         mZMaxField;
    mutable double mCachedFactor;
    mutable double mBz0;
};


#endif
