#ifndef genfit_STARField_h
#define genfit_STARField_h

#include "TVector3.h"
#include "StarMagField/StarMagField.h"
#include "GenFit/AbsBField.h"

//_______________________________________________________________________________________
// Adaptor for STAR magnetic field loaded via StarMagField Maker
class StarFieldAdaptor : public genfit::AbsBField {
  public:
    StarFieldAdaptor() {};

    virtual TVector3 get(const TVector3 &position) const {
        double x[] = {position[0], position[1], position[2]};
        double B[] = {0, 0, 0};

        get( x[0], x[1], x[2], B[0], B[1], B[2] );

        return TVector3(B);
    };

    inline virtual void get(const double &_x, const double &_y, const double &_z, double &Bx, double &By, double &Bz) const {
        double x[] = {_x, _y, _z};
        double B[] = {0, 0, 0};

        if (StarMagField::Instance()){

            float z = x[2];
            // This used to short-circuit the field map inside the TPC volume
            // (|z| < 250 && r < 50) and return a hardcoded Bz = 4.97979927 kG,
            // i.e. the nominal POSITIVE full field as a constant. That ignores
            // the StarMagField scale factor, and so both the polarity and the
            // magnitude of the run being reconstructed. All three FST disks
            // (z = 151-179, r = 5-28) sit inside that box, so for FST hits the
            // map was never consulted at all.
            //
            // It was harmless until now, which is why nothing looked wrong:
            //   - bfc/sim chains run at scale +1 ("Scale factor = 1.000000,
            //     bfield_full_positive_2D.dat, Bz(0) = 4.9798"), so the constant
            //     happened to equal the truth;
            //   - the MuDst afterburner macros never construct a StarMagField at
            //     all, so Instance() is null, this whole block is skipped and B
            //     stays {0,0,0}. The null instance MASKED the hardcode.
            // Once the afterburner is given a real instance (see
            // macro/mudst/fwd_afterburner_db.C) this becomes wrong. Measured at
            // the FST with the 2022 ReversedFullField production:
            //     map -4.990 kG   vs   constant +4.980   (opposite sign)
            // and for a field-off run: map +0.030 kG vs constant +4.980.
            //
            // The map is correct in this region, so just use it. Cost for MC is a
            // 0.06% shift (+4.9798 constant -> +4.9826 from the map).
            if ( fabs(z) > 450 ) {
                // Outside the mapped volume: keep returning zero rather than
                // letting the extrapolator run on edge-extrapolated values.
                B[0] = 0.; B[1] = 0.; B[2] = 0.;
            } else {
                StarMagField::Instance()->Field(x, B);
            }
        }
        

        Bx = B[0];
        By = B[1];
        Bz = B[2];
        return;
    };
};


#endif
