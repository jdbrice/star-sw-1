#ifndef FITTERUTILS_H
#define FITTERUTILS_H

#include "GenFit/KalmanFitter.h"
#include "GenFit/KalmanFitterInfo.h"
#include "GenFit/KalmanFitterRefTrack.h"
#include "GenFit/MaterialEffects.h"
#include "GenFit/PlanarMeasurement.h"
#include "GenFit/RKTrackRep.h"
#include "GenFit/SpacepointMeasurement.h"
#include "GenFit/StateOnPlane.h"
#include "GenFit/TGeoMaterialInterface.h"
#include "GenFit/Track.h"
#include "GenFit/TrackPoint.h"
#include "TVector3.h"


class FitSeedMaker {
    public:
        FitSeedMaker() {}
        virtual ~FitSeedMaker() {}    
        virtual void makeSeed(Seed_t seed, TVector3 &posSeed, TVector3 &momSeed, int &q ) = 0;
};
class ConstFitSeeder : public FitSeedMaker {
    public:
        ConstFitSeeder() {}
        virtual ~ConstFitSeeder() {}
        virtual void makeSeed(Seed_t seed, TVector3 &posSeed, TVector3 &momSeed, int &q ) {
            posSeed.SetXYZ(0,0,0);
            momSeed.SetXYZ(0,0,10);
            q = 1;
        }
};
class GenericFitSeeder : public FitSeedMaker {
    bool isValid = false;
    public:
        GenericFitSeeder() {}
        virtual ~GenericFitSeeder() {}
        // Sentinel returned by computeSignedCurvature() for collinear (undefined-curvature)
        // point triples. Must be excluded in averageCurvature() -- previously that filter
        // checked "!= -1" instead of this sentinel, so collinear triples (guaranteed for
        // genuinely straight tracks, e.g. B=0) leaked through as a bogus near-zero curvature,
        // blowing up the seed pT estimate (pt = K*B/curvature).
        static constexpr double kCollinearCurvature = -1e-7;
        // Simple function to calculate the determinant of a 2x2 matrix
        inline double determinant(double a, double b, double c, double d) {
            return a * d - b * c;
        }
        struct Point {
            double x, y;
        };
        // Function to compute the curvature of a circle given 3 points
        double computeSignedCurvature(const Point& p1, const Point& p2, const Point& p3) {
            // Calculate the lengths of the sides of the triangle
            double A = std::sqrt(std::pow(p2.x - p1.x, 2) + std::pow(p2.y - p1.y, 2));
            double B = std::sqrt(std::pow(p3.x - p2.x, 2) + std::pow(p3.y - p2.y, 2));
            double C = std::sqrt(std::pow(p1.x - p3.x, 2) + std::pow(p1.y - p3.y, 2));

            // Calculate the determinant of the matrix formed by the points
            double det = determinant(p2.x - p1.x, p2.y - p1.y, p3.x - p1.x, p3.y - p1.y);
            // LOG_INFO << "Det: " << det << endm;
            double charge = det > 0 ? -1 : 1;
            // Area of the triangle formed by the three points
            double area = std::abs(det) / 2.0;

            if (area == 0) {
                // Collinear seed triple (routine for low-pT/straight tracks) -- handled via
                // the sentinel return below, not logged (used to print every occurrence,
                // which alone produced ~660k lines/run in dense MC samples).
                return kCollinearCurvature; // Curvature is undefined for collinear points
            }

            // Calculate the radius of the circumcircle using the formula:
            // R = (A * B * C) / (4 * area)
            double radius = (A * B * C) / (4 * area);
            // LOG_INFO << "Radius: " << radius << endm;
            // Curvature is the inverse of the radius
            return charge / radius;
        }
        // Function to compute the average curvature for all combinations of 3 points
        double averageCurvature( const Seed_t points ) {
            // const std::vector<Point>& points;
            int numPoints = points.size();
            if (numPoints < 3) {
                std::cerr << "Not enough points to form a circle." << std::endl;
                return -1;
            }

            double totalCurvature = 0.0;
            int validCombinations = 0;

            // Iterate over all combinations of 3 points
            for (int i = 0; i < numPoints - 2; ++i) {
                for (int j = i + 1; j < numPoints - 1; ++j) {
                    for (int k = j + 1; k < numPoints; ++k) {
                        Point p0 = {points[i]->getX(), points[i]->getY()};
                        Point p1 = {points[j]->getX(), points[j]->getY()};
                        Point p2 = {points[k]->getX(), points[k]->getY()};
                        double curvature = computeSignedCurvature(p0, p1, p2);
                        // printf("Curvature for points (%d, %d, %d): %f\n", i, j, k, curvature);
                        if ( i != 0 || k != numPoints - 1 ) { // skip if not using the first and last point
                            // This improves the seed charge estimate substantially for the beamline / primary tracks
                            // momentum resolution also benefits somewhat
                            // printf("Skipping non extreme points \n");
                            continue;
                        }
                        if (curvature != kCollinearCurvature) {  // Exclude invalid (collinear) combinations
                            totalCurvature += curvature;
                            ++validCombinations;
                        }
                    }
                }
            }

            if (validCombinations == 0) {
                std::cerr << "No valid curvature calculations possible." << std::endl;
                return -1;  // No valid triangles were found
            }

            return totalCurvature / validCombinations;
        }

        template <typename T> int sgn(T val) {
            return (T(0) < val) - (val < T(0));
        }
        virtual void makeSeed(Seed_t seed, TVector3 &posSeed, TVector3 &momSeed, int &q ) {
            const double qc = averageCurvature(seed);
            LOG_DEBUG << "GenericFitSeeder::makeSeed::Curvature: " << qc << endm;
            // posSeed.SetXYZ(seed[0]->getX(), seed[0]->getY(), seed[0]->getZ());
            momSeed.SetXYZ(0,0,10);
        
            // const double BStrength = 0.5; // 0.5 T = 5 Gauss
            // const double C = 0.3 * BStrength; //C depends on the units used for momentum and Bfield (here GeV and Tesla)
            //AAA fix: comment was wrong — K is in units of GeV·cm/kGauss, B=5 is kGauss (5 kG = 0.5 T), not Gauss.
            //         Numerically correct: pt = K*B[kG]/kappa[cm^-1] gives GeV/c.
            //const double K = 0.00029979; // K depends on the units used for Bfield and momentum (here Gauss and GeV)
            const double K = 0.00029979; // momentum in GeV/c, Bfield in kGauss (B=5 kG = 0.5 T)
            // qc == -1 means averageCurvature() found no usable (non-collinear) point triple --
            // e.g. too few points, or a genuinely straight track (guaranteed collinear at B=0).
            // There is no curvature to convert to pT in that case; fall back to a high-pT
            // (~straight-track) default instead of dividing by the sentinel, which previously
            // produced a degenerate ~1.5 MeV seed and crashed downstream field-map lookups
            // with a NaN assertion (StarMagField::Search).
            bool curvatureKnown = (qc != -1.0);
            double pt = curvatureKnown ? fabs((K*5)/qc) : 10.0; // pT from average measured curv
            LOG_DEBUG << "GenericFitSeeder::makeSeed::pt = " << pt << endm;
            // set the momentum seed's transverse momentum
            //AAA fix: this line is dead code — overwritten two lines below by the proper SetXYZ call.
            //momSeed.SetXYZ(pt/sqrt(2.0),pt/sqrt(2.0),10);
            // Use disk 0 and disk 2 (outermost pair) for theta/phi estimation.
            // Using the wider baseline gives a better eta estimate and avoids
            // the degenerate case where disk 0 and disk 1 rasterize to the
            // same strip center (Rxy=0 → tan(theta)=0 → pz blows up).
            TVector3 p0 = TVector3(seed[0]->getX(), seed[0]->getY(), seed[0]->getZ());
            TVector3 p2 = TVector3(seed[2]->getX(), seed[2]->getY(), seed[2]->getZ());
            double dx = (p2.X() - p0.X());
            double dy = (p2.Y() - p0.Y());
            double dz = (p2.Z() - p0.Z());
            double phi = TMath::ATan2(dy, dx);
            double Rxy = sqrt(dx * dx + dy * dy);
            double theta = TMath::ATan2(Rxy, dz);
            if (abs(dx) < 1e-6 || abs(dy) < 1e-6){
                phi = TMath::ATan2( p2.Y(), p2.X() );
            }

            // momSeed.SetPhi(phi);
            // momSeed.SetTheta(theta);
            // Guard the case noted above (Rxy~0 -> theta~0/pi -> tan(theta)~0): dividing pt by
            // it blows up to +-inf/NaN, which previously crashed downstream field-map lookups
            // (StarMagField::Search NaN assertion). Fall back to a fixed high-|pz|, sign matched
            // to the z-direction of travel, instead of dividing by a near-zero tangent.
            double tanTheta = tan(theta);
            double pz = (fabs(tanTheta) > 1e-6) ? (pt / tanTheta) : ((dz >= 0) ? 10.0 : -10.0);
            momSeed.SetXYZ(pt * cos(phi), pt * sin(phi), pz);
            
            // assign charge based on sign of curvature; default to +1 when curvature is unknown
            // (sgn(-1) would otherwise always return -1, silently biasing every degenerate seed negative)
            q = curvatureKnown ? sgn<double>(qc) : 1;
        }
};

#endif