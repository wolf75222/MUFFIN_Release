//
//  QMOMMaxEigenvalue.hpp
//  Muffin
//
//  Largest |eigenvalue| of the 5M and 6M HyQMOM flux (FiveMomentQMOM1D, SixMomentAQMOM1D), with the
//  realizability test that makes getMaxEigenvalue print a warning.
//

#ifndef QMOMMaxEigenvalue_hpp
#define QMOMMaxEigenvalue_hpp

#include <algorithm>
#include <cmath>
#include "../MeshData/CellDataRef.hpp"

// realizable is false where the original code printed "out of the realizability" with q_nStar and r_nStar
inline double qmomMaxEigenvalue(const CellDataRef& u, bool& realizable, double& q_nStar_out, double& r_nStar_out)
{
    using std::max; using std::abs; using std::sqrt; using std::pow;
    const double rho    = u[0];
    const double v      = u[1]/rho;
    const double p      = u[2] - rho*v*v;
    const double q      = u[3] - rho*v*v*v - 3*v*p;
    const double r      = u[4] - rho*v*v*v*v - 6*v*v*p - 4*v*q;
    const double r_Max  =  3*p*p/rho;
    const double K      = r - r_Max;
    const double a_n     = sqrt(p/rho);
    const double q_nStar = q/(rho*pow(a_n,3));
    const double K_nStar = K/(rho*pow(a_n,4));
    const double r_nStar = K_nStar + 3.;

    const double xi_nStar = max(r_nStar - pow(q_nStar,2) - 1., 0.); // It should be positive, otherwise it is numerical error
    realizable  = !(r_nStar - pow(q_nStar,2) - 1. < 0.);
    q_nStar_out = q_nStar;
    r_nStar_out = r_nStar;
    const double A = q_nStar/2.;
    const double B = 1. + xi_nStar + pow(q_nStar,2)/4.;
    const double C = sqrt(xi_nStar*(1 + xi_nStar));

    return max(abs(v + a_n*(A + sqrt(B + C))), max(abs(v + a_n*(A + sqrt(B - C))), max(abs(v + a_n*(A - sqrt(B + C))), abs(v + a_n*(A - sqrt(B - C))))));
}

#endif /* QMOMMaxEigenvalue_hpp */
