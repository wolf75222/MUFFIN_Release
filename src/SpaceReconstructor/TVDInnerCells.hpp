//
//  TVDInnerCells.hpp
//  Muffin
//
//  Limited slopes of the inner cells of one equation, the loop of the TVD2ndOrder*1D reconstructors.
//

#ifndef TVDInnerCells_hpp
#define TVDInnerCells_hpp

#include "Limiter/Limiter.hpp"
#include "Limiter/ThirdOrderLimiter.hpp"
#include "Limiter/minmodLimiter.hpp"

// Left and right values of the cells 1 to n1 - 1 from the cell values P at the centres x: the operations of the
// original loop in the same order. thirdOrder and minmod are inlined, and theta is computed for every cell with a
// denominator that is never zero, then kept only where the original computes it: the loop has neither call nor
// branch, and the compiler vectorises it.
template <class Lim>
inline void tvdInnerCellsInline(const double* __restrict__ P, const double* __restrict__ x, double* __restrict__ R,
                                double* __restrict__ L, int n1)
{
    for (int i = 1; i < n1; i++)
    {
        const double a  = (P[i] - P[i - 1])/(x[i] - x[i - 1]);
        const double b  = (P[i + 1] - P[i])/(x[i + 1] - x[i]);
        const double ba = b/((a != 0.0) ? a : 1.0);
        const double ab = a/((b != 0.0) ? b : 1.0);
        const double limiterR = Lim::value((a != 0.0) ? ba : 0);
        const double limiterL = Lim::value((b != 0.0) ? ab : 0);
        R[i] = P[i] + 0.5*limiterR*(P[i] - P[i - 1]);
        L[i] = P[i] - 0.5*limiterL*(P[i + 1] - P[i]);
    }
}

inline void tvdInnerCells(Limiter& phi, const double* P, const double* x, double* R, double* L, int n1)
{
    if (dynamic_cast<ThirdOrderLimiter*>(&phi)) return tvdInnerCellsInline<ThirdOrderLimiter>(P, x, R, L, n1);
    if (dynamic_cast<minmodLimiter*>(&phi)) return tvdInnerCellsInline<minmodLimiter>(P, x, R, L, n1);
    for (int i = 1; i < n1; i++)
    {
        double a        = (P[i] - P[i - 1])/(x[i] - x[i - 1]);
        double b        = (P[i + 1] - P[i])/(x[i + 1] - x[i]);
        double theta    = (a != 0.0) ? b/a : 0;
        double limiterR = phi(theta);
        theta           = (b != 0.0) ? a/b : 0;
        double limiterL = phi(theta);
        R[i] = P[i] + 0.5*limiterR*(P[i] - P[i - 1]);
        L[i] = P[i] - 0.5*limiterL*(P[i + 1] - P[i]);
    }
}

#endif /* TVDInnerCells_hpp */
