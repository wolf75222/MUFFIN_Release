//
//  TVDInnerCells.hpp
//  Muffin
//
//  Limited slopes of the inner cells, the loop of the TVD2ndOrder*1D reconstructors.
//

#ifndef TVDInnerCells_hpp
#define TVDInnerCells_hpp

#include "Limiter/Limiter.hpp"
#include "Limiter/ThirdOrderLimiter.hpp"
#include "Limiter/minmodLimiter.hpp"
#include "../Parallel.hpp"

// Left and right values of the cells begin to end - 1 from the cell values P at the centres x: the operations of the
// original loop in the same order. thirdOrder and minmod are inlined, theta is computed for every cell with a
// denominator that is never zero, then kept only where the original computes it, and the choices are masks
// (selectValue): the loop has neither call nor branch, and the compiler vectorises it.
template <class Lim>
inline void tvdInnerCellsInline(const double* __restrict__ P, const double* __restrict__ x, double* __restrict__ R,
                                double* __restrict__ L, int begin, int end)
{
    for (int i = begin; i < end; i++)
    {
        const double a  = (P[i] - P[i - 1])/(x[i] - x[i - 1]);
        const double b  = (P[i + 1] - P[i])/(x[i + 1] - x[i]);
        const double ba = b/selectValue(a != 0.0, a, 1.0);
        const double ab = a/selectValue(b != 0.0, b, 1.0);
        const double limiterR = Lim::value(selectValue(a != 0.0, ba, 0.));
        const double limiterL = Lim::value(selectValue(b != 0.0, ab, 0.));
        R[i] = P[i] + 0.5*limiterR*(P[i] - P[i - 1]);
        L[i] = P[i] - 0.5*limiterL*(P[i + 1] - P[i]);
    }
}

// Cells 1 to n1 - 1 of the nEq equations, equation iEq at iEq*nCells in P, R and L; the cells are shared among
// the threads (Parallel.hpp) when the limiter is inlined.
template <class Lim>
inline void tvdInnerCellsEquations(const double* P, const double* x, double* R, double* L, int n1, int nEq, int nCells)
{
    parallelRange(1, n1, [&](int begin, int end, int) {
        for (int iEq = 0; iEq < nEq; ++iEq)
            tvdInnerCellsInline<Lim>(P + iEq*nCells, x, R + iEq*nCells, L + iEq*nCells, begin, end);
    });
}

inline void tvdInnerCells(Limiter& phi, const double* P, const double* x, double* R, double* L, int n1, int nEq,
                          int nCells)
{
    if (dynamic_cast<ThirdOrderLimiter*>(&phi)) return tvdInnerCellsEquations<ThirdOrderLimiter>(P, x, R, L, n1, nEq, nCells);
    if (dynamic_cast<minmodLimiter*>(&phi)) return tvdInnerCellsEquations<minmodLimiter>(P, x, R, L, n1, nEq, nCells);
    for (int iEq = 0; iEq < nEq; ++iEq)
    {
        const double* Pe = P + iEq*nCells;
        double* Re = R + iEq*nCells;
        double* Le = L + iEq*nCells;
        for (int i = 1; i < n1; i++)
        {
            double a        = (Pe[i] - Pe[i - 1])/(x[i] - x[i - 1]);
            double b        = (Pe[i + 1] - Pe[i])/(x[i + 1] - x[i]);
            double theta    = (a != 0.0) ? b/a : 0;
            double limiterR = phi(theta);
            theta           = (b != 0.0) ? a/b : 0;
            double limiterL = phi(theta);
            Re[i] = Pe[i] + 0.5*limiterR*(Pe[i] - Pe[i - 1]);
            Le[i] = Pe[i] - 0.5*limiterL*(Pe[i + 1] - Pe[i]);
        }
    }
}

#endif /* TVDInnerCells_hpp */
