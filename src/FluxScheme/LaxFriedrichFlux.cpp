//
//  LaxFriedrichFlux.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 11/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "LaxFriedrichFlux.hpp"
#include "FluxSchemeRegistrar.hpp"
#include <cmath>
#include <iostream>

using namespace std;

REGISTER_FLUXSCHEME("LaxFriedrich", LaxFriedrichFlux);

vector<double>& LaxFriedrichFlux::operator()(const CellDataRef uL, const CellDataRef uR)
{
    
    double maxEigenVal = max(abs(m_pm->getMaxEigenvalue(uL)),abs(m_pm->getMaxEigenvalue(uR)));
    m_flux_l      = m_pm->computePhysicalFlux(uL);
    // m_eigenvals_l = m_pm->getEigenvalues(uL);
    m_flux_r      = m_pm->computePhysicalFlux(uR);
    // m_eigenvals_r = m_pm->getEigenvalues(uR);

    for (unsigned int iEq = 0; iEq < NBEQS; iEq++){
        
        m_flux[iEq] = 0.5*(m_flux_l[iEq] + m_flux_r[iEq]) - maxEigenVal*0.5*(uR[iEq] - uL[iEq]);
    }

    return m_flux;

}

void LaxFriedrichFlux::fluxWithEigenvalues(const CellDataRef uL, const CellDataRef uR, double lambdaL, double lambdaR, double* F)
{
    const double maxEigenVal = max(lambdaL, lambdaR);
    m_flux_l      = m_pm->computePhysicalFlux(uL);
    m_flux_r      = m_pm->computePhysicalFlux(uR);
    for (unsigned int iEq = 0; iEq < NBEQS; iEq++){
        F[iEq] = 0.5*(m_flux_l[iEq] + m_flux_r[iEq]) - maxEigenVal*0.5*(uR[iEq] - uL[iEq]);
    }
}



// The first line of operator(): the two calls are the arguments of max, evaluated in the order the compiler gives to
// the arguments of a call, the same in both functions.
void LaxFriedrichFlux::eigenvalueCalls(const CellDataRef uL, const CellDataRef uR, bool callL, bool callR)
{
    (void) max(callL ? abs(m_pm->getMaxEigenvalue(uL)) : 0., callR ? abs(m_pm->getMaxEigenvalue(uR)) : 0.);
}
