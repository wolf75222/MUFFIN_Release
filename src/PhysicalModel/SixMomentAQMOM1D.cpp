//
//  SixMomentAQMOM1D.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 26/09/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "SixMomentAQMOM1D.hpp"
#include "QMOMMaxEigenvalue.hpp"
#include "ModelWarnings.hpp"
#include "PhysicalModelRegistrar.hpp"
#include <vector>

REGISTER_PHYSICALMODEL("SixMomentAQMOM1D", SixMomentAQMOM1D);
#include <cmath>
#include <algorithm>    // std::min_element, std::max_element

using namespace std;


vector<double>& SixMomentAQMOM1D::getEigenvalues(const CellDataRef u)
{
    // We compute just the max
    m_eigenvals[0] = SixMomentAQMOM1D::getMaxEigenvalue(u);
    return m_eigenvals;
}

double SixMomentAQMOM1D::getMaxEigenvalue(const CellDataRef u)
{
    bool realizable;
    double q_nStar, r_nStar;
    const double lambda = qmomMaxEigenvalue(u, realizable, q_nStar, r_nStar);
    if(!realizable)
        ModelWarnings::realizability(q_nStar, r_nStar);
    return lambda;
}

bool SixMomentAQMOM1D::getMaxEigenvalueQuiet(const CellDataRef u, double& lambda)
{
    bool realizable;
    double q_nStar, r_nStar;
    lambda = qmomMaxEigenvalue(u, realizable, q_nStar, r_nStar);
    return realizable;
}

vector<double>& SixMomentAQMOM1D::computePhysicalFlux(const CellDataRef u)
{
    const double rho    = u[0];
    const double v      = u[1]/rho;
    const double p      = u[2] - rho*v*v;
    const double q      = u[3] - rho*v*v*v - 3*v*p;
    const double r      = u[4] - rho*v*v*v*v - 6*v*v*p - 4*v*q;
    const double r_Max  =  3*p*p/rho;
    const double K      = r - r_Max;
    m_physFlux[0] = u[1];
    m_physFlux[1] = u[2];
    m_physFlux[2] = u[3];
    m_physFlux[3] = u[4];
    m_physFlux[4] = rho*v*v*v*v*v + 10.*v*v*v*p + 10.*v*v*q + 5.*v*r;
    
    m_physFlux[4] += 2*r*q/p - pow(q,3)/pow(p,2);
    m_physFlux[5] = v*u[5];
    // cout<<"m_physFlux[5]= "<<  m_physFlux[5] <<"u[5] "<< u[5]<< "\n";
    // cout<<"s = "<< s.cast<double>() <<"\n";
    
    return m_physFlux;
}
