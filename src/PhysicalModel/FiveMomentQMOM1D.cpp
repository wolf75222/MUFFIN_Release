//
//  FiveMomentQMOM1D.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 26/09/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "FiveMomentQMOM1D.hpp"
#include "QMOMMaxEigenvalue.hpp"
#include <vector>
#include <cmath>
#include <algorithm>    // std::min_element, std::max_element

#include "PhysicalModelRegistrar.hpp"

REGISTER_PHYSICALMODEL("FiveMomentQMOM1D", FiveMomentQMOM1D);

using namespace std;


vector<double>& FiveMomentQMOM1D::getEigenvalues(const CellDataRef u)
{
    // We compute just the max
    m_eigenvals[0] = FiveMomentQMOM1D::getMaxEigenvalue(u);
    return m_eigenvals;
}

double FiveMomentQMOM1D::getMaxEigenvalue(const CellDataRef u)
{
    bool realizable;
    double q_nStar, r_nStar;
    const double lambda = qmomMaxEigenvalue(u, realizable, q_nStar, r_nStar);
    if(!realizable){
        py::print("WARNING: You are out of the realizability");
        py::print("QStar = ", q_nStar, "  RStar = ", r_nStar);
    }
    return lambda;
}

bool FiveMomentQMOM1D::getMaxEigenvalueQuiet(const CellDataRef u, double& lambda)
{
    bool realizable;
    double q_nStar, r_nStar;
    lambda = qmomMaxEigenvalue(u, realizable, q_nStar, r_nStar);
    return realizable;
}

vector<double>& FiveMomentQMOM1D::computePhysicalFlux(const CellDataRef u)
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
    // cout<<"s = "<< s.cast<double>() <<"\n";
    
    return m_physFlux;
}
