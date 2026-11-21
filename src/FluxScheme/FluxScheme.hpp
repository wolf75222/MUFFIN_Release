//
//  FluxScheme.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 11/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef FluxScheme_hpp
#define FluxScheme_hpp

#include <stdio.h>
#include <algorithm>
#include "../Component.hpp"
#include "../PhysicalModel/PhysicalModel.hpp"

class FluxScheme : public Component
{
public:
    FluxScheme(string name) : Component(name) {}
    virtual ~FluxScheme() {}
    virtual void setPhysicalModel(PhysicalModel* pm){m_pm = pm;}
    virtual vector<double>& operator()(const CellDataRef uL, const CellDataRef uR) = 0;
    // The same flux written to F, given lambdaL = |getMaxEigenvalue(uL)| and lambdaR = |getMaxEigenvalue(uR)|.
    // The schemes that use these values (LaxFriedrich, SixMomentAQMOM1DFlux) do not compute them again.
    virtual void fluxWithEigenvalues(const CellDataRef uL, const CellDataRef uR, double lambdaL, double lambdaR, double* F)
    {
        const vector<double>& f = (*this)(uL, uR);
        std::copy(f.begin(), f.end(), F);
    }
protected:
    PhysicalModel* m_pm;
};

#endif /* FluxScheme_hpp */
