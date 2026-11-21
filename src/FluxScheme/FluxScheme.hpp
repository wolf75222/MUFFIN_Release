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
    // True if eigenvalueCalls gives the calls of getMaxEigenvalue made by operator(), in the same order.
    virtual bool hasEigenvalueCalls() const {return false;}
    // Calls getMaxEigenvalue on the states of operator() for which call is true, in the order of operator(), so
    // that the warnings printed by the model come in the same order.
    virtual void eigenvalueCalls(const CellDataRef uL, const CellDataRef uR, bool callL, bool callR) {}
    // A copy for another thread, working with the model pm (see PhysicalModel::cloneForThread); nullptr if the
    // scheme keeps state between calls.
    virtual std::unique_ptr<FluxScheme> cloneForThread(PhysicalModel* pm) const {return nullptr;}
protected:
    PhysicalModel* m_pm;
};

#endif /* FluxScheme_hpp */
