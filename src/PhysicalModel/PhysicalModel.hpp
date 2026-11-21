//
//  PhysicalModel.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 19/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef PhysicalModel_hpp
#define PhysicalModel_hpp

#include <stdio.h>
#include <memory>
#include <vector>
#include "../Component.hpp"
#include "../CollisionalData/Mixture.hpp"
#include "../DiffusiveFlux/DiffusiveFlux.hpp"


using namespace Parameters;
using namespace std;

class PhysicalModel : public Component
{
public :
    PhysicalModel(string name) : Component(name) {}
    virtual ~PhysicalModel(){};
    void setMixture(Mixture* mixture) {m_mixture = mixture;}
    void setDiffusiveFlux(DiffusiveFlux* diffFlux) {m_diffFlux = diffFlux;}
    Mixture* getMixture(){return m_mixture;}
    virtual vector<double>& getEigenvalues(const CellDataRef u) = 0;
    virtual double getMaxEigenvalue(const CellDataRef u) = 0;
    virtual vector<double>& computePhysicalFlux(const CellDataRef u) = 0;

    // getMaxEigenvalue without side effect, for FVM1D to compute it once per state. A model that has it
    // returns true in hasQuietEigenvalue; getMaxEigenvalueQuiet returns false for a state where
    // getMaxEigenvalue would print, and FVM1D then runs the step as before, prints included.
    virtual bool hasQuietEigenvalue() const {return false;}
    virtual bool getMaxEigenvalueQuiet(const CellDataRef u, double& lambda) {lambda = getMaxEigenvalue(u); return true;}
    // A copy for another thread (MUFFIN_THREADS > 1), for the models whose methods change nothing but their own
    // work arrays; nullptr for the others, which then run on one thread.
    virtual std::unique_ptr<PhysicalModel> cloneForThread() const {return nullptr;}
    
protected :
    Mixture* m_mixture;
    DiffusiveFlux* m_diffFlux;
    
    
};

#endif /* PhysicalModel_hpp */
