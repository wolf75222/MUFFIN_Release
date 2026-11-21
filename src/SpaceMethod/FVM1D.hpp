//
//  FVM1D.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 16/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef FVM1D_hpp
#define FVM1D_hpp

#include <stdio.h>
#include "SpaceMethod.hpp"
#include "../FluxScheme/FluxScheme.hpp"
#include "../SourceTerm/SourceTerm.hpp"
#include "../SpaceReconstructor/SpaceReconstructor.hpp"
#include "../BoundaryCondition/BoundaryCondition.hpp"
#include "../Parameters.hpp"
#include "../MeshData/Cell1D.hpp"
#include <mpi.h>

using namespace Parameters;

class FVM1D : public SpaceMethod
{
public:
    FVM1D(string name) : SpaceMethod(name), m_Fip12(NBEQS), m_Fim12(NBEQS) {}
    ~FVM1D(){}
    void setup();
    void unsetup();
    virtual void setBoundaries();
    virtual void computeRHS();
    void computeRHSReference();             // fluxes, source and time step as in the original code
    virtual void computeSource();
    virtual void computeEField();
    virtual void computeFluxes();
    virtual double getDtOvDx(){
        int world_size, world_rank;
        // We take the min dtOvdx of all the processors
        double globaldtOvdx;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        MPI_Allreduce(&m_dtOvdx, &globaldtOvdx, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

        return globaldtOvdx;
    };
    virtual double getDt(){
        int world_size, world_rank;
        // We take the min dtOvdx of all the processors
        double globaldt;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        MPI_Allreduce(&m_dt, &globaldt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        
        return globaldt;
    };
    
protected:
    unique_ptr<FluxScheme> m_flux;
    unique_ptr<SourceTerm> m_source;
    unique_ptr<SourceTerm> m_EFieldSource;  // Source term for electric field fractional splitting
    unique_ptr<SpaceReconstructor> m_reconstructor;
    unique_ptr<BoundaryCondition> m_InletBC;
    unique_ptr<BoundaryCondition> m_OutletBC;
    CellDataRef m_uInlet;
    CellDataRef m_uOutlet;
    vector<double> m_Fip12;
    vector<double> m_Fim12;

    bool quietEigenvalues();
    void eigenvalueCallsBeforeSource();
    void eigenvalueCallsAfterSource();
    void applySourceFrequency();
    vector<Cell1D>* m_cells = nullptr;
    vector<double>* m_rhs = nullptr;
    vector<double> m_lambdaL;               // |largest eigenvalue| of the left and right states of each cell
    vector<double> m_lambdaR;
    double m_lambdaIn = 0., m_lambdaOut = 0.;
    vector<char> m_warnL;                   // true for the states for which getMaxEigenvalue prints a warning
    vector<char> m_warnR;
    bool m_warnIn = false, m_warnOut = false;
    vector<double> m_faceFlux;              // flux of face i at i*NBEQS, faces 0 to NBCELLS
};

#endif /* FVM1D_hpp */
