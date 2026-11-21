//
//  ForwardEuler.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 17/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "ForwardEuler.hpp"

// Self-register with the factory
#include "TimeMethodRegistrar.hpp"
REGISTER_TIMEMETHOD("ForwardEuler", ForwardEuler);
#include "../MeshData/MeshData.hpp"
#include "../Parallel.hpp"
#include "../MeshData/Cell1D.cpp"
#include <iostream>
#include <ctime>
#include <cmath>
#include <pybind11/pybind11.h>
#include <mpi.h>

void ForwardEuler::takeStep(double dt)
{
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    vector<Cell1D>& cells = MeshData::getInstance().getData<Cell1D>("Cells");
    // TODO: See if I can store the CC values in a reference vector
    vector<double>& rhs = MeshData::getInstance().getData<double>("rhs");
    double* source = MeshData::getInstance().get2DData<double>("source").mutable_data(0,0);
    //// old boundary call = MeshData::getInstance().getData<double>("boundaries");
    
    double& physTime = MeshData::getInstance().getData<double>("physTime")[0];
    
    m_iter++;


    if(my_rank == MPI_WRITER){
        if(m_iter % PRINTRATE == 0)
            py::print("dt = ", dt);
    }
    
    // dt/dx of a cell is computed once for its equations, and the residuals of the equations are summed together,
    // cell after cell: each sum takes its terms in the same order as the loop equation by equation, so the norm
    // is the same to the bit. With several threads the cells are updated in parallel, their squared residuals
    // kept and summed afterwards in this order.
    for(unsigned int iEq = 0; iEq < NBEQS; ++iEq)
        m_norm[iEq] = 0;
    const double dx0 = cells[0].dx;
    auto update = [&](int begin, int end, double* r2) {
        for (int iCell = begin; iCell < end; ++iCell) {
            const double k = dt/cells[iCell].dx;
            for(unsigned int iEq = 0; iEq < NBEQS; ++iEq){
                const double rhsU    = rhs[iEq*NBCELLS + iCell];
                const double S_i     = source[iEq*NBCELLS + iCell];
                cells[iCell].uCC[iEq] = cells[iCell].uCC[iEq] - k*rhsU + S_i*dt;
                const double Residual = S_i*dx0 - rhsU;
                if (r2) r2[iCell*NBEQS + iEq] = Residual*Residual;
                else m_norm[iEq] += Residual*Residual;
            }
        }
    };
    if (muffinThreads() > 1) {
        m_residual2.resize(NBEQS*NBCELLS);
        double* r2 = m_residual2.data();
        parallelRange(0, NBCELLS, [&](int begin, int end, int) {update(begin, end, r2);});
        for (int iCell = 0; iCell < NBCELLS; ++iCell)
            for(unsigned int iEq = 0; iEq < NBEQS; ++iEq)
                m_norm[iEq] += r2[iCell*NBEQS + iEq];
    }
    else update(0, NBCELLS, nullptr);

    physTime += dt;

}
