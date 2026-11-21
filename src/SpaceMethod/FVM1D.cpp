//
//  FVM1D.cpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 16/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#include "FVM1D.hpp"
#include "../PhysicalModel/ModelWarnings.hpp"
#include "../Parallel.hpp"
#include "../FluxScheme/FluxSchemeFactory.hpp"
#include "../SourceTerm/SourceTermFactory.hpp"
#include "../SpaceReconstructor/ReconstructorFactory.hpp"
#include "../BoundaryCondition/BoundaryConditionFactory.hpp"
#include "../MeshData/MeshData.hpp"
#include "../MeshData/Cell1D.hpp"
#include "../CollisionalData/CollisionalData.hpp"
#include <pybind11/pybind11.h>
#include <mpi.h>


namespace py = pybind11;


void FVM1D::setup(){
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    if(my_rank == MPI_WRITER)
        py::print("Setting-up Space discretization class:\t ",this->getName(),"\n");
    
    // TODO: Change this

    
    m_cells = &MeshData::getInstance().getData<Cell1D>("Cells");
    m_rhs   = &MeshData::getInstance().getData<double>("rhs");
    m_lambdaL.assign(NBCELLS, 0.);
    m_lambdaR.assign(NBCELLS, 0.);
    m_warnL.assign(NBCELLS, 0);
    m_warnR.assign(NBCELLS, 0);
    m_faceFlux.assign((NBCELLS + 1)*NBEQS, 0.);

    // Initialize the flux scheme
    m_flux = FluxSchemeFactory::CreateFluxScheme(FLUXSCHEMENAME);
    m_flux->setPhysicalModel(m_pm);

    // Model and flux scheme of each thread (MUFFIN_THREADS > 1): copies for threads 1, 2, ..., if the model and the
    // scheme can be copied; otherwise the fluxes and eigenvalues are computed on one thread.
    m_pmOf.assign(1, m_pm);
    m_fluxOf.assign(1, m_flux.get());
    m_pmThread.clear();
    m_fluxThread.clear();
    m_threads = false;
    if (muffinThreads() > 1 && m_pm->hasQuietEigenvalue()) {
        m_threads = true;
        for (int t = 1; t < muffinThreads() && m_threads; ++t) {
            unique_ptr<PhysicalModel> pm = m_pm->cloneForThread();
            unique_ptr<FluxScheme> flux = pm ? m_flux->cloneForThread(pm.get()) : nullptr;
            m_threads = pm && flux;
            if (m_threads) {
                m_pmOf.push_back(pm.get());
                m_fluxOf.push_back(flux.get());
                m_pmThread.push_back(std::move(pm));
                m_fluxThread.push_back(std::move(flux));
            }
        }
        if (!m_threads) {
            m_pmOf.assign(1, m_pm);
            m_fluxOf.assign(1, m_flux.get());
            m_pmThread.clear();
            m_fluxThread.clear();
        }
    }

    
    // Initialize the source term
    m_source = SourceTermFactory::CreateSourceTerm(SOURCETERMNAME);
    m_source->setPhysicalModel(m_pm);
    m_source->setup();
    
    // Initialize the E-field source term if fractional splitting is enabled
    if (EFIELDFRACTIONAL) {
        // Check if the functionEFieldPhi exists in CollisionalData before creating the source term   
        py::print("Initializing E-field source term for fractional splitting...\n");
        m_EFieldSource = SourceTermFactory::CreateSourceTerm("PythonEFieldSourceTerm");
        m_EFieldSource->setPhysicalModel(m_pm);
        m_EFieldSource->setup();
        py::print("E-field source term initialized.\n");

    }
    
    // Initialize the reconstructor for second order solution in space
    m_reconstructor = ReconstructorFactory::CreateReconstructor(RECONSTRUCTORNAME);
    m_reconstructor->setup(); // initialize the limiter
    
    // TODO: Change this to be set from input file
    m_InletBC  = BoundaryConditionFactory::CreateBoundaryCondition(INLETTYPE,"Left");
    
    m_OutletBC  = BoundaryConditionFactory::CreateBoundaryCondition(OUTLETTYPE,"Right");
    
    if(my_rank == MPI_WRITER){
        py::print("Space discretization using the flux scheme:\t ",m_flux->getName(),"\n");
        py::print("Space discretization using the source term:\t ",m_source->getName(),"\n");
        if (EFIELDFRACTIONAL) {
            py::print("Space discretization using E-field source term:\t ",m_EFieldSource->getName(),"\n");
        }
        py::print("Space discretization using the reconstructor:\t ",m_reconstructor->getName(),"\n");
    }
}

void FVM1D::unsetup(){
    m_reconstructor->unsetup();
}

void FVM1D::setBoundaries(){

    vector<Cell1D>& cells       = MeshData::getInstance().getData<Cell1D>("Cells");
    
    // MPI excahnge of data
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Status status;
    
    int left_neigh = (my_rank>0) ? my_rank-1 : MPI_PROC_NULL;
    int right_neigh = (my_rank<(world_size-1)) ? my_rank+1 : MPI_PROC_NULL;
    

    
    // Send and receive the neighbour information
    vector<double> sendbufL=vector<double>(NBEQS);
    vector<double> sendbufR=vector<double>(NBEQS);
    for(unsigned int iEq = 0; iEq < NBEQS; iEq++){
        sendbufL[iEq] = cells[0].uCC[iEq];
        sendbufR[iEq] = cells[NBCELLS - 1].uCC[iEq];
    }

    MPI_Sendrecv(sendbufL.data(), sendbufL.size(), MPI_DOUBLE, left_neigh, 6, UOUTLET.data(),  UOUTLET.size(), MPI_DOUBLE, right_neigh, 6, MPI_COMM_WORLD, &status);
    MPI_Sendrecv(sendbufR.data(), sendbufR.size(), MPI_DOUBLE, right_neigh, 6, UINLET.data(),  UINLET.size(), MPI_DOUBLE, left_neigh, 6, MPI_COMM_WORLD, &status);
    
    // Here I set-up the m_uInlet
    m_uInlet  = m_InletBC->setBoundary();
    m_uOutlet = m_OutletBC->setBoundary();


// Example of comment to debug the boundaries
//    for (int iP = 0; iP < world_size; ++iP){
//        MPI_Barrier(MPI_COMM_WORLD);
//        if (my_rank == iP) {
//            cout<<"In rank = "<<my_rank<<" boundaries[0] = "<< boundaries[0] <<" boundaries[1] = "<< boundaries[1] <<"\n";
//            cout<<"In rank = "<<my_rank<<" boundaries[2] = "<< boundaries[2] <<" boundaries[3] = "<< boundaries[3] <<"\n";
//            cout<<"In rank = "<<my_rank<<" boundaries[4] = "<< boundaries[4] <<" boundaries[5] = "<< boundaries[5] <<"\n";}
//        MPI_Barrier(MPI_COMM_WORLD);
//    }
}

void FVM1D::computeRHS(){

    // Reconstruction
    m_reconstructor->reconstructField();

    // Each face flux and each eigenvalue once, the same values as computeRHSReference, which computes the
    // flux of a face twice and the eigenvalue of a state four times. When states need the warning of the model,
    // getMaxEigenvalue runs again on them in the order of computeRHSReference, so that the log is the same.
    // computeRHSReference runs instead when the model has no quiet eigenvalue, or when a state needs the warning
    // and the flux scheme does not give the order of its calls.
    if (!m_pm->hasQuietEigenvalue()) {
        computeRHSReference();
        return;
    }
    const bool warnings = !quietEigenvalues();
    if (warnings && !m_flux->hasEigenvalueCalls()) {
        computeRHSReference();
        return;
    }

    const int    n1 = NBCELLS - 1;          // number of cells minus 1
    SourceTerm& sourceterm      = *m_source;
    vector<Cell1D>& cells       = *m_cells;
    vector<double>& rhs         = *m_rhs;
    const double* lL = m_lambdaL.data();
    const double* lR = m_lambdaR.data();
    double* F = m_faceFlux.data();

    // Face i is between cells i - 1 and i; face 0 has the inlet on its left, face NBCELLS the outlet on its right
    auto faces = [&](int begin, int end, int thread) {
        FluxScheme& flux1D = *m_fluxOf[thread];
        for (int i = begin; i < end; ++i) {
            if (i == 0)
                flux1D.fluxWithEigenvalues(m_uInlet, cells[0].uL, m_lambdaIn, lL[0], F);
            else if (i == NBCELLS)
                flux1D.fluxWithEigenvalues(cells[n1].uR, m_uOutlet, lR[n1], m_lambdaOut, F + NBCELLS*NBEQS);
            else
                flux1D.fluxWithEigenvalues(cells[i-1].uR, cells[i].uL, lR[i-1], lL[i], F + i*NBEQS);
        }
    };
    if (m_threads) parallelRange(0, NBCELLS + 1, faces);
    else faces(0, NBCELLS + 1, 0);

    for (unsigned int iEq = 0; iEq < NBEQS; iEq++)
        for (int iCell = 0; iCell < NBCELLS; ++iCell)
            rhs[iEq*NBCELLS + iCell] = F[(iCell+1)*NBEQS + iEq] - F[iCell*NBEQS + iEq];

    // Time step: the sequence of computeRHSReference (running max of the eigenvalues, running min of dt),
    // with the eigenvalues of the states already computed
    double maxEigenvalue = 0;
    maxEigenvalue = max(maxEigenvalue, max(lR[0],max(lL[1], max(m_lambdaIn, lL[0]))));
    m_dtOvdx = CFL/maxEigenvalue;
    m_dt     = m_dtOvdx*cells[0].dx;
    for (int iCell = 1; iCell < n1; ++iCell) {
        maxEigenvalue    = max(maxEigenvalue, max(lR[iCell],max(lL[iCell+1], max(lR[iCell-1], lL[iCell]))));
        m_dtOvdx = min(m_dtOvdx, CFL/maxEigenvalue);
        m_dt     = min(m_dt, m_dtOvdx*cells[iCell].dx);
    }

    if (warnings) eigenvalueCallsBeforeSource();

    // Compute the sourceTerm
    sourceterm.computeSource();

    if (warnings) eigenvalueCallsAfterSource();

    maxEigenvalue    = max(maxEigenvalue, max(lR[n1],max(m_lambdaOut, max(lR[n1-1], lL[n1]))));
    m_dtOvdx = min(m_dtOvdx, CFL/maxEigenvalue);
    m_dt     = min(m_dt, m_dtOvdx*cells[n1].dx);

    applySourceFrequency();
}

// |largest eigenvalue| of every state of the faces, without side effect; false if a state needs the warning
bool FVM1D::quietEigenvalues(){
    vector<Cell1D>& cells = *m_cells;
    double lambda;
    m_warnIn = !m_pm->getMaxEigenvalueQuiet(m_uInlet, lambda);
    m_lambdaIn = abs(lambda);
    m_warnOut = !m_pm->getMaxEigenvalueQuiet(m_uOutlet, lambda);
    m_lambdaOut = abs(lambda);
    auto states = [&](int begin, int end, int thread) {
        PhysicalModel& pm = *m_pmOf[thread];
        double lambda;
        for (int iCell = begin; iCell < end; ++iCell) {
            m_warnL[iCell] = !pm.getMaxEigenvalueQuiet(cells[iCell].uL, lambda);
            m_lambdaL[iCell] = abs(lambda);
            m_warnR[iCell] = !pm.getMaxEigenvalueQuiet(cells[iCell].uR, lambda);
            m_lambdaR[iCell] = abs(lambda);
        }
    };
    if (m_threads) parallelRange(0, NBCELLS, states);
    else states(0, NBCELLS, 0);
    bool ok = !(m_warnIn || m_warnOut);
    for (int iCell = 0; iCell < NBCELLS; ++iCell)
        ok = ok && !(m_warnL[iCell] || m_warnR[iCell]);
    return ok;
}

// The calls of getMaxEigenvalue of computeRHSReference on the states that need the warning, in the same order:
// for each cell the fluxes of its right and left faces, then its four states. The warnings are written at once
// before the source term, which may print.
void FVM1D::eigenvalueCallsBeforeSource(){
    vector<Cell1D>& c = *m_cells;
    const int n1 = NBCELLS - 1;
    const char* wL = m_warnL.data();
    const char* wR = m_warnR.data();
    auto face = [&](const CellDataRef& uL, bool callL, const CellDataRef& uR, bool callR) {
        if (callL || callR) m_flux->eigenvalueCalls(uL, uR, callL, callR);
    };
    auto state = [&](const CellDataRef& u, bool call) {
        if (call) m_pm->getMaxEigenvalue(u);
    };
    ModelWarnings::hold();
    face(c[0].uR, wR[0], c[1].uL, wL[1]);
    face(m_uInlet, m_warnIn, c[0].uL, wL[0]);
    state(c[0].uR, wR[0]); state(c[1].uL, wL[1]); state(m_uInlet, m_warnIn); state(c[0].uL, wL[0]);
    for (int i = 1; i < n1; ++i) {
        face(c[i].uR, wR[i], c[i+1].uL, wL[i+1]);
        face(c[i-1].uR, wR[i-1], c[i].uL, wL[i]);
        state(c[i].uR, wR[i]); state(c[i+1].uL, wL[i+1]); state(c[i-1].uR, wR[i-1]); state(c[i].uL, wL[i]);
    }
    face(c[n1].uR, wR[n1], m_uOutlet, m_warnOut);
    face(c[n1-1].uR, wR[n1-1], c[n1].uL, wL[n1]);
    ModelWarnings::release();
}

void FVM1D::eigenvalueCallsAfterSource(){
    vector<Cell1D>& c = *m_cells;
    const int n1 = NBCELLS - 1;
    ModelWarnings::hold();
    if (m_warnR[n1]) m_pm->getMaxEigenvalue(c[n1].uR);
    if (m_warnOut) m_pm->getMaxEigenvalue(m_uOutlet);
    if (m_warnR[n1-1]) m_pm->getMaxEigenvalue(c[n1-1].uR);
    if (m_warnL[n1]) m_pm->getMaxEigenvalue(c[n1].uL);
    ModelWarnings::release();
}

void FVM1D::applySourceFrequency(){
    vector<Cell1D>& cells = *m_cells;
    const double sourceFrequency = m_source->getFrequency();
    if(abs(sourceFrequency) < 1e-12){ // case where the frequency is close to zero
        const double dtOvdx_convective   = m_dtOvdx;
        const double dt_convective       = m_dt;
        m_dtOvdx = dtOvdx_convective;
        m_dt     = dt_convective;
    }
    else{
        const double Dx                  = cells[0].dx; //Assuming the cells have the same delta_x
        const double dtOvdx_source       = CFL/(sourceFrequency*Dx);
        const double dtOvdx_convective   = m_dtOvdx;
        const double dt_source           = CFL/(sourceFrequency);
        const double dt_convective       = m_dt;

        m_dtOvdx = min(dtOvdx_source, dtOvdx_convective);
        m_dt     = min(dt_source, dt_convective);
    }
}

void FVM1D::computeRHSReference(){
    
    const int    n1 = NBCELLS - 1;          // number of cells minus 1
    
    FluxScheme& flux1D          = *m_flux;      // reference is used to keep the same syntax as before
    SourceTerm& sourceterm      = *m_source;
    vector<Cell1D>& cells       = MeshData::getInstance().getData<Cell1D>("Cells");
    // TODO: See if I can store the CC values in a reference vector
    vector<double>& rhs         = MeshData::getInstance().getData<double>("rhs");
    
    // copy the boundaries to a value
    //not anymore since it is the same ref
    
    // The reconstruction is done by computeRHS
    double maxEigenvalue = 0;

    // Main loop
    m_Fip12         = flux1D(cells[0].uR, cells[1].uL);
    m_Fim12         = flux1D(m_uInlet, cells[0].uL);
    // First cell
    for (unsigned int iEq = 0; iEq < NBEQS; iEq++){
        // compute the flux
        rhs[iEq*NBCELLS]        = m_Fip12[iEq] - m_Fim12[iEq];

    }
    // compute the max eigenvalue
    double eigenval_iR      = abs(m_pm->getMaxEigenvalue(cells[0].uR));
    double eigenval_ip1L    = abs(m_pm->getMaxEigenvalue(cells[1].uL));
    double eigenval_im1R    = abs(m_pm->getMaxEigenvalue(m_uInlet));
    double eigenval_iL      = abs(m_pm->getMaxEigenvalue(cells[0].uL));
    maxEigenvalue           = max(maxEigenvalue, max(eigenval_iR,max(eigenval_ip1L, max(eigenval_im1R, eigenval_iL))));
    m_dtOvdx = CFL/maxEigenvalue;
    m_dt     = m_dtOvdx*cells[0].dx;
    
    // Loop over inner cells
    for (int iCell = 1; iCell < n1; ++iCell) {
        m_Fip12         = flux1D(cells[iCell].uR, cells[iCell+1].uL);
        m_Fim12         = flux1D(cells[iCell-1].uR, cells[iCell].uL);
        for (unsigned int iEq = 0; iEq < NBEQS; iEq++){
            // compute the flux
            rhs[iEq*NBCELLS + iCell]    = m_Fip12[iEq] - m_Fim12[iEq];
        }
        
        // compute the max eigenvalue
        eigenval_iR      = abs(m_pm->getMaxEigenvalue(cells[iCell].uR));
        eigenval_ip1L    = abs(m_pm->getMaxEigenvalue(cells[iCell+1].uL));
        eigenval_im1R    = abs(m_pm->getMaxEigenvalue(cells[iCell-1].uR));
        eigenval_iL      = abs(m_pm->getMaxEigenvalue(cells[iCell].uL));
        maxEigenvalue    = max(maxEigenvalue, max(eigenval_iR,max(eigenval_ip1L, max(eigenval_im1R, eigenval_iL))));
        m_dtOvdx = min(m_dtOvdx, CFL/maxEigenvalue);
        m_dt     = min(m_dt, m_dtOvdx*cells[iCell].dx);
    } // integral flux on internal cells
    // Last cell
    m_Fip12         = flux1D(cells[n1].uR, m_uOutlet);
    m_Fim12         = flux1D(cells[n1-1].uR, cells[n1].uL);
    for (unsigned int iEq = 0; iEq < NBEQS; iEq++){
        // compute the flux
        rhs[iEq*NBCELLS + n1]       = m_Fip12[iEq] - m_Fim12[iEq];
    }
    
    // Compute the sourceTerm
    sourceterm.computeSource();
    
    // compute the max eigenvalue
    eigenval_iR      = abs(m_pm->getMaxEigenvalue(cells[n1].uR));
    eigenval_ip1L    = abs(m_pm->getMaxEigenvalue(m_uOutlet));
    eigenval_im1R    = abs(m_pm->getMaxEigenvalue(cells[n1-1].uR));
    eigenval_iL      = abs(m_pm->getMaxEigenvalue(cells[n1].uL));
    maxEigenvalue    = max(maxEigenvalue, max(eigenval_iR,max(eigenval_ip1L, max(eigenval_im1R, eigenval_iL))));
    m_dtOvdx = min(m_dtOvdx, CFL/maxEigenvalue);
    m_dt     = min(m_dt, m_dtOvdx*cells[n1].dx);

    const double sourceFrequency = sourceterm.getFrequency();
    if(abs(sourceFrequency) < 1e-12){ // case where the frequency is close to zero
        const double dtOvdx_convective   = m_dtOvdx;
        const double dt_convective       = m_dt;
        m_dtOvdx = dtOvdx_convective;
        m_dt     = dt_convective;
    }
    else{
        const double Dx                  = cells[0].dx; //Assuming the cells have the same delta_x
        const double dtOvdx_source       = CFL/(sourceFrequency*Dx);
        const double dtOvdx_convective   = m_dtOvdx;
        const double dt_source           = CFL/(sourceFrequency);
        const double dt_convective       = m_dt;

        m_dtOvdx = min(dtOvdx_source, dtOvdx_convective);
        m_dt     = min(dt_source, dt_convective);
    }

}

void FVM1D::computeFluxes(){

    const int    n1 = NBCELLS - 1;          // number of cells minus 1
    
    FluxScheme& flux1D          = *m_flux;      // reference is used to keep the same syntax as before
    SourceTerm& sourceterm      = *m_source;
    vector<Cell1D>& cells       = MeshData::getInstance().getData<Cell1D>("Cells");
    // TODO: See if I can store the CC values in a reference vector
    vector<double>& rhs         = MeshData::getInstance().getData<double>("rhs");
    
    // copy the boundaries to a value
    //not anymore since it is the same ref
    
    // Reconstruction
    m_reconstructor->reconstructField();
    
    double maxEigenvalue = 0;

    // Main loop
    m_Fip12         = flux1D(cells[0].uR, cells[1].uL);
    m_Fim12         = flux1D(m_uInlet, cells[0].uL);
    // First cell
    for (unsigned int iEq = 0; iEq < NBEQS; iEq++){
        // compute the flux
        rhs[iEq*NBCELLS]        = m_Fip12[iEq] - m_Fim12[iEq];

    }
    // compute the max eigenvalue
    double eigenval_iR      = abs(m_pm->getMaxEigenvalue(cells[0].uR));
    double eigenval_ip1L    = abs(m_pm->getMaxEigenvalue(cells[1].uL));
    double eigenval_im1R    = abs(m_pm->getMaxEigenvalue(m_uInlet));
    double eigenval_iL      = abs(m_pm->getMaxEigenvalue(cells[0].uL));
    maxEigenvalue           = max(maxEigenvalue, max(eigenval_iR,max(eigenval_ip1L, max(eigenval_im1R, eigenval_iL))));
    m_dtOvdx = CFL/maxEigenvalue;
    m_dt     = m_dtOvdx*cells[0].dx;
    
    // Loop over inner cells
    for (int iCell = 1; iCell < n1; ++iCell) {
        m_Fip12         = flux1D(cells[iCell].uR, cells[iCell+1].uL);
        m_Fim12         = flux1D(cells[iCell-1].uR, cells[iCell].uL);
        for (unsigned int iEq = 0; iEq < NBEQS; iEq++){
            // compute the flux
            rhs[iEq*NBCELLS + iCell]    = m_Fip12[iEq] - m_Fim12[iEq];
        }
        
        // compute the max eigenvalue
        eigenval_iR      = abs(m_pm->getMaxEigenvalue(cells[iCell].uR));
        eigenval_ip1L    = abs(m_pm->getMaxEigenvalue(cells[iCell+1].uL));
        eigenval_im1R    = abs(m_pm->getMaxEigenvalue(cells[iCell-1].uR));
        eigenval_iL      = abs(m_pm->getMaxEigenvalue(cells[iCell].uL));
        maxEigenvalue    = max(maxEigenvalue, max(eigenval_iR,max(eigenval_ip1L, max(eigenval_im1R, eigenval_iL))));
        m_dtOvdx = min(m_dtOvdx, CFL/maxEigenvalue);
        m_dt     = min(m_dt, m_dtOvdx*cells[iCell].dx);
    } // integral flux on internal cells
    // Last cell
    m_Fip12         = flux1D(cells[n1].uR, m_uOutlet);
    m_Fim12         = flux1D(cells[n1-1].uR, cells[n1].uL);
    for (unsigned int iEq = 0; iEq < NBEQS; iEq++){
        // compute the flux
        rhs[iEq*NBCELLS + n1]       = m_Fip12[iEq] - m_Fim12[iEq];
    }
    
    // Compute the sourceTerm
    sourceterm.computeSource();
    
    // compute the max eigenvalue
    eigenval_iR      = abs(m_pm->getMaxEigenvalue(cells[n1].uR));
    eigenval_ip1L    = abs(m_pm->getMaxEigenvalue(m_uOutlet));
    eigenval_im1R    = abs(m_pm->getMaxEigenvalue(cells[n1-1].uR));
    eigenval_iL      = abs(m_pm->getMaxEigenvalue(cells[n1].uL));
    maxEigenvalue    = max(maxEigenvalue, max(eigenval_iR,max(eigenval_ip1L, max(eigenval_im1R, eigenval_iL))));
    m_dtOvdx = min(m_dtOvdx, CFL/maxEigenvalue);
    m_dt     = min(m_dt, m_dtOvdx*cells[n1].dx);

    const double sourceFrequency = sourceterm.getFrequency();
    if(abs(sourceFrequency) < 1e-12){ // case where the frequency is close to zero
        const double dtOvdx_convective   = m_dtOvdx;
        const double dt_convective       = m_dt;
        m_dtOvdx = dtOvdx_convective;
        m_dt     = dt_convective;
    }
    else{
        const double Dx                  = cells[0].dx; //Assuming the cells have the same delta_x
        const double dtOvdx_source       = CFL/(sourceFrequency*Dx);
        const double dtOvdx_convective   = m_dtOvdx;
        const double dt_source           = CFL/(sourceFrequency);
        const double dt_convective       = m_dt;

        m_dtOvdx = min(dtOvdx_source, dtOvdx_convective);
        m_dt     = min(dt_source, dt_convective);
    }

}

void FVM1D::computeSource(){
    SourceTerm& sourceterm      = *m_source;
    
    // Compute the sourceTerm
    sourceterm.computeSource();
    const double sourceFrequency = sourceterm.getFrequency();

}

void FVM1D::computeEField(){
    if (EFIELDFRACTIONAL && m_EFieldSource) {
        SourceTerm& EFieldSourceTerm = *m_EFieldSource;
        
        // Compute the E-field source term
        EFieldSourceTerm.computeSource();
    }
}
