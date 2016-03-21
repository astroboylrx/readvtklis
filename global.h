//
//  global.h
//  readvtklis
//
//  Created by Rixin Li on 2/26/15.
//  Copyright (c) 2015 Rixin Li. All rights reserved.
//

#ifndef __readvtklis__global__
#define __readvtklis__global__

#define ENABLE_MPI

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <unistd.h>
#include <getopt.h>
#include <ctime>
#include <cassert>
#include <algorithm>
#ifdef ENABLE_MPI
#include "mpi.h"
#endif

using namespace::std;

/*************************************/
/*****Pre-defined Compiler Macros*****/
/*************************************/

//#define RESERVE_PUSH_BACK
//#define FROM_ARRAY_TO_VECTOR
#define RESIZE_LIST

#define OutflowRate
//#define PeriodicFlux

#define OCTREE
#define QUADTREE
#define BTREE

//#define CorrValue // print correlation value of Mx

/*************************************/
/********Forward Declaration**********/
/*************************************/

// you can forward declaration here, which means if you need global declaration of an function,
// you can declare it here and then define it in another file, also include global.h there.
// But for global object (i.e. class), forward declaration may result in errors like "Member
// access into incomplete type 'XXX'", because definition is needed to access the method of
// the object.

template<typename T>
/*! \fn T ***allocate3d_scalar_array(int *dimensions)
 *  \brief allocate memory to 3d scalar array */
T ***allocate3d_scalar_array(int *dimensions)
{
    // store in x, then y, then stack them as z direction
    T ***data = new T**[dimensions[2]];
    for (int i = 0; i != dimensions[2]; i++) {
        data[i] = new T*[dimensions[1]];
        for (int j = 0; j != dimensions[1]; j++) {
            data[i][j] = new T[dimensions[0]];
            for (int k = 0; k != dimensions[0]; k++) {
                data[i][j][k] = 0;
            }
        }
    }
    return data;
}

template<typename T>
/*! \fn int ***deallocate3d_scalar_array(T ***data, int *dimensions)
 *  \brief deallocate memory of 3d scalar array */
int deallocate3d_scalar_array(T ***data, int *dimensions)
{
    for (int i = 0; i != dimensions[2]; i++) {
        for (int j = 0; j != dimensions[1]; j++) {
            delete [] data[i][j];
        }
        delete [] data[i];
    }
    delete [] data;
    return 0;
}

template<typename T>
/*! \fn T ****allocate3d_vector_array(int *dimensions)
 *  \brief allocate memory to 3d vector array */
T ****allocate3d_vector_array(int *dimensions)
{
    // store in x, then y, then stack them as z direction
    T ****data = new T***[dimensions[2]];
    for (int i = 0; i != dimensions[2]; i++) {
        data[i] = new T**[dimensions[1]];
        for (int j = 0; j != dimensions[1]; j++) {
            data[i][j] = new T*[dimensions[0]];
            for (int k = 0; k != dimensions[0]; k++) {
                data[i][j][k] = new T[3];
                data[i][j][k][0] = 0;
                data[i][j][k][1] = 0;
                data[i][j][k][2] = 0;
            }
        }
    }
    return data;
}

template<typename T>
/*! \fn int ****deallocate3d_vcetor_array(T ****data, int *dimensions)
 *  \brief deallocate memory of 3d vector array */
int deallocate3d_vector_array(T ****data, int *dimensions)
{
    for (int i = 0; i != dimensions[2]; i++) {
        for (int j = 0; j != dimensions[1]; j++) {
            for (int k = 0; k != dimensions[0]; k++) {
                delete [] data[i][j][k];
            }
            delete [] data[i][j];
        }
        delete [] data[i];
    }
    delete [] data;
    return 0;
}

/*! \fn string pvector(T x[3])
 *  \brief print vector */
template<int D, typename T>
string pvector(T x[])
{
    ostringstream oss;
    oss << "[";
    for (int i = 0; i != D; i++) {
        if (i != 0) {
            oss << ",";
        }
        oss << x[i];
    }
    oss << "]";
    return oss.str();
}

/*************************************/
/***********Dividing Line*************/
/*************************************/

/*! \class Paras2probe
 *  \brief parameters to process for the purpose of output
 */
class Paras2probe {
private:
    
public:
    /*********relative to output**********/
    // rule of thumb: if multiple parameters are written to a single file, then gas properties always go first
    
    // output to result_Par.txt (in column-order)
    float *Otime;                                   /*!< orbit time */
    long *N_par;                                    /*!< number of par */
    float *Max_Rhop;                                /*!< max partical density */
    double *RpAV;                                   /*!< <rho_p> */
    double *RpSQ;                                   /*!< <rho_p^2>^0.5 */
    double *RpQU;                                   /*!< <rho_p^4>^0.25 */
    double *Hp;                                     /*!< par scale height */
    float *Hp_in1sigma;                             /*!< par scale height derived from Gaussian 1sigma range */
    double *dSigma;                                 /*!< the change of gas surface density due to outflow */
    
    // output to result_MeanSigma.txt, first sigma_g and then sigma_p
    double **MeanSigma;                             /*!< sigma_g and sigma_p averaged over y */
    double **Sigma;                                 /*!< sigma_g and sigma_p */
    
    // output to result_VpecG.txt, first x, then y, and then z
    double **VpecG;                                 /*!< Vpec_g averaged horizontally at each z, weighted by rho_g */
    
    // output to result_VertRho.txt, first rho_g and then rho_p
    double **VertRho;                               /*!< Vertical Structure of rho_g and rho_p */
    
    float **CorrL;                                  /*!< Correlation Length  */
#ifdef CorrValue
    float **CorrV;                                  /*!< Correlation Value */
#endif
    
    // output to result_GasPar.txt
    ////////////////////////////////////////////////////////////////////
    // For volume-averaged properties, gas momentum or kinetic energy,
    // P_tot = SUM_cells(m_gas*F) = dV_cell * SUM_cells(rho_gas*F)
    // P_tot/V = SUM_cells(rho_gas*F) * dV_cell/V
    //         = SUM_cells(rho_gas*F) * dx*dy*dz / (Nx*dx*Ny*dy*Nz*dz)
    //         = SUM_cells(rho_gas*F) / (Nx*Ny*Nz)
    // here F is u_gas or 0.5*u_gas^2, respectively. Similar formula can
    // be applied to the Particle-in-Mesh momentum or kinetic energy.
    ////////////////////////////////////////////////////////////////////
    // For area-averaged properties, gas momentum or kinetic energy,
    // P_tot/A = SUM_cells(rho_gas*F) * dV_cell/A
    //         = SUM_cells(rho_gas*F) * dx*dy*dz / (Nx*dx*Ny*dy)
    //         = SUM_cells(rho_gas*F) * dz / (Nx*Ny)
    // here F is u_gas or 0.5*u_gas^2, respectively. Similar formula can
    // be applied to the Particle-in-Mesh momentum or kinetic energy.
    ////////////////////////////////////////////////////////////////////
    // There is no point to calculate volume-averaged or area-averaged
    // velocity. Since what is averaged is the mass, velocity is the
    // weight.
    ////////////////////////////////////////////////////////////////////
    // For the particle properties computed directly from all the
    // particles, note that you need to calculate the mass of one single
    // particle in the simulations.
    // Sigma_gas = Integrate[Exp[-z^2/2], {z, -Infinity, Infinity}]
    //           = Sigma_gas_0 = sqrt(2*Pi) = 2.506628274631000502415765
    // m_par = Z * Sigma_gas * A / N_par (if only one species)
    // So,
    // V_par = SUM_pars(v_par)/N_par
    // P_tot/V = SUM_pars(m_par*v_par)/V
    // P_tot/A = SUM_pars(m_par*v_par)/A
    // Ek_par/V = SUM_pars(0.5*m_par*v_par*v_par)/V
    // Ek_par/A = SUM_pars(0.5*m_par*v_par*v_par)/A
    ////////////////////////////////////////////////////////////////////
    // Consider a component of momentum density M = rho*u of either gas
    // or particles. Break into 3 components:
    // M = M_0 + M_1 + M_2, which represent
    // 0: the space and time-averaged motion, i.e. drift;
    // 1: space-averaged motion corrected for drift, i.e. epicycles;
    // 2: ﬂuctuations with vanishing spatial and time averages, i.e.
    //    waves and turbulence.
    // So,
    // M_0 = <M>_VT
    // M_1 = <M>_V - M_0
    // M_2 = M - M_0 - M_1
    // From this point of view, the kinetic energy density
    // e = M^2/(2*rho) = SUM_i(e_i) + SUM_i(SUM_j(M_i*M_j/(2*rho))) i!=j
    // Do we expect the cross terms to be significant? Since <M_2>_V = 0
    // by definition, <M_2/rho>_V would only vanish for a perfectly
    // incompressible fluid. M_0 and M_1 is already spatial averaged, so
    // <e_cross>_V = M_0*<M_2/rho>_V + M_1*<M_2/rho>_V
    // If we want to understand how much energy was contributed by the
    // epicyclic motion we could consider
    // e_1,tot = e_1 + M_1*(M_0+M_2)/rho
    ////////////////////////////////////////////////////////////////////
    
    double **GasHst;                               /*!< volume-averaged gas properties, in the order of p_gas[x,y,z,tot]/V, Ek_gas[x,y,z,tot]/V, then area-averaged gas properties, in the order of p_gas[x,y,z,tot]/A, Ek_gas[x,y,z,tot]/A */
    double **GasHst2;                              /*!< volume-averaged gas properties inside z = +/- 0.1 (the overlap region with the smallest box), in the order of p_gas[x,y,z,tot]/V, Ek_gas[x,y,z,tot]/V, then area-averaged gas properties, in the order of p_gas[x,y,z,tot]/A, Ek_gas[x,y,z,tot]/A */
    double **ParHst;                               /*!< volume-averaged particle properties, in the order of p_par[x,y,z,tot]/V, Ek_par[x,y,z,tot]/V, then area-averaged particle properties, in the order of p_par[x,y,z,tot]/A, Ek_par[x,y,z,tot]/A */
    double **ParLis;                               /*!< particle dynamical properties directly from all the particles, in the order of SUM(v_par[x,y,z,tot])/Npar, SUM(p_par[x,y,z,tot])/V, SUM(p_par[x,y,z,tot])/A, SUM(Ek_par[x,y,z,tot])/V, SUM(Ek_par[x,y,z,tot])/A */
    double **GPME;                                 /*!< Investigate the budget of Gas/Particle Momentum/Energy for Gas, in the below order (<...>_V indicates the average is from definitions, .../V means we average it on purpose.) */
        /* In the output file, it looks like
         between each two dash lines are (Nz+1) lines of data, N_row = n_file
         ---------------------------------------------------------------------
         <M_0>_VT(x,y,z,tot],           <M_0>_VT[x,y,z,tot] per orbit,
         <M_0(z)>_VT[x,y,z,tot] per dz, <M_0(z)>_VT[x,y,z,tot] per orbit
         ---------------------------------------------------------------------
         <M_1(t)>_V[x,y,z,tot]
         <M_1(z,t)>_V[x,y,z,tot] per dz
         ---------------------------------------------------------------------
         std of M_2(t)[x,y,z,tot]/V
         std of M_2(z,t)[x,y,z,tot]/A
         ---------------------------------------------------------------------
         <e_0>_VT[x,y,z,tot],           <e_0>_VT[x,y,z,tot] per orbit,
         <e_0(z)>_VT[x,y,z,tot] per dz, <e_0(z,t)>_VT[x,y,z,tot] per orbit
         ---------------------------------------------------------------------
         <e_1(t)>_V[x,y,z,tot]
         <e_1(z,t)>_V[x,y,z,tot] per dz
         ---------------------------------------------------------------------
         e_2(t)[x,y,z,tot]/V
         e_2(z,t)[x,y,z,tot]/A
         ---------------------------------------------------------------------
         <M_0*M_1/rho>_V(t)[x,y,z,tot]
         <M_0*M_1/rho>_V(z,t)[x,y,z,tot] per dz
         ---------------------------------------------------------------------
         M_0*M_2/rho(t)[x,y,z,tot]/V
         M_0*M_2/rho(z,t)[x,y,z,tot]/A
         ---------------------------------------------------------------------
         M_1*M_2/rho(t)[x,y,z,tot]/V
         M_1*M_2/rho(z,t)[x,y,z,tot]/A
         ---------------------------------------------------------------------
         Note: 
         a, M_0 and e_0 won't consume all the rows, only 1+ceil(n_file/8Pi)
         But in program we use [i][j] to denote [time][...]->N_column=n_file
         b, M_2/V is 0.0 by the definition, so the std of M2/V is used instead
         */
    
    double **GPMEPar;                                /*!< Investigate the budget of Gas/Particle Momentum/Energy for Par, in the order above (<...>_V indicates the average is from definitions, .../V means we average it on purpose.)  */
    
    // output to result_RMPL.txt
    float *RMPL;                                    /*!< Rhop_Max Per Level */
    
    /******relative to calculation********/
    int dimensions[3];                              /*!< the number of cells in each dimension */
    float spacing[3];                               /*!< the spacing of coordinate */
    float *ccx, *ccy, *ccz;                         /*!< cell center coordinates for plots */
    
    /*! \fn int AllocateMemory(int n_file)
     *  \brief allocate memory to parameters (first level)
     */
    int AllocateMemory(int n_file);
    
    /*! \fn int AllocateSubMemory(int n_file, int *dimensions)
     *  \brief allocate memory to parameters (second level)
     */
    int AllocateSubMemory(int n_file, int *dimensions);
    
    /*! \fn ~Paras2probe()
     *  \brief deallocate memory of parameters
     */
    ~Paras2probe();
    
};

#ifdef ENABLE_MPI
/*! \class MPI_info
 *  \brief wrapper of MPI routines
 */
class MPI_info {
private:
    
public:
    int numprocs;                                   /*!< number of processor */
    int myrank;                                     /*!< rank of cpu */
    int master;                                     /*!< rank of master cpu */
    int loop_begin, loop_end, loop_offset;          /*!< begin/end/offset of loop */
    double wait_time;                               /*!< time wasting on waiting */
    MPI::Status status;                             /*!< status of recv function */
    
    Paras2probe paras;                              /*!< contains parameters for the use of MPIAllreduce */
    
    /*! \fn int Initialize(int argc, const char * argv[])
     *  \brief MPI initializaion */
    int Initialize(int argc, const char * argv[]);
    
    /*! \fn int Determine_Loop(int n_file)
     *  \brief determine the begin/end/offset for loop */
    int Determine_Loop(int n_file);
    
    /*! \fn int Barrier()
     *  \brief wrapper of MPI Barrier */
    int Barrier();
    
    /*! \fn double T()
     *  \brief wrapper of MPI WTime */
    double T();
    
    /*! \fn string prank()
     *  \brief return string of Processor myrank: */
    string prank();
    
    /*! \fn int Finalize()
     *  \brief wrapper of MPI Finalize() */
    int Finalize();
    
    /*! \fn ~MPI_info()
     *  \brief destructor */
    ~MPI_info();
    
};

extern MPI_info *myMPI; // declaration is in main function
#endif

#endif /* defined(__readvtklis__global__) */
