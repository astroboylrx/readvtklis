//
//  global.h
//  readvtklis
//
//  Created by Rixin Li on 2/26/15.
//  Copyright (c) 2015 Rixin Li. All rights reserved.
//

#ifndef __readvtklis__global__
#define __readvtklis__global__

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
#include "mpi.h"

using namespace::std;

/*************************************/
/*****Pre-defined Compiler Macros*****/
/*************************************/

#define ENABLE_MPI

//#define RESERVE_PUSH_BACK
//#define FROM_ARRAY_TO_VECTOR
#define RESIZE_LIST

#define OutflowRate
//#define PeriodicFlux

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
    float *Otime;                                  /*!< orbit time */
    long *N_par;                                    /*!< number of par */
    float *Max_Rhop;                               /*!< max partical density */
    float *RpAV;                                   /*!< <rho_p> */
    float *RpSQ;                                   /*!< <rho_p^2>^0.5 */
    float *RpQU;                                   /*!< <rho_p^4>^0.25 */
    float *Hp;                                     /*!< par scale height */
    float *Hp_in1sigma;                            /*!< par scale height derived from Gaussian 1sigma range */
    float *dSigma;                                 /*!< the change of gas surface density due to outflow */
    
    // output to result_MeanSigma.txt, first sigma_g and then sigma_p
    float **MeanSigma;                             /*!< sigma_g and sigma_p averaged over y */
    
    // output to result_VpecG.txt, first x, then y, and then z
    float **VpecG;                                 /*!< Vpec_g averaged horizontally at each z, weighted by rho_g */
    
    // output to result_VertRho.txt, first rho_g and then rho_p
    float **VertRho;                               /*!< Vertical Structure of rho_g and rho_p */
    
    float **CorrL;                                 /*!< Correlation Length  */
#ifdef CorrValue
    float **CorrV;                                 /*!< Correlation Value */
#endif
    
    // output to result_RMPL.txt
    float *RMPL;                                    /*!< Rhop_Max Per Level */
    
    /******relative to calculation********/
    float ****V_gas_0;                              /*!< initial v_gas */
    int dimensions[3];                              /*!< the number of cells in each dimension */
    float spacing[3];                              /*!< the spacing of coordinate */
    float *ccx, *ccy, *ccz;                        /*!< cell center coordinates for plots */
    
    
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
