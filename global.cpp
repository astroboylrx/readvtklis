//
//  global.cpp
//  readvtklis
//
//  Created by Rixin Li on 2/26/15.
//  Copyright (c) 2015 Rixin Li. All rights reserved.
//

#include "global.h"

//template<typename T>
/*! \fn T ***allocate3d_scalar_array(int *dimensions)
 *  \brief allocate memory to 3d scalar array
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
}*/

//template<typename T>
/*! \fn int deallocate3d_scalar_array(T ***data, int *dimensions)
 *  \brief deallocate memory of 3d scalar array
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
}*/
    
//template<typename T>
/*! \fn T ****allocate3d_vector_array(int *dimensions)
 *  \brief allocate memory to 3d vector array
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
}*/

//template<typename T>
/*! \fn int deallocate3d_vcetor_array(T ****data, int *dimensions)
 *  \brief deallocate memory of 3d vector array
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
}*/

/*************************************/
/************Paras2probe**************/
/*************************************/

/*! \fn int Initialize()
 *  \brief allocate memory to parameters (first level)
 */
int Paras2probe::AllocateMemory(int n_file)
{
    Otime = new double[n_file];
    N_par = new long[n_file];
    Max_Rhop = new double[n_file];
    Hp = new double[n_file];
    Hp_in1sigma = new double[n_file];
    dSigma = new double[n_file];
    for (int i = 0; i != n_file; i++) {
        // initialize if you don't assign all of them values but use them for calculation eventually
        Otime[i] = 0;
        N_par[i] = 0;
        Max_Rhop[i] = 0;
        Hp[i] = 0;
        Hp_in1sigma[i] = 0;
        dSigma[i] = 0;
    }
    
    MeanSigma = new double*[n_file];
    VpecG = new double*[n_file];
    VertRho = new double*[n_file];
    CorrL = new double*[n_file];
    return 0;
}

/*! \fn int AllocateSubMemory(int n_file, int *dimensions)
 *  \brief allocate memory to parameters (second level)
 */
int Paras2probe::AllocateSubMemory(int n_file, int *dimensions)
{
    for (int i = 0; i != n_file; i++) {
        MeanSigma[i] = new double[2*dimensions[0]];
        for (int j = 0; j != 2*dimensions[0]; j++) {
            MeanSigma[i][j] = 0;
        }
        VpecG[i] = new double[3*dimensions[2]];
        for (int j = 0; j != 3*dimensions[2]; j++) {
            VpecG[i][j] = 0;
        }
        VertRho[i] = new double[2*dimensions[2]];
        for (int j = 0; j != 2*dimensions[2]; j++) {
            VertRho[i][j] = 0;
        }
        CorrL[i] = new double[dimensions[2]];
        for (int j = 0; j != dimensions[2]; j++) {
            CorrL[i][j] = 0;
        }
    }
    return 0;
}

/*! \fn ~Paras2probe()
 *  \brief deallocate memory of parameters
 */
Paras2probe::~Paras2probe()
{
    delete [] Otime;
    delete [] N_par;
    delete [] Max_Rhop;
    delete [] Hp;
    delete [] Hp_in1sigma;
    delete [] dSigma;
    delete [] MeanSigma;
    delete [] VpecG;
    delete [] VertRho;
    delete [] CorrL;
}


#ifdef ENABLE_MPI

/*************************************/
/************MPI_info*****************/
/*************************************/

/********** Initialization **********/
/*! \fn int Initialize(int argc, const char * argv[])
 *  \brief MPI initializaion */
int MPI_info::Initialize(int argc, const char * argv[])
{
    MPI::Init(argc, (char **&)argv);
    numprocs = MPI::COMM_WORLD.Get_size();
    myrank = MPI::COMM_WORLD.Get_rank();
    master = 0;
    loop_begin = myrank;
    loop_end = myrank;
    loop_offset = numprocs;
    return 0;
}

/********** Determine loop parameter **********/
/*! \fn int Determine_Loop(int n_file)
 *  \brief determine the begin/end/offset for loop */
int MPI_info::Determine_Loop(int n_file) {
    if (n_file < numprocs) {
        if (myrank > n_file - 1) {
            loop_end = -1;
        }
        loop_offset = 1;
    } else {
        // in order to let master processor become available
        // in fact, no special effect, just for future dev
        loop_begin = numprocs - 1 - myrank;
        loop_end = n_file - 1;
    }
    return 0;
}

/********** Barrier **********/
/*! \fn int Barrier()
 *  \brief wrapper of MPI Barrier */
int MPI_info::Barrier() {
    MPI::COMM_WORLD.Barrier();
    return 0;
}

/********** Finalization **********/
/*! \fn int Finalize()
 *  \brief wrapper of MPI Finalize() */
int MPI_info::Finalize() {
    MPI::Finalize();
    return 0;
}

/********** Destructor **********/
/*! \fn ~MPI_info()
 *  \brief destructor */
MPI_info::~MPI_info()
{
    ;
}

#endif // ENABLE_MPI
