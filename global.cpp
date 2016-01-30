//
//  global.cpp
//  readvtklis
//
//  Created by Rixin Li on 2/26/15.
//  Copyright (c) 2015 Rixin Li. All rights reserved.
//

#include "global.h"

/*************************************/
/************Paras2probe**************/
/*************************************/

/*! \fn int Initialize()
 *  \brief allocate memory to parameters (first level)
 */
int Paras2probe::AllocateMemory(int n_file)
{
    Otime = new float[n_file];
    N_par = new long[n_file];
    Max_Rhop = new float[n_file];
    RpAV = new double[n_file];
    RpSQ = new double[n_file];
    RpQU = new double[n_file];
    Hp = new double[n_file];
    Hp_in1sigma = new float[n_file];
    dSigma = new double[n_file];
    for (int i = 0; i != n_file; i++) {
        // initialize if you don't assign all of them values but use them for calculation eventually
        Otime[i] = 0;
        N_par[i] = 0;
        Max_Rhop[i] = 0;
        RpAV[i] = 0;
        RpSQ[i] = 0;
        RpQU[i] = 0;
        Hp[i] = 0;
        Hp_in1sigma[i] = 0;
        dSigma[i] = 0;
    }
    
    MeanSigma = new double*[n_file];
    VpecG = new double*[n_file];
    VertRho = new double*[n_file];
    CorrL = new float*[n_file];
#ifdef CorrValue
    CorrV = new float*[n_file];
#endif
    GasHst = new double*[n_file];
    GasHst2 = new double*[n_file];
    ParHst = new double*[n_file];
    ParLis = new double*[n_file];
    GPME = new double*[n_file];
    GPMEPar = new double*[n_file];
    for (int i = 0; i != n_file; i++) {
        GasHst[i] = new double[16];
        GasHst2[i] = new double[16];
        ParHst[i] = new double[16];
        ParLis[i] = new double[20];
        for (int j = 0; j != 16; j++) {
            GasHst[i][j] = 0;
            GasHst2[i][j] = 0;
            ParHst[i][j] = 0;
            ParLis[i][j] = 0;
        }
        for (int j = 16; j != 20; j++) {
            ParLis[i][j] = 0;
        }
    }
    
    return 0;
}

/*! \fn int AllocateSubMemory(int n_file, int *dimensions)
 *  \brief allocate memory to parameters (second level)
 */
int Paras2probe::AllocateSubMemory(int n_file, int *dimensions)
{
    Sigma = new double*[4*dimensions[0]];
    for (int i = 0; i != 4*dimensions[0]; i++) {
        Sigma[i] = new double[dimensions[1]];
        for (int j = 0; j != dimensions[1]; j++) {
            Sigma[i][j] = 0;
        }
    }
    for (int i = 0; i != n_file; i++) {
        MeanSigma[i] = new double[4*dimensions[0]];
        for (int j = 0; j != 4*dimensions[0]; j++) {
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
        CorrL[i] = new float[3*dimensions[2]];
        for (int j = 0; j != 3*dimensions[2]; j++) {
            CorrL[i][j] = 0;
        }
#ifdef CorrValue
        int Nz = dimensions[2], Nlines = Nz * (dimensions[0]/2+1);
        CorrV[i] = new float[3*Nlines];
        for (int j = 0; j != 3*Nlines; j++) {
            CorrV[i][j] = 0;
        }
#endif
        GPME[i] = new double[4*(dimensions[2]+1)*9];
        for (int j = 0; j != 4*(dimensions[2]+1)*9; j++) {
            GPME[i][j] = 0;
        }
        GPMEPar[i] = new double[4*(dimensions[2]+1)*9];
        for (int j = 0; j != 4*(dimensions[2]+1)*9; j++) {
            GPMEPar[i][j] = 0;
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
    delete [] RpAV;
    delete [] RpSQ;
    delete [] RpQU;
    delete [] Hp;
    delete [] Hp_in1sigma;
    delete [] dSigma;
    delete [] MeanSigma;
    delete [] Sigma;
    delete [] VpecG;
    delete [] VertRho;
    delete [] CorrL;
#ifdef CorrValue
    delete [] CorrV;
#endif
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
    wait_time = 0;
    return 0;
}

/********** Determine loop parameter **********/
/*! \fn int Determine_Loop(int n_file)
 *  \brief determine the begin/end/offset for loop */
int MPI_info::Determine_Loop(int n_file)
{
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
int MPI_info::Barrier()
{
    MPI::COMM_WORLD.Barrier();
    return 0;
}

/********** T **********/
/*! \fn double T()
 *  \brief wrapper of MPI WTime */
double MPI_info::T()
{
    return MPI::Wtime();
}

/********** prank **********/
/*! \fn string prank()
 *  \brief return string of Processor myrank: */
string MPI_info::prank()
{
    ostringstream oss;
    oss << "Processor " << myrank << ": ";
    return oss.str();
}

/********** Finalization **********/
/*! \fn int Finalize()
 *  \brief wrapper of MPI Finalize() */
int MPI_info::Finalize()
{
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
