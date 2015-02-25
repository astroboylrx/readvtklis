//
//  fop.h
//  readvtklis
//
//  Created by Rixin Li on 1/15/15.
//  Copyright (c) 2015 Rixin Li. All rights reserved.
//

#ifndef __readvtklis__fop__
#define __readvtklis__fop__

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

#define ENABLE_MPI

//#define RESERVE_PUSH_BACK
//#define FROM_ARRAY_TO_VECTOR
#define RESIZE_LIST

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
    
    // below is the varialbes used for MPI_Allreduce
    double *orbit_time;                             /*!< orbit time */
    double *max_rho_par;                            /*!< max par density */
    double *Hp;                                     /*!< par scale height */
    double *Hp_in1sigma;              /*!< par scale height deviation from Gaussian distribution */
    long *n_par;                                    /*!< number of par */
    long **cpuid_dist;                              /*!< cpuid distribution */
    double **Sigma_par_y;                           /*!< sigma_p averaged over y */
    double **Vturb_gas;                             /*!< V_turb_g averaged equatorially */
    
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
    
    /*! \fn int NewVars(int n_file, int n_cpu)
     *  \brief Allocate space for vars for MPI_Allreduce */
    int NewVars(int n_file, int n_cpu);
    
    /*! \fn ~MPI_info()
     *  \brief destructor */
    ~MPI_info();
    
};

extern MPI_info *myMPI; // declaration is in main function
#endif

/*! \class FileIO
 *  \brief Information about I/O
 */
class FileIO {
private:
    
public:
    string data_path;                               /*!< data path for read */
    string data_basename;                           /*!< the basename for data file */
    string post_name;                               /*!< the post name for data file */
    string output_path_name;                        /*!< the name for output file */
    string output_cpuid_path_name;                  /*!< the name for cpuid output file */
    string output_sigma_path_name;                  /*!< the name for output of <Sigma_par>_y */
    string output_vturbg_path_name;                 /*!< the name for output of <V_turb_gas>_xy */
    
    vector<string> lis_filenames;                   /*!< the vecotr for lis filenames */
    vector<string> vtk_filenames;                   /*!< the vector for vtk filenames */
    
    int start_no, end_no;                           /*!< the start_number/end_number for file */
    int n_file;                                     /*!< the number of file */
    int n_cpu;                                      /*!< the number of processors */
    
    int ParNum_flag,                                /*!< flag: total particle number */
        RhoParMax_flag,                             /*!< flag: maximum of particle density */
        HeiPar_flag,                                /*!< flag: particle scale height */
        CpuID_flag,                                 /*!< flag: particle's cpuid */
        SigmaParY_flag,                             /*!< flag: particle column density averaged azimuthally and vertically */
        VturbGas_flag,                              /*!< flag: gas turbulence velocity averaged equatorially, weighted by rho_g */
        //New_flag,                                 /*!< flag: example of new flag */
        UselessEnd_flag;                            /*!< flag: just in order to add flag conveniently */
    // after adding a flag:
    // MPI_info: add the var in MPI_info, add new statement in NewVars(), add delete in ~MPI_info()
    // FileIO: add var in option and new statement in Initialize(), add if statement in Check_Input_Path_Filename(), add delete in ~FileIO(), add output in Output_Data()
    //
    float mratio;                                   /*!< cst: total particle to gas mass ratio */
    float Omega_K;                                  /*!< cst: Omega_K = 1 */
    
    double *orbit_time;                             /*!< the orbital time */
    double *Hp;                                     /*!< the scale height of particles */
    double *Hp_in1sigma;              /*!< par scale height deviation from Gaussian distribution */
    double *max_rho_par;                            /*!< the max density of particles */
    long *n_par;                                    /*!< the number of particles */
    long **CpuID_dist;                              /*!< cpuid distribution */
    double **Sigma_par_y;                           /*!< <Sigma_par>_y */
    double **Vturb_gas;                             /*!< <V_turb_g>_xy */
    float ****V_gas_0;                             /*!< initial v_gas */
    int dimensions[3];                              /*!< the number of cells in each dimension */
    double spacing[3];                              /*!< the spacing of coordinate */
    double *ccx, *ccy, *ccz;                        /*!< cell center coordinates for plots */
    
    FileIO();                                       /*!< constructor */
    ~FileIO();                                      /*!< destructor */
    
    /*! \fn int Initialize(int argc, const char * argv[])
     *  \brief initialization */
    int Initialize(int argc, const char * argv[]);
    
    /*! \fn int Generate_Filename()
     *  \brief generate filenames for processing */
    int Generate_Filename();
    
    /*! \fn int Check_Input_Path_Filename()
     *  \brief check path and filename */
    int Check_Input_Path_Filename();
    
    /*! \fn int Print_Stars(string info)
     *  \brief print stars with info */
    int Print_Stars(string info);

    /*! \fn int Output_Data()
     *  \brief output data */
    int Output_Data();
    
    /*! \fn int Print_Usage(const char *progname)
     *  \brief print usage */
    int Print_Usage(const char *progname);
};
extern FileIO *fio;

#endif /* defined(__readvtklis__fop__) */
