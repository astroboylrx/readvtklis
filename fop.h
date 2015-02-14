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
    
    /*! \fn int Initialize(int argc, const char * argv[])
     *  \brief MPI initializaion */
    int Initialize(int argc, const char * argv[]) {
        MPI::Init(argc, (char **&)argv);
        numprocs = MPI::COMM_WORLD.Get_size();
        myrank = MPI::COMM_WORLD.Get_rank();
        master = 0;
        loop_begin = myrank;
        loop_end = myrank;
        loop_offset = numprocs;
        return 0;
    }
    
    /*! \fn int Determine_Loop(int n_file)
     *  \brief determine the begin/end/offset for loop */
    int Determine_Loop(int n_file) {
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
    
    /*! \fn int Barrier()
     *  \brief wrapper of MPI Barrier */
    int Barrier() {
        MPI::COMM_WORLD.Barrier();
        return 0;
    }
    
    /*! \fn int Finalize()
     *  \brief wrapper of MPI Finalize() */
    int Finalize() {
        MPI::Finalize();
        return 0;
    }
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
    
    vector<string> lis_filenames;                   /*!< the vecotr for lis filenames */
    vector<string> vtk_filenames;                   /*!< the vector for vtk filenames */
    
    int start_no, end_no;                           /*!< the start_number/end_number for file */
    int n_file;                                     /*!< the number of file */
    int n_cpu;                                      /*!< the number of processors */
    
    int ParNum_flag,                                /*!< flag: total particle number */
        RhoParMax_flag,                             /*!< flag: maximum of particle density */
        HeiPar_flag,                                /*!< flag: particle scale height */
        CpuID_flag,                                 /*!< flag: particle's cpuid */
        //New_flag,                                 /*!< flag: example of new flag */
        UselessEnd_flag;                            /*!< flag: just in order to add flag conveniently */
    // time, H_p, max_mp, n_par
    double *orbit_time;                             /*!< the orbital time */
    double *Hp;                                     /*!< the scale height of particles */
    double *max_rho_par;                            /*!< the max density of particles */
    long *n_par;                                    /*!< the number of particles */
    float mratio;                                   /*!< total particle to gas mass ratio */
    long **CpuID_dist;                                /*!< cpuid distribution */
    
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
