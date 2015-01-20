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
#include <ctime>
#include <cassert>
#include <algorithm>
#include "mpi.h"

using namespace::std;

#define ENABLE_MPI

// After tests, reserve and push_back is fastest
//#define RESERVE_PUSH_BACK
//#define FROM_ARRAY_TO_VECTOR
#define RESIZE_LIST

#ifdef ENABLE_MPI
class MPI_info {
private:
    
public:
    int numprocs, myrank, master;
    int loop_begin, loop_end, loop_offset;
    
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
    int Barrier() {
        MPI::COMM_WORLD.Barrier();
        return 0;
    }
    int Finalize() {
        MPI::Finalize();
        return 0;
    }
};

extern MPI_info *myMPI; // declaration is in main function
#endif

class FileIO {
private:
    
public:
    // all sorts of path and name
    string data_path;
    string data_basename;
    string post_name;
    string output_path_name;
    
    // lis and vtk filenames
    vector<string> lis_filenames;
    vector<string> vtk_filenames;
    
    // the start number and the end number to process
    int start_no, end_no, n_file;
    
    // time, H_p, max_mp
    double *orbit_time, *Hp, *max_rho_par;
    // total particle to gas mass ratio
    float mratio;
    
    // constructor and destructor
    FileIO(int argc, const char * argv[]);
    ~FileIO();
    
    // generate file name in order
    int Generate_Filename();
    
    // check path and filename
    int Check_Path_Filename();
    
    // print stars contain info
    int print_stars(string info);
    
    // output data to file
    int output_data();

};


#endif /* defined(__readvtklis__fop__) */
