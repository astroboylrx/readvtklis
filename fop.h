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

using namespace::std;

// After tests, reserve and push_back is fastest
//#define RESERVE_PUSH_BACK
//#define FROM_ARRAY_TO_VECTOR
#define RESIZE_LIST

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
    int start_no, end_no;
    
    // constructor and destructor
    FileIO(int argc, const char * argv[]);
    ~FileIO();
    
    // generate file name in order
    int Generate_Filename();
    
    // print stars contain info
    int print_stars(string info);

};


#endif /* defined(__readvtklis__fop__) */
