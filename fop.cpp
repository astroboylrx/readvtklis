//
//  fop.cpp
//  readvtklis
//
//  Created by Rixin Li on 1/15/15.
//  Copyright (c) 2015 Rixin Li. All rights reserved.
//

#include "fop.h"

/********** Print stars contain info **********/
int FileIO::print_stars(string info)
{   cout << endl;
    cout << setw(10) << setfill('*') << "*";
    cout << " " << info << " ";
    cout << setw(10) << setfill('*') << "*" << endl;
    //cout << endl;
    return 0;
}

/********** Constructor **********/
FileIO::FileIO(int argc, const char * argv[])
{
    mratio = 0.02;
    int temp;
    int iflag = 0, bflag = 0, sflag = 0, fflag = 0, oflag = 0;
    if (argc < 11) {
        cout << "USAGE: " << argv[0] << " -i <data_path> -b <data_basename> -s <post_name> -f <# (range(f1:f2))> -o <output_path_name>\n" << endl;
        cout << "Example: ./readvtklis -i comb -b Cout -s all -f 0:100 -o result.txt" << endl;
        exit(1);
    } else {
        //print_stars("Check Path");
        while ((temp = getopt(argc, (char **)argv, "i:b:s:f:o:")) != -1) {
            switch (temp) {
                case 'i': {
                    data_path.assign(optarg);
                    //cout << "data_path is " << data_path << endl;
                    iflag = 1;
                    break;
                }
                case 'b': {
                    data_basename.assign(optarg);
                    //cout << "data_basename is " << data_basename << endl;
                    bflag = 1;
                    break;
                }
                case 's': {
                    post_name.assign(optarg);
                    //cout << "post_name is " << post_name << endl;
                    sflag = 1;
                    break;
                }
                case 'f': {
                    //sscanf(optarg,"%d:%d", &start_no, &end_no);
                    // obviously, the line above can replace the four line below
                    istringstream ifs;
                    ifs.str(optarg);
                    char tempchar;
                    ifs >> start_no >> tempchar >> end_no;
                    
                    if (start_no < 0) {
                        cout << "The start number should be larger than 0. (Auto fix to 0)" << endl;
                        start_no = 0;
                    }
                    if (end_no < start_no) {
                        cout << "The end number should be larger than the start number. (Auto fix to start number + 1)." << endl;
                        end_no += start_no;
                    }
                    //cout << "start_no=" << start_no << ", end_no=" << end_no << endl;
                    n_file = end_no - start_no + 1;
                    fflag = 1;
                    break;
                }
                case 'o': {
                    output_path_name.assign(optarg);
                    //cout << "output_path_name is " << output_path_name << endl;
                    oflag = 1;
                    break;
                }
                case '?': {
                    if (optopt == 'i' || optopt == 'b' || optopt == 's' || optopt == 'f' || optopt == 'o')
                        fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                    else if (isprint (optopt))
                        fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                    else
                        fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
                    exit(2);
                }
                default: {
                    cout << "Argument wrong." << endl;
                    abort();
                }
            }
        }
        // Since argc needs to > 11, so the five IF statements should not be called at all.
        if (iflag == 0) {
            cout << "We need data_path. Abort." << endl;
            abort();
        }
        if (bflag == 0) {
            cout << "We need data_basename. Abort." << endl;
            abort();
        }
        if (sflag == 0) {
            cout << "We need post_name. Abort." << endl;
            abort();
        }
        if (fflag == 0) {
            cout << "We need file range. Abort." << endl;
            abort();
        }
        if (oflag == 0) {
            cout << "We need output_path_name. Abort." << endl;
            abort();
        }
    }
    orbit_time = new double[end_no-start_no+1];
    Hp = new double[end_no-start_no+1];
    max_rho_par = new double[end_no-start_no+1];
}

/********** Generate file name in order **********/
int FileIO::Generate_Filename()
{
    if (*data_path.rbegin() != '/') {
        data_path.push_back('/');
    }
    for (int i = start_no; i <= end_no; i++) {
        stringstream ss;
        ss << setw(4) << setfill('0') << i;
        string file_no = ss.str();
        lis_filenames.push_back(data_path+data_basename+'.'+file_no+'.'+post_name+".lis");
        vtk_filenames.push_back(data_path+data_basename+'.'+file_no+".vtk");
    }
    /*
    print_stars("Check Filenames");
    cout << "We generate " << lis_filenames.size() << " lis_filenames in total." << endl;
    cout << "The first one is " << *lis_filenames.begin() << endl;
    cout << "We generate " << vtk_filenames.size() << " vtk_filenames in total." << endl;
    cout << "The first one is " << *vtk_filenames.begin() << endl;
     */
    
    return 0;
}

/********** Check path and filename **********/
int FileIO::Check_Path_Filename()
{
    print_stars("Check Path");
    cout << "data_path is " << data_path << endl;
    cout << "data_basename is " << data_basename << endl;
    cout << "post_name is " << post_name << endl;
    cout << "start_no=" << start_no << ", end_no=" << end_no << endl;
    cout << "output_path_name is " << output_path_name << endl;
    print_stars("Check Filenames");
    cout << "We generate " << lis_filenames.size() << " lis_filenames in total." << endl;
    cout << "The first one is " << *lis_filenames.begin() << endl;
    cout << "We generate " << vtk_filenames.size() << " vtk_filenames in total." << endl;
    cout << "The first one is " << *vtk_filenames.begin() << endl;
    return 0;
}


/********** Destructor **********/
FileIO::~FileIO()
{
    vector<string> temp1, temp2;
    lis_filenames.swap(temp1);
    vtk_filenames.swap(temp2);
    
    /***
     By clear() method, a reallocation is not guaranteed to happen, and the vector capacity is not guaranteed to change due to calling this function. A typical alternative that forces a reallocation is to use swap:
     vector<T>().swap(x);   // clear x reallocating
     In Mac OS, using Xcode, test verifies that clear() won't guarantee the free of memory.
     */
    delete [] orbit_time;
    delete [] Hp;
    delete [] max_rho_par;
    
}

/********** Output data to file **********/
int FileIO::output_data()
{
    // shouldn't happen
    if (output_path_name.length() == 0) {
        cout << "Error: No output path/name. " << endl;
        return 1;
    }
    ofstream file;
    file.open(output_path_name.c_str(), ofstream::out);
    if (!file.is_open()) {
        cout << "Failed to open " << output_path_name << endl;
        return 1;
    }
    file << setw(15) << setfill(' ') << "#orbit_time";
    file << setw(15) << setfill(' ') << "max_rho_par";
    file << setw(15) << setfill(' ') << "H_p";
    file << endl;
    for (int i = 0; i != end_no-start_no+1; i++) {
        file << setw(15) << scientific << orbit_time[i];
        file << setw(15) << scientific << max_rho_par[i];
        file << setw(15) << scientific << Hp[i];
        file << endl;
    }
    file.close();
    
    return 0;
}
