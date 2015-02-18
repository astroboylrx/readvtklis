//
//  fop.cpp
//  readvtklis
//
//  Created by Rixin Li on 1/15/15.
//  Copyright (c) 2015 Rixin Li. All rights reserved.
//

#include "fop.h"

#ifdef ENABLE_MPI
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

/********** New vars for MPI_Allreduce **********/
/*! \fn int NewVars(int n_file, int n_cpu)
 *  \brief Allocate space for vars for MPI_Allreduce */
int MPI_info::NewVars(int n_file, int n_cpu)
{
    orbit_time = new double[fio->n_file];
    max_rho_par = new double[fio->n_file];
    Hp = new double[fio->n_file];
    n_par = new long[fio->n_file];
    cpuid_dist = new long*[fio->n_file];
    Sigma_par_y = new double*[fio->n_file];
    for (int i = 0; i != fio->n_file; i++) {
        // initialize if you don't assign all of them values but use them for calculation
        orbit_time[i] = 0;
        max_rho_par[i] = 0;
        Hp[i] = 0;
        n_par[i] = 0;
        cpuid_dist[i] = new long[fio->n_cpu];
        for (int j = 0; j != fio->n_cpu; j++) {
            cpuid_dist[i][j] = 0;
        }
    }
    return 0;
}

/********** Destructor **********/
/*! \fn ~MPI_info()
 *  \brief destructor */
MPI_info::~MPI_info()
{
    delete [] orbit_time;
    delete [] max_rho_par;
    delete [] Hp;
    delete [] cpuid_dist;
    delete [] Sigma_par_y;
}

#endif // ENABLE_MPI


/********** Print stars contain info **********/
/*! \fn int Print_Stars(string info)
 *  \brief print stars with info */
int FileIO::Print_Stars(string info)
{   cout << endl;
    cout << setw(10) << setfill('*') << "*";
    cout << " " << info << " ";
    cout << setw(10) << setfill('*') << "*" << endl;
    //cout << endl;
    return 0;
}

/********** Print usage **********/
/*! \fn int Print_Usage(const char *progname)
 *  \brief print usage */
int FileIO::Print_Usage(const char *progname)
{
    cout << "USAGE: " << progname << " -c <n_cpu> -i <data_path> -b <data_basename> -s <post_name> -f <# (range(f1:f2))> -o <output_path_name> [--ParNum --RhoParMax --HeiPar --CpuID --SigmaParY]\n" << endl;
    cout << "Example: ./readvtklis -c 16 -i comb -b Cout -s all -f 0:100 -o result.txt --ParNum" << endl;
    return 0;
}

/********** Constructor **********/
/*! \fn FileIO::FileIO()
 *  \brief constructor */
FileIO::FileIO()
{
    ;
}

/********** Initialization **********/
/*! \fn int Initialize(int argc, const char * argv[])
 *  \brief initialization */
int FileIO::Initialize(int argc, const char * argv[])
{
    mratio = 0.02;
    int temp;
    // initialize flags
    ParNum_flag = 0;
    RhoParMax_flag = 0;
    HeiPar_flag = 0;
    CpuID_flag = 0;
    
    //Specifying the expected options
    static struct option long_options[] = {
        // These options set a flag
        {"ParNum", no_argument, &ParNum_flag, 1},
        {"RhoParMax", no_argument, &RhoParMax_flag, 1},
        {"HeiPar", no_argument, &HeiPar_flag, 1},
        {"CpuID", no_argument, &CpuID_flag, 1},
        {"SigmaParY", no_argument, &SigmaParY_flag, 1},
        // These options don't set a flag
        {"input", required_argument, 0, 'i'},
        {"basename", required_argument, 0, 'b'},
        {"postname", required_argument, 0, 's'},
        {"filenumber", required_argument, 0, 'f'},
        {"output", required_argument, 0, 'o'},
        {"ncpu", required_argument, 0, 'c'},
        // End
        {0,0,0,0}
    };
    
    if (argc < 12) {
#ifdef ENABLE_MPI
        if (myMPI->myrank == myMPI->master) {
#endif
            Print_Usage(argv[0]);
#ifdef ENABLE_MPI
        }
#endif
        exit(1);
    } else {
        /* for debug
#ifdef ENABLE_MPI
        if (myMPI->myrank == myMPI->master) {
#endif
            cout << "argc = " << argc << endl;
            for (int i = 0; i != argc; i++) {
                cout << "argv[" << i << "]: " << argv[i] << endl;
            }
#ifdef ENABLE_MPI
        }
#endif
         */

        while (1) {
            // getopt_long stores the option
            int option_index = 0;
            temp = getopt_long(argc, (char *const *)argv, "i:b:s:f:o:c:", long_options, &option_index);
            if (temp == -1) {
                break;
            }
            
            switch (temp) {
                case 0: {
                    // if this option set a flag, do nothing else now
                    if (long_options[option_index].flag != 0) {
                        break;
                    }
                    cout << "option " << long_options[option_index].name;
                    if (optarg) {
                        cout << " with arg " << optarg;
                    }
                    cout << endl;
                    break;
                }
                case 'i': {
                    data_path.assign(optarg);
                    //cout << "data_path is " << data_path << endl;
                    break;
                }
                case 'b': {
                    data_basename.assign(optarg);
                    //cout << "data_basename is " << data_basename << endl;
                    break;
                }
                case 's': {
                    post_name.assign(optarg);
                    //cout << "post_name is " << post_name << endl;
                    break;
                }
                case 'c': {
                    istringstream ifs;
                    ifs.str(optarg);
                    ifs >> n_cpu;
                    //cout << "numproc is " << n_cpu << endl;
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
                    break;
                }
                case 'o': {
                    output_path_name.assign(optarg);
                    //cout << "output_path_name is " << output_path_name << endl;
                    break;
                }
                case '?': {
#ifdef ENABLE_MPI
                    if (myMPI->myrank == myMPI->master) {
#endif
                        if (optopt == 'i' || optopt == 'b' || optopt == 's' || optopt == 'f' || optopt == 'o' || optopt == 'c')
                            fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                        else if (isprint (optopt))
                            fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                        else
                            fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
#ifdef ENABLE_MPI
                    }
#endif
                    exit(2);
                }
                default: {
#ifdef ENABLE_MPI
                    if (myMPI->myrank == myMPI->master) {
#endif
                        cout << temp << endl;
                        cout << "Argument wrong." << endl;
#ifdef ENABLE_MPI
                    }
#endif
                    abort();
                }
            }
        }
#ifdef ENABLE_MPI
        if (myMPI->myrank == myMPI->master) {
#endif
            // print any remaining command line arguments (not options)
            if (optind < argc) {
                cout << "Non-option ARGV-elements: ";
                while (optind < argc) {
                    cout << argv[optind++];
                }
                cout << endl;
            }
#ifdef ENABLE_MPI
        }
#endif
        
    }
    orbit_time = new double[n_file];
    n_par = new long[n_file];
    max_rho_par = new double[n_file];
    Hp = new double[n_file];
    CpuID_dist = new long*[n_file];
    Sigma_par_y = new double*[n_file];
    
    return 0;
}

/********** Generate file name in order **********/
/*! \fn int Generate_Filename()
 *  \brief generate filenames for processing */
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
     Print_Stars("Check Filenames");
     cout << "We generate " << lis_filenames.size() << " lis_filenames in total." << endl;
     cout << "The first one is " << *lis_filenames.begin() << endl;
     cout << "We generate " << vtk_filenames.size() << " vtk_filenames in total." << endl;
     cout << "The first one is " << *vtk_filenames.begin() << endl;
     */
    if (output_path_name.length() == 0) {
        cout << "Error: No output path/name. " << endl;
        return 1;
    }
    if (CpuID_flag) {
        output_cpuid_path_name = output_path_name.substr(0, output_path_name.find_last_of('.'))+"_CpuID.txt";
    }
    if (SigmaParY_flag) {
        output_sigma_path_name = output_path_name.substr(0, output_path_name.find_last_of('.'))+"_SigmaParY.txt";
    }
    return 0;
}

/********** Check path and filename **********/
/*! \fn int Check_Input_Path_Filename()
 *  \brief check path and filename */
int FileIO::Check_Input_Path_Filename()
{
    Print_Stars("Check Input");
    cout << "data_path is " << data_path << endl;
    cout << "data_basename is " << data_basename << endl;
    cout << "post_name is " << post_name << endl;
    cout << "start_no=" << start_no << ", end_no=" << end_no << endl;
    cout << "output_path_name is " << output_path_name << endl;
    cout << "n_cpu is " << n_cpu << endl;
    Print_Stars("Check Filenames");
    cout << "We generate " << lis_filenames.size() << " lis_filenames in total." << endl;
    cout << "The first one is " << *lis_filenames.begin() << endl;
    cout << "We generate " << vtk_filenames.size() << " vtk_filenames in total." << endl;
    cout << "The first one is " << *vtk_filenames.begin() << endl;
    Print_Stars("Check Output");
    if (ParNum_flag || RhoParMax_flag || HeiPar_flag || CpuID_flag || SigmaParY_flag) {
        cout << "Output includes: " << endl;
        if (ParNum_flag) {
            cout << "particle numbers" << endl;
        }
        if (RhoParMax_flag) {
            cout << "Max particle density" << endl;
        }
        if (HeiPar_flag) {
            cout << "Particle scale height" << endl;
        }
        if (CpuID_flag) {
            cout << "CPU ID" << endl;
        }
        if (SigmaParY_flag) {
            cout << "<Sigma_p>_y" << endl;
        }
    } else {
        cout << "Need output choices." << endl;
        exit(1);
    }
    return 0;
}


/********** Destructor **********/
/*! \fn FileIO::~FileIO()
 *  \brief destructor */
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
    delete [] CpuID_dist;
    
    /* memory problem, don't know why
    if (CpuID_flag) {
        for (int i = 0; i != n_file; i++) {
            delete [] CpuID_dist[i];
        }
        delete [] CpuID_dist;
    }
     */
    delete [] Sigma_par_y;
}

/********** Output data to file **********/
/*! \fn int Output_Data()
 *  \brief output data */
int FileIO::Output_Data()
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
    if (ParNum_flag) {
        file << setw(15) << setfill(' ') << "n_par";
    }
    if (RhoParMax_flag) {
        file << setw(15) << setfill(' ') << "max_rho_par";
    }
    if (HeiPar_flag) {
        file << setw(15) << setfill(' ') << "H_p";
    }
    file << endl;
    for (int i = 0; i != n_file; i++) {
        file << setw(15) << scientific << orbit_time[i];
        if (ParNum_flag) {
            file << setw(15) << n_par[i];
        }
        if (RhoParMax_flag) {
            file << setw(15) << scientific << max_rho_par[i];
        }
        if (HeiPar_flag) {
            file << setw(15) << scientific << Hp[i];
        }
        file << endl;
    }
    file.close();
    
    // cpuid part
    if (CpuID_flag) {
        ofstream file_CpuID;
        file_CpuID.open(output_cpuid_path_name.c_str(), ofstream::out);
        if (!file_CpuID.is_open()) {
            cout << "Failed to open " << output_cpuid_path_name << endl;
            return 1;
        }
        file_CpuID << setw(15) << setfill(' ') << "#orbit_time";
        for (int i = 0; i != n_cpu; i++) {
            std::stringstream temp_out;
            temp_out << i;
            file_CpuID << setw(15) << "CPU"+temp_out.str() << setfill(' ');
        }
        file_CpuID << endl;
        for (int i = 0; i != n_file; i++) {
            file_CpuID << setw(15) << scientific << orbit_time[i];
            for (int j = 0; j != n_cpu; j++) {
                file_CpuID << setw(15) << CpuID_dist[i][j];
            }
            file_CpuID << endl;
        }
        file_CpuID.close();
    }
    
    // <Sigma_p>_y part
    if (SigmaParY_flag) {
        ofstream file_SigmaParY;
        file_SigmaParY.open(output_sigma_path_name.c_str(), ofstream::out);
        if (!file_SigmaParY.is_open()) {
            cout << "Failed to open " << output_sigma_path_name << endl;
            return 1;
        }
        file_SigmaParY << setw(15) << setfill(' ') << "#The first row of data is orbit time. The first column is x (radial direction) coordinate. Others is data.";
        file_SigmaParY << endl;
        file_SigmaParY << setw(15) << setfill(' ') << 0.0;
        for (int i = 0; i != n_file; i++) {
            file_SigmaParY << setw(15) << scientific << orbit_time[i];
        }
        file_SigmaParY << endl;
        for (int i = 0; i != dimensions[0]; i++) {
            file_SigmaParY << setw(15) << scientific << ccx[i];
            for (int j = 0; j != n_file; j++) {
                file_SigmaParY << setw(15) << Sigma_par_y[j][i];
            }
            file_SigmaParY << endl;
        }
        file_SigmaParY.close();
    }

    return 0;
}
