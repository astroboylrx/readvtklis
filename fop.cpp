//
//  fop.cpp
//  readvtklis
//
//  Created by Rixin Li on 1/15/15.
//  Copyright (c) 2015 Rixin Li. All rights reserved.
//

#include "fop.h"



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
    cout << "USAGE: " << progname << " -c <n_cpu> -i <data_path> -b <data_basename> -s <post_name> -f <# (range(f1:f2))> -o <output_path_name> [--ParNum --RhoParMax --HeiPar --MeanSigma --VpecG --VertRho]\n" << endl;
    cout << "Example: ./readvtklis -c 16 -i comb -b Cout -s all -f 0:100 -o result.txt --ParNum  --RhoParMax --HeiPar --MeanSigma --VpecG --VertRho" << endl;
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
    MeanSigma_flag = 0;
    VpecG_flag = 0;
    VertRho_flag = 0;
    
    //Specifying the expected options
    static struct option long_options[] = {
        // These options set a flag
        {"ParNum", no_argument, &ParNum_flag, 1},
        {"RhoParMax", no_argument, &RhoParMax_flag, 1},
        {"HeiPar", no_argument, &HeiPar_flag, 1},
        {"MeanSigma", no_argument, &MeanSigma_flag, 1},
        {"VpecG", no_argument, &VpecG_flag, 1},
        {"VertRho", no_argument, &VertRho_flag, 1},
        // These options don't set a flag
        {"ncpu", required_argument, 0, 'c'},
        {"input", required_argument, 0, 'i'},
        {"basename", required_argument, 0, 'b'},
        {"postname", required_argument, 0, 's'},
        {"filenumber", required_argument, 0, 'f'},
        {"output", required_argument, 0, 'o'},
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
                    iof.data_path.assign(optarg);
                    //cout << "data_path is " << data_path << endl;
                    break;
                }
                case 'b': {
                    iof.data_basename.assign(optarg);
                    //cout << "data_basename is " << data_basename << endl;
                    break;
                }
                case 's': {
                    iof.post_name.assign(optarg);
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
                    iof.output_path_name.assign(optarg);
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
    paras.AllocateMemory(n_file);
    
    return 0;
}

/********** Generate file name in order **********/
/*! \fn int Generate_Filename()
 *  \brief generate filenames for processing */
int FileIO::Generate_Filename()
{
    if (*iof.data_path.rbegin() != '/') {
        iof.data_path.push_back('/');
    }
    for (int i = start_no; i <= end_no; i++) {
        stringstream ss;
        ss << setw(4) << setfill('0') << i;
        string file_no = ss.str();
        lis_filenames.push_back(iof.data_path+iof.data_basename+'.'+file_no+'.'+iof.post_name+".lis");
        vtk_filenames.push_back(iof.data_path+iof.data_basename+'.'+file_no+".vtk");
    }
    /*
     Print_Stars("Check Filenames");
     cout << "We generate " << lis_filenames.size() << " lis_filenames in total." << endl;
     cout << "The first one is " << *lis_filenames.begin() << endl;
     cout << "We generate " << vtk_filenames.size() << " vtk_filenames in total." << endl;
     cout << "The first one is " << *vtk_filenames.begin() << endl;
     */
    if (iof.output_path_name.length() == 0) {
        cout << "Error: No output path/name. " << endl;
        return 1;
    }
    if (MeanSigma_flag) {
        iof.output_sigma_path_name = iof.output_path_name.substr(0, iof.output_path_name.find_last_of('.'))+"_MeanSigma.txt";
    }
    if (VpecG_flag) {
        iof.output_vpecg_path_name = iof.output_path_name.substr(0, iof.output_path_name.find_last_of('.'))+"_VturbGas.txt";
    }
    if (VertRho_flag) {
        iof.output_vertrho_path_name = iof.output_path_name.substr(0, iof.output_path_name.find_last_of('.'))+"_VertRho.txt";
    }
    return 0;
}

/********** Check path and filename **********/
/*! \fn int Check_Input_Path_Filename()
 *  \brief check path and filename */
int FileIO::Check_Input_Path_Filename()
{
    Print_Stars("Check Input");
    cout << "data_path is " << iof.data_path << endl;
    cout << "data_basename is " << iof.data_basename << endl;
    cout << "post_name is " << iof.post_name << endl;
    cout << "start_no=" << start_no << ", end_no=" << end_no << endl;
    cout << "output_path_name is " << iof.output_path_name << endl;
    cout << "n_cpu is " << n_cpu << endl;
    Print_Stars("Check Filenames");
    cout << "We generate " << lis_filenames.size() << " lis_filenames in total." << endl;
    cout << "The first one is " << *lis_filenames.begin() << endl;
    cout << "We generate " << vtk_filenames.size() << " vtk_filenames in total." << endl;
    cout << "The first one is " << *vtk_filenames.begin() << endl;
    Print_Stars("Check Output");
    if (ParNum_flag || RhoParMax_flag || HeiPar_flag || MeanSigma_flag || VpecG_flag || VertRho_flag) {
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
        if (MeanSigma_flag) {
            cout << "<Sigma_g> and <Sigma_p>" << endl;
        }
        if (VpecG_flag) {
            cout << "Gas peculiar velocity components" << endl;
        }
        if (VertRho_flag) {
            cout << "Vertical Rho_g and Rho_p" << endl;
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
}

/********** Output data to file **********/
/*! \fn int Output_Data()
 *  \brief output data */
int FileIO::Output_Data()
{
    // shouldn't happen
    if (iof.output_path_name.length() == 0) {
        cout << "Error: No output path/name. " << endl;
        return 1;
    }
    
    // result.txt part
    if (ParNum_flag || RhoParMax_flag || HeiPar_flag) {
        ofstream file;
        file.open(iof.output_path_name.c_str(), ofstream::out);
        
        if (!file.is_open()) {
            cout << "Failed to open " << iof.output_path_name << endl;
            return 1;
        }
        
        file << setw(15) << setfill(' ') << "#orbit_time";
        if (ParNum_flag) {
            file << setw(15) << setfill(' ') << "N_par";
        }
        if (RhoParMax_flag) {
            file << setw(15) << setfill(' ') << "Max_Rhop";
        }
        if (HeiPar_flag) {
            file << setw(15) << setfill(' ') << "H_p";
            file << setw(15) << setfill(' ') << "Hp_in1sigma";
        }
        file << endl;
        for (int i = 0; i != n_file; i++) {
            file << setw(15) << scientific << paras.Otime[i];
            if (ParNum_flag) {
                file << setw(15) << paras.N_par[i];
            }
            if (RhoParMax_flag) {
                file << setw(15) << scientific << paras.Max_Rhop[i];
            }
            if (HeiPar_flag) {
                file << setw(15) << scientific << paras.Hp[i];
                file << setw(15) << scientific << paras.Hp_in1sigma[i];
            }
            file << endl;
        }
        file.close();
    }
    
    // <Sigma_p>_y part
    if (MeanSigma_flag) {
        ofstream file_MeanSigma;
        file_MeanSigma.open(iof.output_sigma_path_name.c_str(), ofstream::out);
        if (!file_MeanSigma.is_open()) {
            cout << "Failed to open " << iof.output_sigma_path_name << endl;
            return 1;
        }
        file_MeanSigma << "#The first row of data is orbit time. The first column is x (radial direction) coordinate. Others is data, but divided into two blocks, first is for gas, then for particles.";
        file_MeanSigma << endl;
        file_MeanSigma << setw(15) << setfill(' ') << 0.0;
        for (int i = 0; i != n_file; i++) {
            file_MeanSigma << setw(15) << scientific << paras.Otime[i];
        }
        file_MeanSigma << endl;
        for (int i = 0; i != 2*paras.dimensions[0]; i++) {
            file_MeanSigma << setw(15) << scientific << paras.ccx[i%paras.dimensions[0]];
            for (int j = 0; j != n_file; j++) {
                file_MeanSigma << setw(15) << paras.MeanSigma[j][i];
            }
            file_MeanSigma << endl;
        }
        file_MeanSigma.close();
    }
    
    if (VpecG_flag) {
        ofstream file_VpecG;
        file_VpecG.open(iof.output_vpecg_path_name.c_str(), ofstream::out);
        if (!file_VpecG.is_open()) {
            cout << "Failed to open " << iof.output_vpecg_path_name << endl;
            return 1;
        }
        file_VpecG << "#The first row is orbit time, The first column is z  (vertical direction) coordinate. Others is data, but divided into three blocks, first is Vx, then Vy and then Vz.";
        file_VpecG << endl;
        file_VpecG << setw(15) << setfill(' ') << 0.0;
        for (int i = 0; i != n_file; i++) {
            file_VpecG << setw(15) << scientific << paras.Otime[i];
        }
        file_VpecG << endl;
        for (int i = 0; i != 3*paras.dimensions[2]; i++) {
            file_VpecG << setw(15) << scientific << paras.ccz[i%paras.dimensions[2]];
            for (int j = 0; j != n_file; j++) {
                file_VpecG << setw(15) << scientific << paras.VpecG[j][i];
            }
            file_VpecG << endl;
        }
        file_VpecG.close();
    }
    
    if (VertRho_flag) {
        ofstream file_VertRho;
        file_VertRho.open(iof.output_vertrho_path_name.c_str(), ofstream::out);
        if (!file_VertRho.is_open()) {
            cout << "Failed to open " << iof.output_vertrho_path_name << endl;
            return 1;
        }
        file_VertRho << "#The first row of data is orbit time. The first column is z (vertical direction) coordinate. Others is data, but divided into two blocks, first is for gas, then for particles.";
        file_VertRho << endl;
        file_VertRho << setw(15) << setfill(' ') << 0.0;
        for (int i = 0; i != n_file; i++) {
            file_VertRho << setw(15) << scientific << paras.Otime[i];
        }
        file_VertRho << endl;
        for (int i = 0; i != 2*paras.dimensions[2]; i++) {
            file_VertRho << setw(15) << scientific << paras.ccz[i%paras.dimensions[2]];
            for (int j = 0; j != n_file; j++) {
                file_VertRho << setw(15) << scientific << paras.VertRho[j][i];
            }
            file_VertRho << endl;
        }
        file_VertRho.close();
    }

    return 0;
}
