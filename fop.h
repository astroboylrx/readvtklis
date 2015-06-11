//
//  fop.h
//  readvtklis
//
//  Created by Rixin Li on 1/15/15.
//  Copyright (c) 2015 Rixin Li. All rights reserved.
//

#ifndef __readvtklis__fop__
#define __readvtklis__fop__

#include "global.h"

/*! \class IO_Filenames
 *  \brief contains all I/O filenames
 */
class IO_Filenames {
private:
    
public:
    string data_path;                               /*!< data path for read */
    string data_basename;                           /*!< the basename for data file */
    string post_name;                               /*!< the post name for data file */
    string output_path_name;                        /*!< the name for output file */
    string output_sigma_path_name;                  /*!< the name for output of MeanSigma */
    string output_vpecg_path_name;                  /*!< the name for output of VpecG */
    string output_vertrho_path_name;                /*!< the name for output of VertRho */
    string output_corrl_path_name;                  /*!< the name for output of Correlation Length */
#ifdef CorrValue
    string output_corrv_path_name;                  /*!< the name for output of Correlation Value */
#endif
    string data_level;                              /*!< the level of data */
    string data_domain;                             /*!< the domain of data */

};

/*! \class FileIO
 *  \brief Information about I/O
 */
class FileIO {
private:
    
public:
    IO_Filenames iof;                               /*!< contains all filenames for input/output */
    
    vector<string> lis_filenames;                   /*!< the vecotr for lis filenames */
    vector<string> vtk_filenames;                   /*!< the vector for vtk filenames */
    
    int start_no, end_no, interval;               /*!< the start_number/end_number/interval for file */
    int n_file;                                     /*!< the number of file */
    int n_cpu;                                      /*!< the number of processors */
    
    int ParNum_flag,                                /*!< flag: total particle number */
        RhoParMax_flag,                             /*!< flag: maximum of particle density */
        HeiPar_flag,                                /*!< flag: particle scale height */
        dSigma_flag,                                /*!< flag: change of gas surface density */

        MeanSigma_flag,                             /*!< flag: gas/particle column density averaged azimuthally and vertically */
        VpecG_flag,                                 /*!< flag: gas peculiar velocity components averaged equatorially, weighted by rho_g */
        VertRho_flag,                               /*!< flag: vertical structure of rho_g and rho_p */
        CorrL_flag,                                 /*!< flag: correlation length of Vgas */
        //New_flag,                                 /*!< flag: example of new flag */
        UselessEnd_flag;                            /*!< flag: just in order to add flag conveniently */
    // after adding a flag:
    // MPI_info: add the var in MPI_info, add new statement in NewVars(), add delete in ~MPI_info()
    // FileIO: add var in option and new statement in Initialize(), add if statement in Check_Input_Path_Filename(), add delete in ~FileIO(), add output in Output_Data()
    
    // some constants here, if needed, can be added into input parameters
    float mratio;                                   /*!< cst: total particle to gas mass ratio */
    float Omega_K;                                  /*!< cst: Omega_K = 1 */
    
    Paras2probe paras;                              /*!< contains parameters for the use of MPIAllreduce */
    
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
