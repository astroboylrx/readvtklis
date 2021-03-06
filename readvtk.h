//
//  readvtk.h
//  readvtklis
//
//  Created by Rixin Li on 1/14/15.
//  Copyright (c) 2015 Rixin Li. All rights reserved.
//

#ifndef __readvtklis__readvtk__
#define __readvtklis__readvtk__

#include "global.h"
#include "fop.h"

// Due to the need for high performance, we use pointer based multi-dimensional arrays here.

/*! \class CellData_Scalar
 *  \brief Contains one kind of scalar cell data in vtk file
 */
class CellData_Scalar{
private:
    
public:
    string dataname;                            /*!< data name */
    string datatype;                            /*!< data type */
    int numcomp;                                /*!< number of components, default 1 */
    string tablename;                           /*!< table name, default: "default" */
    long pos;                                   /*!< the position of data, offset from SEEK_SET */
    int dimensions[3];                          /*!< the number of cells in each dimension */
    float ***data;                              /*!< data array, now we are only dealing with float type */
    
    CellData_Scalar(int *dimensions_been_told); /*!< constructor */
    ~CellData_Scalar();                         /*!< destructor */
    
    /*! \fn int Initialize_Data(int *dimensions_been_told)
     *  \brief initialize data */
    int Initialize_Data(int *dimensions_been_told);
    
    /* types should be read from right to left. For example, this is (starting from the rightmost *) a pointer to a constant pointer to an int.
     int * const *x;
     so, (int * & dimension) is correct order, otherwise, (int & * dimension) would be a pointer to a reference, which is not possible.
     ...passing a pointer is easier, passing a reference to a pointer comes across error while binding a tempoaray to the reference, which I can't understand for now
     */
    
    /*! \fn int Read_Scalar_Data(string filename)
     *  \brief read scalar data */
    int Read_Scalar_Data(string filename);
    
    /*! \fn int Free_Data()
     *  \brief free the data memory */
    int Free_Data();
};

/*! \class CellData_Vector
 *  \brief Contains one kind of vector cell data in vtk file
 */
class CellData_Vector{
private:
    
public:
    string dataname;                            /*!< data name */
    string datatype;                            /*!< data type */
    long pos;                                   /*!< the position of data, offset from SEEK_SET */
    int dimensions[3];                          /*!< the number of cells in each dimension */
    float ****data;                             /*!< data array, now we are only dealing with float type */
    
    CellData_Vector(int *dimensions_been_told); /*!< constructor */
    ~CellData_Vector();                         /*!< destructor */
    
    /*! \fn int Initialize_Data(int *dimensions_been_told)
     *  \brief initialize data */
    int Initialize_Data(int *dimensions_been_told);
    
    /*! \fn int Read_Vector_Data(string filename)
     *  \brief read vector data */
    int Read_Vector_Data(string filename);
    
    /*! \fn int Free_Data()
     *  \brief free the data memory */
    int Free_Data();
    
};

/*! \class VtkFile
 *  \brief Information about the entire vtk file
 */
class VtkFile {
private:
    
public:
    string Version;                             /*!< vtk version */
    string Header;                              /*!< vtk header */
    float time;                                 /*!< current time */
    string FileFormat;                          /*!< file format, ACSII or BINARY */
    string DatasetStructure;                    /*!< dataset structure, now we are only dealing with STRUCTURED_POINTS */
    // DATASET STRUCTURED_POINTS
    // DIMENSIONS n_x n_y n_z
    // ORIGIN x y z
    // SPACING s_x s_y s_z
    int dimensions[3];                          /*!< the number of cells in each dimension */
    int kps, kpe;                               /*!< the begin/end indices of the cells that correspond to the lower/upper initial 4*Hp=0.1 */
    float origin[3];                            /*!< the coordinate of origin point */
    float spacing[3];                           /*!< the spacing of coordinate */
    float L[3];                                 /*!< the length of box */
    float cell_volume;                          /*!< the cell volume */
    float ****cell_center;                      /*!< the coordinate cell center */
    long n_CellData;                            /*!< number of CELL_DATA, should be equal to the product of dimensions */
    // in fact, now we are only dealing with CELL_DATA
    vector<CellData_Scalar> cd_scalar;          /*!< vector of cell data */
    vector<CellData_Vector> cd_vector;          /*!< vector of vector data */
    double Sigma_gas_0;                         /*!< initial column density of gas */
    double Sigma_gas_0_inbox;                   /*!< initial column density of gas truncated by vertical box size */
    double Sigma_gas_0_mid;                     /*!< initial gas density in the midplane's 2 cells */
    double Sigma_gas_0_in2Hp;                   /*!< initial gas density inside the particle scale height */
    float m_gas;                                /*!< total gas mass */
    float m_par;                                /*!< total particle mass */
    float Max_Rhog;                             /*!< maximum density of gas */
    float Max_Rhop;                             /*!< maximum density of particle */
    double RpAV;                                /*!< <rho_p> */
    double RpSQ;                                /*!< <rho_p^2>^0.5 */
    double RpQU;                                /*!< <rho_p^4>^0.25 */
    double dSigma;                              /*!< change of gas surface density due to outflow */
    double dynscal[32];                         /*!< volume-averaged gas properties, in the order of p_gas[x,y,z,tot]/V, Ek_gas[x,y,z,tot]/V, then area-averaged gas properties, in the order of p_gas[x,y,z,tot]/A, Ek_gas[x,y,z,tot]/A, then volume-averaged particle properties, in the order of p_par[x,y,z,tot]/V, Ek_par[x,y,z,tot]/V, then area-averaged particle properties, in the order of p_par[x,y,z,tot]/A, Ek_par[x,y,z,tot]/A */
    double dynscal2[16];                        /*!< volume-averaged gas properties inside z = +/- 0.1 (the overlap region with the smallest box), in the order of p_gas[x,y,z,tot]/V, Ek_gas[x,y,z,tot]/V, then area-averaged gas properties, in the order of p_gas[x,y,z,tot]/A, Ek_gas[x,y,z,tot]/A */
    double *GPME;                               /*!< see the definition in global.h */
    double *GPMEPar;                            /*!< see the definition in global.h */
    
    VtkFile();                                  /*!< constructor */
    ~VtkFile();                                 /*!< destructor */
    
    /*! \fn int Construct_Coor(int *dimensions_been_told)
     *  \brief consturct coordinate grid */
    int Construct_Coor(int *dimensions_been_told);

    /*! \fn int Read_Header_Record_Pos(string filename)
     *  \brief read header and record data position */
    int Read_Header_Record_Pos(string filename);
    
    /*! \fn int Read_Data(string filename)
     *  \brief read data */
    int Read_Data(string filename);
    
    /*! \fn int Print_File_Info()
     *  \brief print file info */
    int Print_File_Info();
    
    /*! \fn int Calculate_Mass_Find_Max()
     *  \brief calculate mass and find maximum and mass loss */
    int Calculate_Mass_Find_Max();
    
    /*! \fn int VpecG(double *VpecG)
     *  \brief return gas peculiar velocity components averaged horizontally at each z, weighted by rho_g 
     */
    int VpecG(double *VpecG);
    
    /*! \fn int MeanSigma(double *MeanSigma)
     *  \brief calculate sigma_g and sigma_p averaged over y */
    int MeanSigma(double *MeanSigma);
    
    /*! \fn int Sigma(double **Sigma)
     *  \brief calculate the real sigma_g and sigma_p */
    int Sigma(double **Sigma);
    
    /*! \fn int VertRho(double *VertRho)
     *  \brief calculate vertical rho_g and rho_p */
    int VertRho(double *VertRho);
    
#ifdef ENABLE_FFT
    double *fft1d_data;
    fftw_complex *spectra;
    fftw_plan plan;
    
    /*! \fn int PowerSpectra()
     *  \brief apply fft to data */
    int PowerSpectra1DX(CellData_Scalar const &s_data, double *ps);
    int PowerSpectra1DX(CellData_Vector const &v_data, int const component, double *ps);
#endif // ENABLE_FFT
    
    /*! \fn int CorrLen(float *CoorL)
     *  \brief calculate the correlation length */
    int CorrLen(float *CorrL
#ifdef CorrValue
                , float *CorrV
#endif
                );
    
    /*! \fn int GasPar()
     *  \brief calculate the basic dynamic info */
    int GasPar();
    
};

#endif /* defined(__readvtklis__readvtk__) */
