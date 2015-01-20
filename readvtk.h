//
//  readvtk.h
//  readvtklis
//
//  Created by Rixin Li on 1/14/15.
//  Copyright (c) 2015 Rixin Li. All rights reserved.
//

#ifndef __readvtklis__readvtk__
#define __readvtklis__readvtk__

#include "fop.h"

/********** vtk file format **********/
// due to the need for high performance, we use pointer based multi-dimensional arrays here
class CellData_Scalar{
private:
    
public:
    string dataname, datatype;
    int numcomp;
    // numcomp: number of components, default 1
    string tablename;
    long pos; // the position of data, offset from SEEK_SET
    int dimensions[3];
    float ***data;
    // in fact, now we are only dealing with float type
    
    // constructor and destructor
    CellData_Scalar(int *dimensions_been_told);
    ~CellData_Scalar();
    
    // construct data
    int Initialize_Data(int *dimensions_been_told);
    /* types should be read from right to left. For example, this is (starting from the rightmost *) a pointer to a constant pointer to an int.
     int * const *x;
     so here, (int * & dimension) is correct order, otherwise, (int & * dimension) would be a pointer to a reference, which is not possible.
     ...passing a pointer is easier, passing a reference to a pointer comes across error while binding a tempoaray to the reference, which I can't understand for now
     */
    
    // read scalar data
    int Read_Scalar_Data(string filename);
    
    // free data memory
    int Free_Data();
};

class CellData_Vector{
private:
    
public:
    string dataname, datatype;
    long pos; // the position of data, offset from SEEK_SET
    int dimensions[3];
    float ****data;
    // in fact, now we are only dealing with float type
    
    // constructor and destructor
    CellData_Vector(int *dimensions_been_told);
    ~CellData_Vector();
    
    // construct data
    int Initialize_Data(int *dimensions_been_told);
    
    // read vector data
    int Read_Vector_Data(string filename);
    
    // free data memory
    int Free_Data();
    
};


class VtkFile {
private:
    
public:
    string Version;
    string Header;
    double time;
    string FileFormat; // ACSII or BINARY
    string DatasetStructure;
    // in fact, now we are only dealing with STRUCTURED_POINTS
    // DATASET STRUCTURED_POINTS
    // DIMENSIONS n_x n_y n_z
    // ORIGIN x y z
    // SPACING s_x s_y s_z
    int dimensions[3];
    double origin[3], spacing[3], ****cell_center; // cell center coor
    long n_CellData; // number of CELL_DATA, should be equal to the product of dimensions
    long n_PointData;
    // in fact, now we are only dealing with CELL_DATA
    vector<CellData_Scalar> cd_scalar;
    vector<CellData_Vector> cd_vector;
    // relative total gas mass and particle mass
    // and maximum density of gas and par(rescaled)
    double m_gas, m_par, max_rho_gas, max_rho_par;
    
    // constructor and destructor
    VtkFile();
    ~VtkFile();
    
    // construct coordinate grid
    int Construct_Coor();
    
    // read header and record data position
    int Read_Header_Record_Pos(string filename);
    
    // read data
    int Read_Data(string filename);
    
    // print file info
    int Print_File_Info();
    
    // calculate mass and find maximum
    int Calculate_Mass_Find_Max();
};



#endif /* defined(__readvtklis__readvtk__) */
