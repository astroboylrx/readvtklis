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
class CellData_Scaler{
private:
    
public:
    string dataname, datatype, numcomp;
    // numcomp: number of components, default 1
    string tablename;
    long pos; // the position of data, offset from SEEK_SET
    int dimensions[3];
    float ***data;
    // in fact, now we are only dealing with float type
    
    // constructor and destructor
    CellData_Scaler();
    ~CellData_Scaler();
    
    
    // read scaler data
    int Read_Scaler_Data(string filename);
};

class CellData_Vector{
private:
    
public:
    string DataName, DataType;
    long pos; // the position of data, offset from SEEK_SET
    int dimensions[3];
    float ****data;
    // in fact, now we are only dealing with float type
    
    // constructor and destructor
    CellData_Vector();
    ~CellData_Vector();
    
    // read vector data
    int Read_Vector_Data(string filename);
    
};


class VtkFile {
private:
    
public:
    string Version;
    string Header;
    string FileFormat; // ACSII or BINARY
    string DatasetStructure;
    // in fact, now we are only dealing with STRUCTURED_POINTS
    // DATASET STRUCTURED_POINTS
    // DIMENSIONS n_x n_y n_z
    // ORIGIN x y z
    // SPACING s_x s_y s_z
    int dimensions[3];
    float origin[3], spacing[3];
    long CellData; // number of CELL_DATA, should be equal to the product of dimensions
    long PointData;
    // in fact, now we are only dealing with CELL_DATA
    vector<CellData_Scaler> cd_scaler;
    vector<CellData_Vector> cd_vector;
    
    // constructor and destructor
    VtkFile();
    ~VtkFile();
    
    // read header and record data position
    int Read_Header_Record_Pos(string filename);
    
};



#endif /* defined(__readvtklis__readvtk__) */
