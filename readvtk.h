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

class CellData_Scaler{
private:
    
public:
    string dataname, datatype, numcomp; // number of components, default 1
    string tablename;
    vector<float> data;
};

class CellData_Vector{
private:
    
public:
    string DataName, DataType;
    
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
    
    
    
    
    
};



#endif /* defined(__readvtklis__readvtk__) */
