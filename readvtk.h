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


/********** Cell data **********/
class VtkFile {
private:
    
public:
    int Version;
    int nx1, nx2, nx3, i, j, k;
    
};

class CellData_Scaler{
private:
    
public:
    int i;
};

class CellData_Vector{
private:
    
public:
    int i;
};



#endif /* defined(__readvtklis__readvtk__) */
