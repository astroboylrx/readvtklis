//
//  readvtk.cpp
//  readvtklis
//
//  Created by Rixin Li on 1/14/15.
//  Copyright (c) 2015 Rixin Li. All rights reserved.
//

#include "readvtk.h"

/*************************************/
/*********CellData_Scaler*************/
/*************************************/

/********** Constructor **********/
CellData_Scaler::CellData_Scaler()
{
    ;
}

/********** Destructor **********/
CellData_Scaler::~CellData_Scaler()
{
    for (int i = 0; i != dimensions[0]; i++) {
        for (int j = 0; j != dimensions[1]; j++) {
            delete [] data[i][j];
        }
        delete [] data[i];
    }
    delete [] data;
}

/********** Read scaler data **********/
int CellData_Scaler::Read_Scaler_Data(string filename)
{
    return 0;
}

/*************************************/
/*********CellData_Vector*************/
/*************************************/

/********** Constructor **********/
CellData_Vector::CellData_Vector()
{
    ;
}

/********** Destructor **********/
CellData_Vector::~CellData_Vector()
{
    for (int i = 0; i != dimensions[0]; i++) {
        for (int j = 0; j != dimensions[1]; j++) {
            for (int k = 0; k != dimensions[2]; k++) {
                delete [] data[i][j][k];
            }
            delete [] data[i][j];
        }
        delete [] data[i];
    }
    delete [] data;
}

/********** Read scaler data **********/
int CellData_Vector::Read_Vector_Data(string filename)
{
    return 0;
}

/*************************************/
/************VtkFile******************/
/*************************************/

/********** Constructor **********/
VtkFile::VtkFile()
{
    ;
}

/********** Destructor **********/
VtkFile::~VtkFile()
{
    if (cd_scaler.size() > 0) {
        vector<CellData_Scaler> temp;
        cd_scaler.swap(temp);
    }
    if (cd_vector.size() > 0) {
        vector<CellData_Vector> temp;
        cd_vector.swap(temp);
    }
}

/********** Read header and record data position **********/
int VtkFile::Read_Header_Record_Pos(string filename)
{
    
    return 0;
}


