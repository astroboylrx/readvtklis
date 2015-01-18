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
    /*
    for (int i = 0; i != dimensions[0]; i++) {
        for (int j = 0; j != dimensions[1]; j++) {
            delete [] data[i][j];
        }
        delete [] data[i];
    }
    delete [] data;
     */
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
    string tempstring;
    ifstream file (filename.c_str(), ios::in);
    if (!file.is_open()) {
        cout << "Failed to open " << filename << endl;
        return 1;
    }
    
    getline(file, Version);
    if (Version.compare("# vtk DataFile Version 3.0") != 0 && Version.compare("# vtk DataFile Version 2.0") != 0) {
        cout << "First line of " << filename << " is " << Version << endl;
        return 1;
    }
    //cout << "Version: " << Version << endl;
    
    getline(file, Header);
    if (Header.find("CONSERVED") != string::npos) {
        size_t time_pos = Header.find("time= ");
        // stod() reads to the end of the number, so we can ignore those afterwards
        time = stod(Header.substr(time_pos+6));
        
    }
    //cout << setprecision(6) << scientific << "time: " << time << endl;
    
    getline(file, FileFormat);
    if (FileFormat.compare("BINARY") != 0) {
        cout << "Unsopported file format: " << FileFormat << endl;
        return 1;
    }
    //cout << FileFormat << endl;
    
    getline(file, DatasetStructure);
    if (DatasetStructure.compare("DATASET STRUCTURED_POINTS") != 0) {
        cout << "Unsopported dataset structure: " << DatasetStructure << endl;
        return 1;
    }
    //cout << DatasetStructure << endl;
    
    getline(file, tempstring, ' ');
    if (tempstring.compare("DIMENSIONS") == 0) {
        istringstream iss;
        getline(file, tempstring);
        iss.str(tempstring);
        iss >> dimensions[0] >> dimensions[1] >> dimensions[2];
        //We want to store the number of grid cells, not the number of grid cell corners
        for (int i = 0; i != 3; i++) {
            dimensions[i]--;
        }
        /* alternative:
        getline(file, tempstring, ' ');
        dimensions[0] = stoi(tempstring);
        getline(file, tempstring, ' ');
        dimensions[1] = stoi(tempstring);
        getline(file, tempstring);
        dimensions[2] = stoi(tempstring);
         */
    } else {
        cout << "No dimensions info: " << endl;
        return 1;
    }
    //cout << dimensions[0] << " " << dimensions[1] << " " << dimensions[2] << endl;
    
    getline(file, tempstring, ' ');
    if (tempstring.compare("ORIGIN") == 0) {
        istringstream iss;
        getline(file, tempstring);
        iss.str(tempstring);
        iss >> origin[0] >> origin[1] >> origin[2];
        
    } else {
        cout << "No origin info: " << endl;
        return 1;
    }
    //cout << origin[0] << " " << origin[1] << " " << origin[2] << endl;
    
    getline(file, tempstring, ' ');
    if (tempstring.compare("SPACING") == 0) {
        istringstream iss;
        getline(file, tempstring);
        iss.str(tempstring);
        iss >> spacing[0] >> spacing[1] >> spacing[2];
        
    } else {
        cout << "No spacing info: " << endl;
        return 1;
    }
    //cout << spacing[0] << " " << spacing[1] << " " << spacing[2] << endl;
    getline(file, tempstring, ' ');
    if (tempstring.compare("CELL_DATA") == 0) {
        istringstream iss;
        getline(file, tempstring);
        iss.str(tempstring);
        iss >> n_CellData;
        if (n_CellData != dimensions[0]*dimensions[1]*dimensions[2]) {
            cout << "Nx*Ny*Nz = " << dimensions[0]*dimensions[1]*dimensions[2] << "!= n_CellData = " << n_CellData;
            return 1;
        }
    } else {
        cout << "No info about the number of CELL_DATA" << endl;
        return 1;
    }
    //cout << n_CellData << endl;
    
    //while (!file.eof()) {
        getline(file, tempstring, ' ');
        if (tempstring.compare("SCALARS") == 0) {
            // create a new CellData_Scaler
            cd_scaler.push_back(*new CellData_Scaler);
            // fetch data from file to it
            getline(file, cd_scaler.back().dataname, ' ');
            getline(file, tempstring);
            size_t ws_pos = tempstring.find_first_of(' ');
            // in case of the existence of numComp
            if (ws_pos != string::npos) {
                istringstream iss;
                iss.str(tempstring);
                iss >> cd_scaler.back().datatype >> std::ws >> cd_scaler.back().numcomp;
            } else {
                cd_scaler.back().datatype = tempstring;
                cd_scaler.back().numcomp = 1;
            }
            // final check
            getline(file, tempstring);
            if (tempstring.compare("LOOKUP_TABLE default") != 0) {
                cout << "Expected \"LOOKUP_TABLE default\", unsupportted file" << endl;
                return 1;
            }
            
            
            
        } else if (tempstring.compare("POINTS") == 0) {
            ;
        } else {
            cout << "No info about SCALARS or POINTS" << endl;
            return 1;
        }
    //}
    cout << cd_scaler.back().datatype << " " << cd_scaler.back().numcomp << endl;
    return 0;
}


