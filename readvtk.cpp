//
//  readvtk.cpp
//  readvtklis
//
//  Created by Rixin Li on 1/14/15.
//  Copyright (c) 2015 Rixin Li. All rights reserved.
//

#include "readvtk.h"

/*************************************/
/*********CellData_Scalar*************/
/*************************************/

/********** Overload != for dimension comparison **********/
inline bool operator!=(int (&dimensions)[3], CellData_Scalar &cd_scalar)
{
    if (dimensions[2] != cd_scalar.dimensions[2]) {
        return true;
    } else if (dimensions[1] != cd_scalar.dimensions[1]) {
        return true;
    } else if (dimensions[0] != cd_scalar.dimensions[0]) {
        return true;
    }
    return false;
}

/********** Constructor **********/
CellData_Scalar::CellData_Scalar(int *dimensions_been_told)
{
    Initialize_Data(dimensions_been_told);
}

/********** Destructor **********/
CellData_Scalar::~CellData_Scalar()
{
    ;
}

/********** Construct data **********/
int CellData_Scalar::Initialize_Data(int *dimensions_been_told)
{
    dimensions[2] = dimensions_been_told[2];
    dimensions[1] = dimensions_been_told[1];
    dimensions[0] = dimensions_been_told[0];
    data = NULL;
    // store in x, then y, then stack them as z direction
    data = new float**[dimensions[2]];
    for (int i = 0; i != dimensions[2]; i++) {
        data[i] = new float*[dimensions[1]];
        for (int j = 0; j != dimensions[1]; j++) {
            data[i][j] = new float[dimensions[0]];
        }
    }
    data[0][0][0] = 0;
    return 0;
}

/********** Read scalar data **********/
int CellData_Scalar::Read_Scalar_Data(string filename)
{
    ifstream file (filename.c_str(), ios::binary);
    file.seekg(pos, ios::beg);
    //FILE *fp;
    //fp = fopen(filename.c_str(), "r");
    //fseek(fp, pos, SEEK_SET);
    for (int i = 0; i != dimensions[2]; i++) {
        for (int j = 0; j != dimensions[1]; j++) {
            for (int k = 0; k != dimensions[0]; k++) {
                unsigned char temp[sizeof(float)];
                file.read(reinterpret_cast<char *>(temp), sizeof(float));
                unsigned char t;
                t = temp[0];
                temp[0] = temp[3];
                temp[3] = t;
                t = temp[1];
                temp[1] = temp[2];
                temp[2] = t;
                data[i][j][k] = reinterpret_cast<float&>(temp);
                // You can read it directly due to Endianness
                // Big Endian should be converted to samll endian
                //file.read((char *)(&data[i][j][k]), sizeof(float));
                //fread(&data[i][j][k], sizeof(float), 1, fp);
            }
        }
    }
    file.close();
    //fclose(fp);
    return 0;
}

/********** Free data memory **********/
int CellData_Scalar::Free_Data()
{
    for (int i = 0; i != dimensions[2]; i++) {
        for (int j = 0; j != dimensions[1]; j++) {
            delete [] data[i][j];
        }
        delete [] data[i];
    }
    delete [] data;
    dimensions[0] = 0;
    dimensions[1] = 0;
    dimensions[2] = 0;
    return 0;
}

/*************************************/
/*********CellData_Vector*************/
/*************************************/

/********** Overload != for dimension comparison **********/
inline bool operator!=(int (&dimensions)[3], CellData_Vector &cd_vector)
{
    if (dimensions[2] != cd_vector.dimensions[2]) {
        return true;
    } else if (dimensions[1] != cd_vector.dimensions[1]) {
        return true;
    } else if (dimensions[0] != cd_vector.dimensions[0]) {
        return true;
    }
    return false;
}

/********** Constructor **********/
CellData_Vector::CellData_Vector(int *dimensions_been_told)
{
    Initialize_Data(dimensions_been_told);
}

/********** Destructor **********/
CellData_Vector::~CellData_Vector()
{
    ;
}

/********** Construct data **********/
int CellData_Vector::Initialize_Data(int *dimensions_been_told)
{
    dimensions[2] = dimensions_been_told[2];
    dimensions[1] = dimensions_been_told[1];
    dimensions[0] = dimensions_been_told[0];
    data = NULL;
    // store in x, y and then stack them in z direction
    data = new float***[dimensions[2]];
    for (int i = 0; i != dimensions[2]; i++) {
        data[i] = new float**[dimensions[1]];
        for (int j = 0; j != dimensions[1]; j++) {
            data[i][j] = new float*[dimensions[0]];
            for (int k = 0; k != dimensions[0]; k++) {
                data[i][j][k] = new float[3];
            }
        }
    }
    return 0;
}

/********** Free data **********/
int CellData_Vector::Free_Data()
{
    for (int i = 0; i != dimensions[2]; i++) {
        for (int j = 0; j != dimensions[1]; j++) {
            for (int k = 0; k != dimensions[0]; k++) {
                delete [] data[i][j][k];
            }
            delete [] data[i][j];
        }
        delete [] data[i];
    }
    delete [] data;
    dimensions[0] = 0;
    dimensions[1] = 0;
    dimensions[2] = 0;
    return 0;
}

/********** Read scalar data **********/
int CellData_Vector::Read_Vector_Data(string filename)
{
    ifstream file (filename.c_str(), ios::binary);
    file.seekg(pos, ios::beg);
    for (int i = 0; i != dimensions[2]; i++) {
        for (int j = 0; j != dimensions[1]; j++) {
            for (int k = 0; k != dimensions[0]; k++) {
                for (int l = 0; l != 3; l++) {
                    unsigned char temp[sizeof(float)];
                    file.read(reinterpret_cast<char *>(temp), sizeof(float));
                    unsigned char t;
                    t = temp[0];
                    temp[0] = temp[3];
                    temp[3] = t;
                    t = temp[1];
                    temp[1] = temp[2];
                    temp[2] = t;
                    data[i][j][k][l] = reinterpret_cast<float&>(temp);
                }
                //file.read((char *)(data[i][j][k]), 3*sizeof(float));
            }
        }
    }
    file.close();
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
    if (cd_scalar.size() > 0) {
        for (int i = 0; i != cd_scalar.size(); i++) {
            cd_scalar[i].Free_Data();
        }
        vector<CellData_Scalar> temp;
        cd_scalar.swap(temp);
    }
    if (cd_vector.size() > 0) {
        for (int i = 0; i != cd_vector.size(); i++) {
            cd_vector[i].Free_Data();
        }
        vector<CellData_Vector> temp;
        cd_vector.swap(temp);
    }
    if (cell_center != NULL) {
        for (int i = 0; i != dimensions[2]; i++) {
            for (int j = 0; j != dimensions[1]; j++) {
                for (int k = 0; k != dimensions[0]; k++) {
                    delete [] cell_center[i][j][k];
                }
                delete [] cell_center[i][j];
            }
            delete [] cell_center[i];
        }
        delete [] cell_center;
    }
    
}

/********** Construct coordinate grid **********/
int VtkFile::Construct_Coor()
{
    // actually, here what we construct is the center coordinate of each cell, along [z][y][x], from ORIGIN
    
    cell_center = new double***[dimensions[2]];
    for (int i = 0; i != dimensions[2]; i++) {
        cell_center[i] = new double**[dimensions[1]];
        for (int j = 0; j != dimensions[1]; j++) {
            cell_center[i][j] = new double*[dimensions[0]];
            for (int k = 0; k != dimensions[0]; k++) {
                cell_center[i][j][k] = new double[3];
                cell_center[i][j][k][0] = origin[0]+spacing[0]*(k+0.5);
                cell_center[i][j][k][1] = origin[1]+spacing[1]*(j+0.5);
                cell_center[i][j][k][2] = origin[2]+spacing[2]*(i+0.5);
            }
        }
    }
    return 0;
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
    //cout << "File format: " << FileFormat << endl;
    
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
    //cout << "DIMENSIONS " << dimensions[0] << " " << dimensions[1] << " " << dimensions[2] << endl;
    
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
    //cout << "ORIGIN " << origin[0] << " " << origin[1] << " " << origin[2] << endl;
    
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
    //cout << "SPACING " << spacing[0] << " " << spacing[1] << " " << spacing[2] << endl;
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
    //cout << "CELL_DATA " << n_CellData << endl;
    
    int n_cd_scalar = 0, n_cd_vector = 0;
    //long filepos1, filepos2;
    while (!file.eof()) {
        getline(file, tempstring, ' ');
        if (tempstring.compare("SCALARS") == 0) {
            n_cd_scalar++;
            // if vector has elements, no need to push_back a new one
            if (cd_scalar.size() >= n_cd_scalar) {
                if (dimensions != cd_scalar[n_cd_scalar-1]) {
                    cd_scalar[n_cd_scalar-1].Initialize_Data(dimensions);
                }
            } else {
                // create a new CellData_Scalar
                cd_scalar.push_back(CellData_Scalar(dimensions));
            }
            // fetch data from file to it
            getline(file, cd_scalar[n_cd_scalar-1].dataname, ' ');
            getline(file, tempstring);
            size_t ws_pos = tempstring.find_first_of(' ');
            // in case of the existence of numComp
            if (ws_pos != string::npos) {
                istringstream iss;
                iss.str(tempstring);
                iss >> cd_scalar[n_cd_scalar-1].datatype >> std::ws >> cd_scalar[n_cd_scalar-1].numcomp;
            } else {
                cd_scalar[n_cd_scalar-1].datatype = tempstring;
                cd_scalar[n_cd_scalar-1].numcomp = 1;
            }
            if (cd_scalar[n_cd_scalar-1].datatype.compare("float") != 0) {
                cout << "Expected float format, found " << tempstring << endl;
                return 1;
            }
            
            // final check
            getline(file, tempstring);
            if (tempstring.compare("LOOKUP_TABLE default") != 0) {
                cout << "Expected \"LOOKUP_TABLE default\", unsupportted file" << endl;
                return 1;
            }
            
            cd_scalar[n_cd_scalar-1].tablename = "default";
            //cout << "Found scalar " << cd_scalar[n_cd_scalar-1].dataname << " " << cd_scalar[n_cd_scalar-1].datatype << endl;
            cd_scalar[n_cd_scalar-1].pos = file.tellg();
            file.seekg(sizeof(float)*n_CellData, file.cur);
            // for debug
            /*
            filepos1 = file.tellg();
            file.seekg(0, ios::end);
            filepos2 = file.tellg();
            cout << "distance to eof: " << filepos2 - filepos1 << "bytes" << endl;
            file.seekg(filepos1, ios::beg);
             */
        } else if (tempstring.compare("VECTORS") == 0) {
            n_cd_vector++;
            // if vector has elements, no need to push_back a new one
            if (cd_vector.size() >= n_cd_vector) {
                if (dimensions != cd_vector[n_cd_vector-1]) {
                    cd_vector[n_cd_vector-1].Initialize_Data(dimensions);
                }
            } else {
                // create a new CellData_Scalar
                cd_vector.push_back(CellData_Vector(dimensions));
            }
            // fetch data from file to it
            getline(file, cd_vector[n_cd_vector-1].dataname, ' ');
            getline(file, cd_vector[n_cd_vector-1].datatype);
            //cout << "Found vector " << cd_vector[n_cd_vector-1].dataname << " " << cd_vector[n_cd_vector-1].datatype << endl;
            
            cd_vector[n_cd_vector-1].pos = file.tellg();
            //cout << cd_vector[n_cd_vector-1].pos << endl;
            file.seekg(sizeof(float)*n_CellData*3, file.cur);
            // for debug
            /*
            filepos1 = file.tellg();
            file.seekg(0, ios::end);
            filepos2 = file.tellg();
            cout << "distance to eof: " << filepos2 - filepos1 << "bytes" << endl;
            file.seekg(filepos1, ios::beg);
             */
        } else {
            // it seems vtk file has a new empty line in the end
            if (tempstring.length() != 0 ) {
                cout << "No info about SCALARS or VECTORS, it is " << tempstring << tempstring.length() << endl;
                return 1;
            }
        }
    }
    file.close();

    return 0;
}

/********** Read data **********/
int VtkFile::Read_Data(string filename)
{
    for (int i = 0; i != cd_scalar.size(); i++) {
        cd_scalar[i].Read_Scalar_Data(filename);
    }
    for (int i = 0; i != cd_vector.size(); i++) {
        cd_vector[i].Read_Vector_Data(filename);
    }
    return 0;
}

/********** Print file info **********/
int VtkFile::Print_File_Info()
{
    if (Version.length() == 0) {
        cout << "No file info for now" << endl;
        return 1;
    }
    cout << "\nVersion: " << Version << endl;
    cout << setprecision(6) << scientific << "time: " << time << endl;
    cout << "File format: " << FileFormat << endl;
    cout << DatasetStructure << endl;
    cout << "DIMENSIONS " << dimensions[0] << " " << dimensions[1] << " " << dimensions[2] << endl;
    cout << "ORIGIN " << origin[0] << " " << origin[1] << " " << origin[2] << endl;
    cout << "SPACING " << spacing[0] << " " << spacing[1] << " " << spacing[2] << endl;
    cout << "CELL_DATA " << n_CellData << endl;
    return 0;
}

/********** Calculate mass and find maximum **********/
int VtkFile::Calculate_Mass_Find_Max()
{
    // here we use the specific features of Athena output
    // i.e. the first scalar is density, and then  particle_density
    // and we assume the volume of one spacing cube is 1, so
    // the total mass = total density*1
    if (cd_scalar.size() != 2) {
        cout << "The size of Scalar vector is wrong." << endl;
        return 1;
    }
    for (vector<CellData_Scalar>::iterator it = cd_scalar.begin(); it != cd_scalar.end(); it++) {
        double m_temp = 0, maximum = 0, tempdata;
        for (int k = 0; k != dimensions[2]; k++) {
            for (int j = 0; j != dimensions[1]; j++) {
                for (int i = 0; i != dimensions[0]; i++) {
                    tempdata = it->data[i][j][k];
                    m_temp += tempdata;
                    if (tempdata > maximum) {
                        maximum = tempdata;
                    }
                }
            }
        }
        if (it->dataname.compare("density") == 0) {
            m_gas = m_temp;
            max_rho_gas = maximum;
        } else if (it->dataname.compare("particle_density") == 0) {
            m_par = m_temp;
            max_rho_par = maximum;
        } else {
            cout << "Unkonwn data name: " << it->dataname << endl;
        }
    }
    return 0;
}

