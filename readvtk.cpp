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
/*! \fn int Initialize_Data(int *dimensions_been_told)
 *  \brief initialize data */
int CellData_Scalar::Initialize_Data(int *dimensions_been_told)
{
    dimensions[2] = dimensions_been_told[2];
    dimensions[1] = dimensions_been_told[1];
    dimensions[0] = dimensions_been_told[0];
    
    data = allocate3d_scalar_array<float>(dimensions);

    return 0;
}

/********** Read scalar data **********/
/*! \fn int Read_Scalar_Data(string filename)
 *  \brief read scalar data */
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
/*! \fn int Free_Data()
 *  \brief free the data memory */
int CellData_Scalar::Free_Data()
{
    deallocate3d_scalar_array(data, dimensions);
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
/*! \fn int Initialize_Data(int *dimensions_been_told)
 *  \brief initialize data */
int CellData_Vector::Initialize_Data(int *dimensions_been_told)
{
    dimensions[2] = dimensions_been_told[2];
    dimensions[1] = dimensions_been_told[1];
    dimensions[0] = dimensions_been_told[0];
    data = NULL;
    
    data = allocate3d_vector_array<float>(dimensions);
    
    return 0;
}

/********** Free data **********/
/*! \fn int Free_Data()
 *  \brief free the data memory */
int CellData_Vector::Free_Data()
{
    deallocate3d_vector_array(data, dimensions);
    
    dimensions[0] = 0;
    dimensions[1] = 0;
    dimensions[2] = 0;
    return 0;
}

/********** Read scalar data **********/
/*! \fn int Read_Vector_Data(string filename)
 *  \brief read vector data */
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

/********** Overload != for dimension comparison **********/
inline bool operator!=(int (&dimensions)[3], VtkFile &vf)
{
    if (dimensions[2] != vf.dimensions[2]) {
        return true;
    } else if (dimensions[1] != vf.dimensions[1]) {
        return true;
    } else if (dimensions[0] != vf.dimensions[0]) {
        return true;
    }
    return false;
}

/********** Constructor **********/
VtkFile::VtkFile()
{
    ;
}

/********** Destructor **********/
VtkFile::~VtkFile()
{
    /*
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
    deallocate3d_vector_array(cell_center, dimensions);
     */
    
}

/********** Construct coordinate grid **********/
/*! \fn int Construct_Coor(int *dimensions_been_told)
 *  \brief consturct coordinate grid */
int VtkFile::Construct_Coor(int *dimensions_been_told)
{
    cell_center = allocate3d_vector_array<double>(dimensions);
    // actually, here what we construct is the center coordinate of each cell, along [z][y][x], from ORIGIN
    
    cell_center = allocate3d_vector_array<double>(dimensions);
    for (int i = 0; i != dimensions[2]; i++) {
        for (int j = 0; j != dimensions[1]; j++) {
            for (int k = 0; k != dimensions[0]; k++) {
                cell_center[i][j][k][0] = origin[0]+spacing[0]*(k+0.5);
                cell_center[i][j][k][1] = origin[1]+spacing[1]*(j+0.5);
                cell_center[i][j][k][2] = origin[2]+spacing[2]*(i+0.5);
            }
        }
    }
    //if (fio->paras.ccz == NULL) {
    fio->paras.ccz = new double[dimensions[2]];
    //}
    for (int k = 0; k != dimensions[2]; k++) {
        fio->paras.ccz[k] = origin[2]+spacing[2]*(k+0.5);
    }
    //if (fio->paras.ccy == NULL) {
    fio->paras.ccy = new double[dimensions[1]];
    //}
    for (int j = 0; j != dimensions[1]; j++) {
        fio->paras.ccy[j] = origin[1]+spacing[1]*(j+0.5);
    }
    //if (fio->paras.ccx == NULL) {
    fio->paras.ccx = new double[dimensions[0]];
    //}
    for (int i = 0; i != dimensions[0]; i++) {
        fio->paras.ccx[i] = origin[0]+spacing[0]*(i+0.5);
    }
    
    return 0;
}

/********** Read header and record data position **********/
/*! \fn int Read_Header_Record_Pos(string filename)
 *  \brief read header and record data position */
int VtkFile::Read_Header_Record_Pos(string filename)
{
    Sigma_gas_0 = 2.50662827463100050241576528481;
    string tempstring;
    ifstream file (filename.c_str(), ios::in);
    if (!file.is_open()) {
        cout << "Failed to open " << filename << "\n";
        return 1;
    }
    
    getline(file, Version);
    if (Version.compare("# vtk DataFile Version 3.0") != 0 && Version.compare("# vtk DataFile Version 2.0") != 0) {
        cout << "First line of " << filename << " is " << Version << "\n";
        return 1;
    }
    //cout << "Version: " << Version << "\n";
    
    getline(file, Header);
    if (Header.find("CONSERVED") != string::npos) {
        size_t time_pos = Header.find("time= ");
        // stod() reads to the end of the number, so we can ignore those afterwards
        // time = stod(Header.substr(time_pos+6));
        time = strtod((Header.substr(time_pos+6, 12)).c_str(), NULL);
    }
    //cout << setprecision(6) << scientific << "time: " << time << "\n";
    
    getline(file, FileFormat);
    if (FileFormat.compare("BINARY") != 0) {
        cout << "Unsopported file format: " << FileFormat << "\n";
        return 1;
    }
    //cout << "File format: " << FileFormat << "\n";
    
    getline(file, DatasetStructure);
    if (DatasetStructure.compare("DATASET STRUCTURED_POINTS") != 0) {
        cout << "Unsopported dataset structure: " << DatasetStructure << "\n";
        return 1;
    }
    //cout << DatasetStructure << "\n";
    
    getline(file, tempstring, ' ');
    if (tempstring.compare("DIMENSIONS") == 0) {
        istringstream iss;
        getline(file, tempstring);
        iss.str(tempstring);
        iss >> dimensions[0] >> dimensions[1] >> dimensions[2];
        //We want to store the number of grid cells, not the number of grid cell corners
        for (int i = 0; i != 3; i++) {
            dimensions[i]--;
            fio->paras.dimensions[i] = dimensions[i];
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
        cout << "No dimensions info: " << "\n";
        return 1;
    }
    //cout << "DIMENSIONS " << dimensions[0] << " " << dimensions[1] << " " << dimensions[2] << "\n";
    
    getline(file, tempstring, ' ');
    if (tempstring.compare("ORIGIN") == 0) {
        istringstream iss;
        getline(file, tempstring);
        iss.str(tempstring);
        iss >> origin[0] >> origin[1] >> origin[2];
        
    } else {
        cout << "No origin info: " << "\n";
        return 1;
    }
    //cout << "ORIGIN " << origin[0] << " " << origin[1] << " " << origin[2] << "\n";
    
    getline(file, tempstring, ' ');
    if (tempstring.compare("SPACING") == 0) {
        istringstream iss;
        getline(file, tempstring);
        iss.str(tempstring);
        iss >> spacing[0] >> spacing[1] >> spacing[2];
        cell_volume = spacing[0] * spacing[1] * spacing[2];
        fio->paras.spacing[0] = spacing[0];
        fio->paras.spacing[1] = spacing[1];
        fio->paras.spacing[2] = spacing[2];
    } else {
        cout << "No spacing info: " << "\n";
        return 1;
    }
    //cout << "SPACING " << spacing[0] << " " << spacing[1] << " " << spacing[2] << "\n";
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
        cout << "No info about the number of CELL_DATA" << "\n";
        return 1;
    }
    //cout << "CELL_DATA " << n_CellData << "\n";
    
    int n_cd_scalar = 0, n_cd_vector = 0;
    if (cell_center == NULL) {
        Construct_Coor(dimensions);
    }
    //long filepos1, filepos2;
    while (!file.eof()) {
        getline(file, tempstring, ' ');
        if (tempstring.compare("SCALARS") == 0) {
            n_cd_scalar++;// cout << "n_cd_scalar=" << n_cd_scalar << "\n";
            // if vector has elements, no need to push_back a new one
            if (cd_scalar.size() >= n_cd_scalar) {
                if (dimensions != cd_scalar[n_cd_scalar-1]) {
                    cd_scalar[n_cd_scalar-1].Free_Data();
                    cd_scalar[n_cd_scalar-1].Initialize_Data(dimensions);
                    //cout << "Real new data" << "\n";
                }
                //cout << "Initialize_Data" << "\n";
            } else {
                // create a new CellData_Scalar
                cd_scalar.push_back(CellData_Scalar(dimensions));
                //cout << "Push_back" << "\n";
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
                cout << "Expected float format, found " << tempstring << "\n";
                return 1;
            }
            
            // final check
            getline(file, tempstring);
            if (tempstring.compare("LOOKUP_TABLE default") != 0) {
                cout << "Expected \"LOOKUP_TABLE default\", unsupportted file" << "\n";
                return 1;
            }
            
            cd_scalar[n_cd_scalar-1].tablename = "default";
            //cout << "Found scalar " << cd_scalar[n_cd_scalar-1].dataname << " " << cd_scalar[n_cd_scalar-1].datatype << "\n";
            cd_scalar[n_cd_scalar-1].pos = file.tellg();
            file.seekg(sizeof(float)*n_CellData, file.cur);
            // for debug
            /*
            filepos1 = file.tellg();
            file.seekg(0, ios::end);
            filepos2 = file.tellg();
            cout << "distance to eof: " << filepos2 - filepos1 << "bytes" << "\n";
            file.seekg(filepos1, ios::beg);
             */
        } else if (tempstring.compare("VECTORS") == 0) {
            n_cd_vector++;
            // if vector has elements, no need to push_back a new one
            if (cd_vector.size() >= n_cd_vector) {
                if (dimensions != cd_vector[n_cd_vector-1]) {
                    cd_vector[n_cd_vector-1].Free_Data();
                    cd_vector[n_cd_vector-1].Initialize_Data(dimensions);
                }
            } else {
                // create a new CellData_Scalar
                cd_vector.push_back(CellData_Vector(dimensions));
            }
            // fetch data from file to it
            getline(file, cd_vector[n_cd_vector-1].dataname, ' ');
            getline(file, cd_vector[n_cd_vector-1].datatype);
            //cout << "Found vector " << cd_vector[n_cd_vector-1].dataname << " " << cd_vector[n_cd_vector-1].datatype << "\n";
            
            cd_vector[n_cd_vector-1].pos = file.tellg();
            //cout << cd_vector[n_cd_vector-1].pos << "\n";
            file.seekg(sizeof(float)*n_CellData*3, file.cur);
            // for debug
            /*
            filepos1 = file.tellg();
            file.seekg(0, ios::end);
            filepos2 = file.tellg();
            cout << "distance to eof: " << filepos2 - filepos1 << "bytes" << "\n";
            file.seekg(filepos1, ios::beg);
             */
        } else {
            // it seems vtk file has a new empty line in the end
            if (tempstring.length() != 0 ) {
                cout << "No info about SCALARS or VECTORS, it is " << tempstring << tempstring.length() << "\n";
                return 1;
            }
        }
    }
    file.close();

    return 0;
}

/********** Read data **********/
/*! \fn int Read_Data(string filename)
 *  \brief read data */
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
/*! \fn int Print_File_Info()
 *  \brief print file info */
int VtkFile::Print_File_Info()
{
    if (Version.length() == 0) {
        cout << "No file info for now" << "\n";
        return 1;
    }
    cout << "\nVersion: " << Version << "\n";
    cout << setprecision(6) << scientific << "time: " << time << "\n";
    cout << "File format: " << FileFormat << "\n";
    cout << DatasetStructure << "\n";
    cout << "DIMENSIONS " << dimensions[0] << " " << dimensions[1] << " " << dimensions[2] << "\n";
    cout << "ORIGIN " << origin[0] << " " << origin[1] << " " << origin[2] << "\n";
    cout << "SPACING " << spacing[0] << " " << spacing[1] << " " << spacing[2] << "\n";
    cout << "CELL_DATA " << n_CellData << "\n";
    return 0;
}

/********** Calculate mass and find maximum **********/
/*! \fn int Calculate_Mass_Find_Max()
 \brief calculate mass and find maximum and mass loss */
int VtkFile::Calculate_Mass_Find_Max()
{
    // here we use the specific features of Athena output
    // i.e. the first scalar is density, and then  particle_density
    // and we assume the volume of one spacing cube is 1, so
    // the total mass = total density*1
    if (cd_scalar.size() != 2) {
        cout << "The size of Scalar vector is wrong." << "\n";
        return 1;
    }
    for (vector<CellData_Scalar>::iterator it = cd_scalar.begin(); it != cd_scalar.end(); it++) {
        double m_temp = 0, maximum = 0, tempdata;
        for (int k = 0; k != dimensions[2]; k++) {
            for (int j = 0; j != dimensions[1]; j++) {
                for (int i = 0; i != dimensions[0]; i++) {
                    tempdata = it->data[k][j][i];
                    m_temp += tempdata;
                    if (tempdata > maximum) {
                        maximum = tempdata;
                    }
                }
            }
        }
        if (it->dataname.compare("density") == 0) {
            m_gas = m_temp*cell_volume;
            Max_Rhog = maximum;
        } else if (it->dataname.compare("particle_density") == 0) {
            m_par = m_temp*cell_volume;
            Max_Rhop = maximum;
        } else {
            cout << "Unkonwn data name: " << it->dataname << "\n";
        }
    }
    
    for (vector<CellData_Vector>::iterator it = cd_vector.begin(); it != cd_vector.end(); it++) {
        if (it->dataname.compare("momentum") == 0) {
            dSigma = 0; //int inflow_count = 0;
            for (int iy = 0; iy != dimensions[1]; iy++) {
                for (int ix = 0; ix != dimensions[0]; ix++) {
                    if (cd_vector[0].data[dimensions[2]-1][iy][ix][2]  > 0) {
                        dSigma += cd_vector[0].data[dimensions[2]-1][iy][ix][2];
                    } else {
                        dSigma += 0.5 * cd_vector[0].data[dimensions[2]-1][iy][ix][2];
                        //cout << "inflow: ix=" << ix << ", iy=" << iy << "\n";
                        //inflow_count++;
                    }
                    if (cd_vector[0].data[0][iy][ix][2] < 0) {
                        dSigma -= cd_vector[0].data[0][iy][ix][2];
                    } else {
                        dSigma -= 0.5 * cd_vector[0].data[0][iy][ix][2];
                        //cout << "inflow: ix=" << ix << ", iy=" << iy << "\n";
                        //inflow_count++;
                    }
                }
            }
            //cout << "in/out=" << float(inflow_count)/(2*dimensions[0]*dimensions[1]) << "\n";
            dSigma /= (dimensions[0]*dimensions[1]);
        }
    }
    return 0;
}

/********** Calculate gas peculiar velocity **********/
/*! \fn int VpecG(double *VpecGx)
 *  \brief return gas peculiar velocity components averaged horizontally at each z, weighted by rho_g */
int VtkFile::VpecG(double *VpecG)
{
    double temp_rhog = 0, temp_momentum[3];
    for (int iz = 0; iz != dimensions[2]; iz++) {
        temp_rhog = 0;
        temp_momentum[0] = 0; temp_momentum[1] = 0; temp_momentum[2] = 0;
        for (int iy = 0; iy != dimensions[1]; iy++) {
            for (int ix = 0; ix != dimensions[0]; ix++) {
                temp_rhog += cd_scalar[0].data[iz][iy][ix];
                temp_momentum[0] += cd_vector[0].data[iz][iy][ix][0];
                temp_momentum[1] += cd_vector[0].data[iz][iy][ix][1];
                temp_momentum[2] += cd_vector[0].data[iz][iy][ix][2];
            }
        }
        VpecG[iz] = temp_momentum[0]/temp_rhog;
        VpecG[iz+dimensions[2]] = temp_momentum[1]/temp_rhog;
        VpecG[iz+2*dimensions[2]] = temp_momentum[2]/temp_rhog;
    }
    return 0;
}

/********** Calculate mean column density **********/
/*! \fn int MeanSigma(double *MeanSigma)
 *  \brief calculate sigma_g and sigma_p averaged over y */
int VtkFile::MeanSigma(double *MeanSigma)
{
    for (int ix = 0; ix != dimensions[0]; ix++) {
        MeanSigma[ix] = 0;
        MeanSigma[ix+dimensions[0]] = 0;
        for (int iz = 0; iz != dimensions[2]; iz++) {
            for (int iy = 0; iy != dimensions[1]; iy++) {
                MeanSigma[ix] += cd_scalar[0].data[iz][iy][ix];
                MeanSigma[ix+dimensions[0]] += cd_scalar[1].data[iz][iy][ix] * cell_volume;
            }
        }
        MeanSigma[ix] /= Sigma_gas_0_inbox;
        MeanSigma[ix+dimensions[0]] /= (Sigma_gas_0 * spacing[0] * spacing[1] * dimensions[1]);
    }
    return 0;
}

/********** Calculate vertical rho_g and rho_p **********/
/*! \fn int VertRho(double *VertRho)
 *  \brief calculate vertical rho_g and rho_p */
int VtkFile::VertRho(double *VertRho)
{
    for (int iz = 0; iz != dimensions[2]; iz++) {
        VertRho[iz] = 0;
        VertRho[iz+dimensions[2]] = 0;
        for (int iy = 0; iy != dimensions[1]; iy++) {
            for (int ix = 0; ix != dimensions[0]; ix++) {
                VertRho[iz] += cd_scalar[0].data[iz][iy][ix];
                VertRho[iz+dimensions[2]] += cd_scalar[1].data[iz][iy][ix];
            }
        }
        VertRho[iz] /= dimensions[0]*dimensions[1];
        VertRho[iz+dimensions[2]] /= dimensions[0]*dimensions[1];
    }
    return 0;
}

/********** Calculate correlation length **********/
/*! \fn int CorrLen(double *CoorL)
 *  \brief calculate the correlation length */
int VtkFile::CorrLen(double *CorrL)
{
    int l = 0, lcMx, lcVx, lcRho;
    double MxMx, RhoRho;
    double *corrMx = new double[dimensions[0]];
    double *corrVx = new double[dimensions[0]];
    double *corrRho = new double[dimensions[0]];
    for (int iz = 0; iz != dimensions[2]; iz++) {
        // Initilization
        CorrL[iz] = 0;
        CorrL[iz+dimensions[2]] = 0;
        CorrL[iz+2*dimensions[2]] = 0;
        lcMx = dimensions[0]-1;
        lcVx = lcMx;
        lcRho = lcMx;
        
        for (l = 0; l != dimensions[0]; l ++) {
            corrMx[l] = 0;
            corrVx[l] = 0;
            corrRho[l] = 0;
            for (int iy = 0; iy != dimensions[1]; iy++) {
                for (int ix = 0; ix != dimensions[0]-l; ix++) {
                    MxMx = (cd_vector[0].data[iz][iy][ix][0] * cd_vector[0].data[iz][iy][ix+l][0]);
                    corrMx[l] += MxMx;
                    RhoRho = (cd_scalar[0].data[iz][iy][ix] * cd_scalar[0].data[iz][iy][ix+l]);
                    corrRho[l] += RhoRho;
                    corrVx[l] += MxMx/RhoRho;
                }
            }
        }
        for (l = dimensions[0]-1; l >= 0; l--) {
            corrMx[l] /= corrMx[0];
            corrVx[l] /= corrVx[0];
            corrRho[l] /= corrRho[0];
            if (fabs(corrMx[l]-0.5) < fabs(corrMx[lcMx]-0.5)) {
                lcMx = l;
            }
            if (fabs(corrVx[l]-0.5) < fabs(corrVx[lcVx]-0.5)) {
                lcVx = l;
            }
            if (fabs(corrRho[l]-0.5) < fabs(corrRho[lcRho]-0.5)) {
                lcRho = l;
            }
        }
        CorrL[iz] = lcMx * spacing[0];
        CorrL[iz+dimensions[2]] = lcVx * spacing[0];
        CorrL[iz+2*dimensions[2]] = lcRho * spacing[0];
    }
    
    delete [] corrMx;
    delete [] corrVx;
    delete [] corrRho;
    return 0;
}







