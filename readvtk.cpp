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
    cell_center = NULL;
    GPME = NULL;
    GPMEPar = NULL;
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
    cell_center = allocate3d_vector_array<float>(dimensions);
    // actually, here what we construct is the center coordinate of each cell, along [z][y][x], from ORIGIN
    
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
    fio->paras.ccz = new float[dimensions[2]];
    //}
    for (int k = 0; k != dimensions[2]; k++) {
        fio->paras.ccz[k] = origin[2]+spacing[2]*(k+0.5);
    }
    //if (fio->paras.ccy == NULL) {
    fio->paras.ccy = new float[dimensions[1]];
    //}
    for (int j = 0; j != dimensions[1]; j++) {
        fio->paras.ccy[j] = origin[1]+spacing[1]*(j+0.5);
    }
    //if (fio->paras.ccx == NULL) {
    fio->paras.ccx = new float[dimensions[0]];
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
        // time = stod(Header.substr(time_pos+6));
        time = strtod((Header.substr(time_pos+6, 12)).c_str(), NULL);
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
        cell_volume = spacing[0] * spacing[1] * spacing[2];
        fio->paras.spacing[0] = spacing[0];
        fio->paras.spacing[1] = spacing[1];
        fio->paras.spacing[2] = spacing[2];
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
    if (cell_center == NULL) { //cout << myMPI->myrank << " is here." << endl;
        Construct_Coor(dimensions);
    }
    kps = dimensions[2]/2 - int(round(0.1/spacing[2]));
    kpe = dimensions[2]/2 + int(round(0.1/spacing[2]));
    for (int i = 0; i != 3; i++) {
        L[i] = dimensions[i] * spacing[i];
    }
    
    if (GPME == NULL) {
        GPME = new double[4*dimensions[2]];
        for (int i = 0; i != 4*dimensions[2]; i++) {
            GPME[i] = 0;
        }
        //cout << "GPME initial allocation done!\n";
    }
    if (GPMEPar == NULL) {
        GPMEPar = new double[4*dimensions[2]];
        for (int i = 0; i != 4*dimensions[2]; i++) {
            GPMEPar[i] = 0;
        }
        //cout << "GPMEPar initial allocation done!\n";
    }
    
    //long filepos1, filepos2;
    while (!file.eof()) {
        getline(file, tempstring, ' ');
        if (tempstring.compare("SCALARS") == 0) {
            n_cd_scalar++;// cout << "n_cd_scalar=" << n_cd_scalar << endl;
            // if vector has elements, no need to push_back a new one
            if (cd_scalar.size() >= n_cd_scalar) {
                if (dimensions != cd_scalar[n_cd_scalar-1]) {
                    cd_scalar[n_cd_scalar-1].Free_Data();
                    cd_scalar[n_cd_scalar-1].Initialize_Data(dimensions);
                    //cout << "Real new data" << endl;
                }
                //cout << "Initialize_Data" << endl;
            } else {
                // create a new CellData_Scalar
                cd_scalar.push_back(CellData_Scalar(dimensions));
                //cout << "Push_back" << endl;
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
/*! \fn int Calculate_Mass_Find_Max()
 \brief calculate mass and find maximum and mass loss */
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
    
    /* RL: move these calculations to octree.cpp
    for (vector<CellData_Scalar>::iterator it = cd_scalar.begin(); it != cd_scalar.end(); it++) {
        if (it->dataname.compare("density") == 0) {
            //m_gas = m_temp*cell_volume;
            //Max_Rhog = maximum;
        } else if (it->dataname.compare("particle_density") == 0) {
            float m_temp = 0, maximum = 0, tempdata;
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
            m_par = m_temp*cell_volume;
            Max_Rhop = maximum;
            
            RpAV = 0; RpSQ = 0; RpQU = 0;
            float tempSQ = 0, tempV = 0;
            for (int k = kps; k != kpe; k++) {
                for (int j = 0; j != dimensions[1]; j++) {
                    for (int i = 0; i != dimensions[0]; i++) {
                        tempdata = it->data[k][j][i];
                        RpAV += tempdata;
                        tempSQ = tempdata * tempdata;
                        RpSQ += tempSQ;
                        RpQU += tempSQ * tempSQ;
                    }
                }
            }
            tempV = (kpe-kps)*dimensions[1]*dimensions[0];
            RpAV /= tempV;
            RpSQ = sqrt(RpSQ/tempV);
            RpQU = sqrt(sqrt(RpQU/tempV));
        } else {
            cout << "Unkonwn data name: " << it->dataname << endl;
        }
    } //*/
    
    //int inflow_count = 0; float inflow = 0.0, outflow = 0.0;
    for (vector<CellData_Vector>::iterator it = cd_vector.begin(); it != cd_vector.end(); it++) {
        if (it->dataname.compare("momentum") == 0) {
            dSigma = 0; //int inflow_count = 0;
            for (int iy = 0; iy != dimensions[1]; iy++) {
                for (int ix = 0; ix != dimensions[0]; ix++) {
#ifdef OutflowRate
                    if (cd_vector[0].data[dimensions[2]-1][iy][ix][2]  > 0) {
                        dSigma += cd_vector[0].data[dimensions[2]-1][iy][ix][2];
                        //outflow += abs(cd_vector[0].data[dimensions[2]-1][iy][ix][2]);
                    } else {
                        dSigma += 0.5 * cd_vector[0].data[dimensions[2]-1][iy][ix][2];
                        //cout << "inflow: ix=" << ix << ", iy=" << iy << endl;
                        //inflow += abs(cd_vector[0].data[dimensions[2]-1][iy][ix][2]);
                        //inflow_count++;
                    }
                    if (cd_vector[0].data[0][iy][ix][2] < 0) {
                        dSigma -= cd_vector[0].data[0][iy][ix][2];
                        //outflow += abs(cd_vector[0].data[0][iy][ix][2]);
                    } else {
                        dSigma -= 0.5 * cd_vector[0].data[0][iy][ix][2];
                        //cout << "inflow: ix=" << ix << ", iy=" << iy << endl;
                        //inflow += abs(cd_vector[0].data[0][iy][ix][2]);
                        //inflow_count++;
                    }
#endif // OutflowRate
#ifdef PeriodicFlux
                    dSigma += cd_vector[0].data[dimensions[2]-1][iy][ix][2];
#endif
                }
            }
#ifdef OutflowRate
            //cout << "in/out=" << float(inflow_count)/(2*dimensions[0]*dimensions[1]) << " " << inflow/outflow << endl;
#endif
            dSigma /= (dimensions[0]*dimensions[1]);
#ifdef PeriodicFlux
            dSigma = abs(dSigma);
#endif
        }
    }
    return 0;
}

/********** Calculate gas peculiar velocity **********/
/*! \fn int VpecG(float *VpecGx)
 *  \brief return gas peculiar velocity components averaged horizontally at each z, weighted by rho_g */
int VtkFile::VpecG(float *VpecG)
{
    float temp_rhog = 0, temp_momentum[3];
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
/*! \fn int MeanSigma(float *MeanSigma)
 *  \brief calculate sigma_g and sigma_p averaged over y */
int VtkFile::MeanSigma(float *MeanSigma)
{
    int mid[2], TwoNx;
    mid[0] = dimensions[2]/2-1;
    mid[1] = dimensions[2]/2;
    TwoNx = dimensions[0]*2;
    int Hp[2], ThreeNx;
    Hp[0] = dimensions[2]/2-8;
    Hp[1] = dimensions[2]/2+7;
    ThreeNx = dimensions[0]*3;
    
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

        // Since the MeanSigma for gas over entire height contains too much gas structures, we want to
        // check what if we only calculate the average near midplane
        for (int iy = 0; iy != dimensions[1]; iy++) {
            MeanSigma[ix+TwoNx] += (cd_scalar[0].data[mid[0]][iy][ix]+cd_scalar[0].data[mid[1]][iy][ix]);
        }
        for (int iz = Hp[0]; iz != Hp[1]; iz++) {
            for (int iy = 0; iy != dimensions[1]; iy++) {
                MeanSigma[ix+ThreeNx] += cd_scalar[0].data[iz][iy][ix];
            }
        }
        MeanSigma[ix+TwoNx] /= Sigma_gas_0_mid;
        MeanSigma[ix+ThreeNx] /= Sigma_gas_0_in2Hp;
    }
    return 0;
}

/********** Calculate vertical rho_g and rho_p **********/
/*! \fn int VertRho(float *VertRho)
 *  \brief calculate vertical rho_g and rho_p */
int VtkFile::VertRho(float *VertRho)
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
/*! \fn int CorrLen(float *CoorL, float *CorrV)
 *  \brief calculate the correlation length */
int VtkFile::CorrLen(float *CorrL
#ifdef CorrValue
                    , float *CorrV
#endif
                     )
{
    int Nx = dimensions[0], Ny = dimensions[1], Nz = dimensions[2];
    int l = 0, Nl = Nx/2+1, lcMx, lcVx, lcRho;
    float MxMx, RhoRho, MxBar, RhoBar, VxBar;
    float *corrMx = new float[Nl];
    float *corrVx = new float[Nl];
    float *corrRho = new float[Nl];
    float NxNy = Nx*Ny;
    for (int iz = 0; iz != Nz; iz++) {
        // Initilization
        CorrL[iz] = 0;      lcMx = Nl-1;    MxBar = 0;
        CorrL[iz+Nz] = 0;   lcVx = Nl-1;    VxBar = 0;
        CorrL[iz+2*Nz] = 0; lcRho = Nl-1;   RhoBar = 0;
        // Calculate Mean Value
        for (int iy = 0; iy != Ny; iy++) {
            for (int ix = 0; ix != Nx; ix++) {
                MxBar += cd_vector[0].data[iz][iy][ix][0];
                VxBar += cd_vector[0].data[iz][iy][ix][0]/cd_scalar[0].data[iz][iy][ix];
                RhoBar += cd_scalar[0].data[iz][iy][ix];
            }
        }
        // RL: the explicit method
        //MxBar /= NxNy;
        //VxBar /= NxNy;
        //RhoBar /= NxNy;
        
        // Calculate Correlation Value
        for (l = 0; l != Nl; l ++) {
            corrMx[l] = 0; MxBar = 0;
            corrVx[l] = 0; VxBar = 0;
            corrRho[l] = 0; RhoBar = 0;
            for (int iy = 0; iy != Ny; iy++) {
                for (int ix = 0; ix != Nx; ix++) {
                    // RL: the explicit method
                    //MxMx = ((cd_vector[0].data[iz][iy][ix][0]-MxBar) * (cd_vector[0].data[iz][iy][(ix+l)%Nx][0]-MxBar));
                    //RhoRho = ((cd_scalar[0].data[iz][iy][ix]-RhoBar) * (cd_scalar[0].data[iz][iy][(ix+l)%Nx]-RhoBar));
                    //corrVx[l] += ((cd_vector[0].data[iz][iy][ix][0]/cd_scalar[0].data[iz][iy][ix] - VxBar) * (cd_vector[0].data[iz][iy][(ix+l)%Nx][0]/cd_scalar[0].data[iz][iy][(ix+l)%Nx] - VxBar));
                    
                    // RL: the more effective way to do it
                    MxMx = (cd_vector[0].data[iz][iy][ix][0] * cd_vector[0].data[iz][iy][(ix+l)%Nx][0]);
                    corrMx[l] += MxMx;
                    RhoRho = (cd_scalar[0].data[iz][iy][ix] * cd_scalar[0].data[iz][iy][(ix+l)%Nx]);
                    corrRho[l] += RhoRho;
                    corrVx[l] += MxMx/RhoRho;
                }
            }
        }
        // Substract Mean Square (RL: more effective way)
        corrMx[l] -= (MxBar * MxBar / NxNy);
        corrVx[l] -= (VxBar * VxBar / NxNy);
        corrRho[l] -= (RhoBar * RhoBar / NxNy);
        
        for (l = Nl-1; l >= 0; l--) {
            corrMx[l] /= corrMx[0];
            corrVx[l] /= corrVx[0];
            corrRho[l] /= corrRho[0];
#ifdef CorrValue
            int Nlines = Nz * Nl;
            CorrV[Nl*iz+l] = corrMx[l];
            CorrV[Nlines+Nl*iz+l] = corrVx[l];
            CorrV[2*Nlines+Nl*iz+l] = corrRho[l];            
#endif
        }
        
        for (l = 1; l < Nl; l++) {
            // original criteria
            if (corrMx[l] < 0.5 && corrMx[l-1] > 0.5) {
                lcMx = (corrMx[l-1]-0.5)<(0.5-corrMx[l]) ? l-1:l;
                break;
            }
            //*/
            /* This criteria is trying to find the inflection point
            if ((corrMx[l-1]-corrMx[l]) > (corrMx[l]-corrMx[l+1])) {
                lcMx = l;
                break;
            }
            //*/
        }
        for (l = 1; l < Nl; l++) {
            if (corrRho[l] < 0.5 && corrRho[l-1] > 0.5) {
                lcRho = (corrRho[l-1]-0.5)<(0.5-corrRho[l]) ? l-1:l;
                break;
            }
        }
        for (l = 1; l < Nl; l++) {
            if (corrVx[l] < 0.5 && corrVx[l-1] > 0.5) {
                lcVx = (corrVx[l-1]-0.5)<(0.5-corrVx[l]) ? l-1:l;
                break;
            }
        }
        
        CorrL[iz] = lcMx * spacing[0];
        CorrL[iz+Nz] = lcVx * spacing[0];
        CorrL[iz+2*Nz] = lcRho * spacing[0];
    }
    
    delete [] corrMx;
    delete [] corrVx;
    delete [] corrRho;
    return 0;
}

/********** GasPar **********/
/*! \fn int GasPar()
 *  \brief calculate the basic dynamic info */
int VtkFile::GasPar()
{
    float temp_p_gas_i, temp_p_par_i;
    // initialization to zero since they only use +=
    for (int i = 0; i != 32; i++) {
        dynscal[i] = 0;
    }
    for (int i = 0; i != 16; i++) {
        dynscal2[i] = 0;
    }
    for (int i = 0; i != 4*dimensions[2]; i++) {
        GPME[i] = 0;
        GPMEPar[i] = 0;
    }
    
    
    // serious compuation for dynscal
    for (int iz = 0; iz != dimensions[2]; iz++) {
        for (int iy = 0; iy != dimensions[1]; iy++) {
            for (int ix = 0; ix != dimensions[0]; ix++) {
                //std::copy(&(cd_vector[0].data[iz][iy][ix][0]), &(cd_vector[0].data[iz][iy][ix][4]), &(dynscal[3]));
                for (int i = 0; i != 3; i++) {
                    temp_p_gas_i = cd_vector[0].data[iz][iy][ix][i];
                    dynscal[i] += temp_p_gas_i; // p_gas
                    GPME[i*dimensions[2]+iz] += temp_p_gas_i;
                    dynscal[i+4] += 0.5*temp_p_gas_i*temp_p_gas_i/cd_scalar[0].data[iz][iy][ix]; // Ek_gas
                    temp_p_par_i = cd_vector[1].data[iz][iy][ix][i];
                    GPMEPar[i*dimensions[2]+iz] += temp_p_par_i;
                    dynscal[i+16] += temp_p_par_i; // p_par
                    // notice rho_par might be ~0 (or smaller)
                    if (cd_scalar[1].data[iz][iy][ix] > 1e-15 || cd_scalar[1].data[iz][iy][ix] < -1e-15) {
                        dynscal[i+20] += 0.5*temp_p_par_i*temp_p_par_i/cd_scalar[1].data[iz][iy][ix]; // Ek_par
                    }
                }
            }
        }
        for (int i = 0; i != 3; i++) {
            GPME[i*dimensions[2]+iz] /= (dimensions[1]*dimensions[0]);
            GPMEPar[i*dimensions[2]+iz] /= (dimensions[1]*dimensions[0]);
        }
        GPME[3*dimensions[2]+iz] = sqrt(GPME[iz]*GPME[iz]+GPME[dimensions[2]+iz]*GPME[dimensions[2]+iz]+GPME[2*dimensions[2]+iz]*GPME[2*dimensions[2]+iz]);
        GPMEPar[3*dimensions[2]+iz] = sqrt(GPMEPar[iz]*GPMEPar[iz]+GPMEPar[dimensions[2]+iz]*GPMEPar[dimensions[2]+iz]+GPMEPar[2*dimensions[2]+iz]*GPMEPar[2*dimensions[2]+iz]);
    }
    
    /* Calculate sigma of
    float rho_g[256] = {0}, sigma_rho_g[256] = {0};
    for (int iz = 0; iz != dimensions[2]; iz++) {
        for (int iy = 0; iy != dimensions[1]; iy++) {
            for (int ix = 0; ix != dimensions[0]; ix++) {
                rho_g[iz] += cd_scalar[0].data[iz][iy][ix];
            }
        }
        rho_g[iz] /= (dimensions[1]*dimensions[0]);
        for (int iy = 0; iy != dimensions[1]; iy++) {
            for (int ix = 0; ix != dimensions[0]; ix++) {
                sigma_rho_g[iz] += (cd_scalar[0].data[iz][iy][ix]-rho_g[iz])*(cd_scalar[0].data[iz][iy][ix]-rho_g[iz]);
            }
        }
        sigma_rho_g[iz]/= (dimensions[1]*dimensions[0]);
        sigma_rho_g[iz] = sqrt(sigma_rho_g[iz]);
    }
    for (int i = 0; i != 256; i++) {
        cout << -0.4+0.003125*(i+0.5) << " " << sigma_rho_g[i] << endl;
    }
    exit(0);
     */ 
    
    dynscal[3] = sqrt(dynscal[0]*dynscal[0]+dynscal[1]*dynscal[1]+dynscal[2]*dynscal[2]); // p_gas_tot
    dynscal[7] = dynscal[4]+dynscal[5]+dynscal[6]; // Ek_gas_tot
    dynscal[19] = sqrt(dynscal[16]*dynscal[16]+dynscal[17]*dynscal[17]+dynscal[18]*dynscal[18]); // p_par_tot
    dynscal[23] = dynscal[20]+dynscal[21]+dynscal[22]; // Ek_par_tot
    
    // See global.h
    // volume-averaged: *1/Nx/Ny/Nz
    // area-averaged: *dz/Nx/Ny
    
    for (int i = 0; i != 8; i++) {
        dynscal[i+8] = dynscal[i] * spacing[2] / dimensions[1] / dimensions[0];
        dynscal[i] /= (dimensions[2]*dimensions[1]*dimensions[0]);
        dynscal[i+24] = dynscal[i+16] * spacing[2] / dimensions[1] / dimensions[0];
        dynscal[i+16] /= (dimensions[2]*dimensions[1]*dimensions[0]);
    }
    
    // serious compuation for dynscal2
    for (int iz = kps; iz != kpe; iz++) {
        for (int iy = 0; iy != dimensions[1]; iy++) {
            for (int ix = 0; ix != dimensions[0]; ix++) {
                //std::copy(&(cd_vector[0].data[iz][iy][ix][0]), &(cd_vector[0].data[iz][iy][ix][4]), &(dynscal[3]));
                for (int i = 0; i != 3; i++) {
                    temp_p_gas_i = cd_vector[0].data[iz][iy][ix][i];
                    dynscal2[i] += temp_p_gas_i; // p_gas
                    dynscal2[i+4] += 0.5*temp_p_gas_i*temp_p_gas_i/cd_scalar[0].data[iz][iy][ix]; // Ek_gas
                }
            }
        }
    }
    
    dynscal2[3] = sqrt(dynscal[0]*dynscal[0]+dynscal[1]*dynscal[1]+dynscal[2]*dynscal[2]); // p_gas_tot
    dynscal2[7] = dynscal[4]+dynscal[5]+dynscal[6]; // Ek_gas_tot
    
    // See global.h
    // volume-averaged: *1/Nx/Ny/Nz
    // area-averaged: *dz/Nx/Ny
    
    for (int i = 0; i != 8; i++) {
        dynscal2[i+8] = dynscal[i] * spacing[2] / dimensions[1] / dimensions[0];
        dynscal2[i] /= (dimensions[2]*dimensions[1]*dimensions[0]);
    }
    
    return 0;
}







