//
//  main.cpp
//  readvtklis
//
//  Created by Rixin Li on 1/14/15.
//  Copyright (c) 2015 Rixin Li. All rights reserved.
//

#include "fop.h"
#include "readlis.h"
#include "readvtk.h"

using namespace std;

// global variables
#ifdef ENABLE_MPI
MPI_info *myMPI = new MPI_info;
#endif
FileIO *fio = new FileIO;

int main(int argc, const char * argv[]) {
    
    clock_t begin_t, end_t;
    double elapsed_secs, mpi_begin_t, mpi_end_t;
    begin_t = clock();
#ifdef ENABLE_MPI
    myMPI->Initialize(argc, argv);
    mpi_begin_t = MPI::Wtime();
#endif
    fio->Initialize(argc, argv);
    fio->Generate_Filename();
    ParticleList *pl = new ParticleList;
    VtkFile *vf = new VtkFile;

#ifdef ENABLE_MPI
    if (myMPI->myrank == myMPI->master) {
#endif
        fio->Check_Input_Path_Filename();
        fio->Print_Stars("Master: Begin to Process Data");
#ifdef ENABLE_MPI
    }
    myMPI->Determine_Loop(fio->n_file);
    myMPI->NewVars(fio->n_file, fio->n_cpu);
    myMPI->Barrier();
    // for debug
    cout << "Processor " << myMPI->myrank << ": " << myMPI->loop_begin << " " << myMPI->loop_end << " " << myMPI->loop_offset << endl;
    
    // for reading V_gas_0
    if (fio->VturbGas_flag) {
        
        //if (myMPI->myrank == myMPI->master) {
            if (vf->Read_Header_Record_Pos(fio->vtk_filenames[0])) {
                cout << "Having problem reading header..." << endl;
                exit(1);
            }
            vf->Read_Data(fio->vtk_filenames[0]);
            fio->V_gas_0 = new float***[vf->dimensions[2]];
            for (int i = 0; i != vf->dimensions[2]; i++) {
                fio->V_gas_0[i] = new float**[vf->dimensions[1]];
                for (int j = 0; j != vf->dimensions[1]; j++) {
                    fio->V_gas_0[i][j] = new float*[vf->dimensions[0]];
                    for (int k = 0; k != vf->dimensions[0]; k++) {
                        fio->V_gas_0[i][j][k] = new float[3];
                        for (int l = 0; l != 3; l++) {
                            fio->V_gas_0[i][j][k][l] = vf->cd_vector[0].data[i][j][k][l]/vf->cd_scalar[0].data[i][j][k];;
                        }
                    }
                }
            }
        //} else {
            
        //}
        myMPI->Barrier();
    }
    
    for (int i = myMPI->loop_begin ; i <= myMPI->loop_end; i += myMPI->loop_offset) {
#else
    for (int i = 0; i < fio->n_file; i++) {
#endif
        // vtk part
        cout
#ifdef ENABLE_MPI
        << "Processor " << myMPI->myrank << ": "
#endif
        << "Reading " << fio->data_basename+"." << setw(4) << setfill('0') << i+fio->start_no << endl; //" " << fio->ParNum_flag << fio->RhoParMax_flag << fio->HeiPar_flag << endl;
        
        if (vf->Read_Header_Record_Pos(fio->vtk_filenames[i])) {
            cout << "Having problem reading header..." << endl;
            exit(1);
        }
#ifdef ENABLE_MPI
        // some setup in first step
        // since some initialization happens inside the loop
        // we require n_file should larger than numprocs
        if (i == myMPI->loop_begin) {
            if (fio->SigmaParY_flag) {
                for (int j = 0; j != fio->n_file; j++) {
                    myMPI->Sigma_par_y[j] = new double[vf->dimensions[0]];
                    for (int k = 0; k != vf->dimensions[0]; k++) {
                        myMPI->Sigma_par_y[j][k] = 0;
                    }
                }
            }
            if (fio->VturbGas_flag) {
                for (int j = 0; j != fio->n_file; j++) {
                    myMPI->Vturb_gas[j] = new double[vf->dimensions[2]];
                    for (int k = 0; k != vf->dimensions[2]; k++) {
                        myMPI->Vturb_gas[j][k] = 0;
                    }
                }
            }
            
        }
#else
        //vf->Print_File_Info();
#endif
        if (fio->RhoParMax_flag || fio->SigmaParY_flag || fio->VturbGas_flag) {
            vf->Read_Data(fio->vtk_filenames[i]);
            vf->Calculate_Mass_Find_Max();
            //cout << "m_gas = " << vf->m_gas << "; m_par = " << vf->m_par << endl;
            // I have checked the total gas mass and par mass, which is corresponding to mratio = 0.02
        }
        if (fio->ParNum_flag || fio->HeiPar_flag || fio->CpuID_flag) {
            // lis part
            pl->ReadLis(fio->lis_filenames[i]);
        }
                
        // recording data
#ifdef ENABLE_MPI
        myMPI->orbit_time[i] = vf->time;
        if (fio->ParNum_flag) {
            myMPI->n_par[i] = pl->n;
        }
        if (fio->RhoParMax_flag) {
            myMPI->max_rho_par[i] = vf->max_rho_par;
        }
        if (fio->HeiPar_flag) {
            pl->ScaleHeight(myMPI->Hp[i], myMPI->Hp_in1sigma[i]);
        }
        if (fio->CpuID_flag) {
            pl->CpuID();
            delete [] myMPI->cpuid_dist[i];
            myMPI->cpuid_dist[i] = pl->CpuID_dist;
        }
        if (fio->SigmaParY_flag) {
            for (int ix = 0; ix != vf->dimensions[0]; ix++) {
                for (int iz = 0; iz != vf->dimensions[2]; iz++) {
                    for (int iy = 0; iy != vf->dimensions[1]; iy++) {
                        myMPI->Sigma_par_y[i][ix] += vf->cd_scalar[1].data[iz][iy][ix] * vf->cell_volume;
                    }
                }
                myMPI->Sigma_par_y[i][ix] /= (vf->Sigma_gas_0 * vf->spacing[0] * vf->spacing[1] * vf->dimensions[1]);
            }
        }
        if (fio->VturbGas_flag) {
            vf->Calculate_Vturb_Gas(myMPI->Vturb_gas[i], fio->V_gas_0);
        }
        
#else
        fio->orbit_time[i] = vf->time;
        if (fio->ParNum_flag) {
            fio->n_par[i] = pl->n;
        }
        if (fio->RhoParMax_flag) {
            fio->max_rho_par[i] = vf->max_rho_par;
        }
        if (fio->HeiPar_flag) {
            pl->ScaleHeight(fio->Hp[i], fio->Hp_in1sigma[i]);
        }
        if (fio->CpuID_flag) {
            pl->CpuID();
            fio->CpuID_dist[i] = pl->CpuID_dist;
        }
        if (fio->SigmaParY_flag) {
            fio->Sigma_par_y[i] = new double[vf->dimensions[0]];
            for (int ix = 0; ix != vf->dimensions[0]; ix++) {
                fio->Sigma_par_y[i][ix] = 0;
                for (int iz = 0; iz != vf->dimensions[2]; iz++) {
                    for (int iy = 0; iy != vf->dimensions[1]; iy++) {
                        fio->Sigma_par_y[i][ix] += vf->cd_scalar[1].data[iz][iy][ix] * vf->cell_volume;
                    }
                }
                fio->Sigma_par_y[i][ix] /= (vf->Sigma_gas_0 * vf->spacing[0] * vf->spacing[1] * vf->dimensions[1]);
            }
        }
        
        // for temp output test
        /*
        ofstream file_ccpos;
        string file_ccpos_name = fio->output_path_name;
        file_ccpos.open((file_ccpos_name.append(to_string(i))).c_str(), ofstream::out);
        if (!file_ccpos.is_open()) {
            cout << "Failed to open " << (char)i+fio->output_path_name.c_str() << endl;
            return 1;
        }
        int temp_z = 32;
        for (int j = 0; j != fio->dimensions[1]; j++) {
            for (int k = 0; k != fio->dimensions[0]; k++) {
                file_ccpos << setw(15) << scientific << vf->cd_vector[0].data[temp_z][k][j][2]/vf->cd_scalar[0].data[k][j][temp_z];
            }
            file_ccpos << endl;
        }
        file_ccpos.close();
         */
        if (fio->VturbGas_flag) {
            if (i == 0) {
                if (fio->VturbGas_flag) {
                    fio->V_gas_0 = new float***[vf->dimensions[2]];
                    for (int i = 0; i != vf->dimensions[2]; i++) {
                        fio->V_gas_0[i] = new float**[vf->dimensions[1]];
                        for (int j = 0; j != vf->dimensions[1]; j++) {
                            fio->V_gas_0[i][j] = new float*[vf->dimensions[0]];
                            for (int k = 0; k != vf->dimensions[0]; k++) {
                                fio->V_gas_0[i][j][k] = new float[3];
                                for (int l = 0; l != 3; l++) {
                                    fio->V_gas_0[i][j][k][l] = vf->cd_vector[0].data[i][j][k][l]/vf->cd_scalar[0].data[i][j][k];
                                }
                            }
                        }
                    }
                }
            }
            fio->Vturb_gas[i] = new double[vf->dimensions[2]];
            vf->Calculate_Vturb_Gas(fio->Vturb_gas[i], fio->V_gas_0);
        }

#endif
#ifndef RESIZE_LIST
        pl->InitializeList();
#endif
    }
    
#ifdef ENABLE_MPI
    myMPI->Barrier();
    MPI::COMM_WORLD.Allreduce(myMPI->orbit_time, fio->orbit_time, fio->n_file, MPI::DOUBLE, MPI::SUM);
    if (fio->ParNum_flag) {
        MPI::COMM_WORLD.Allreduce(myMPI->n_par, fio->n_par, fio->n_file, MPI::LONG, MPI::SUM);
    }
    if (fio->RhoParMax_flag) {
        MPI::COMM_WORLD.Allreduce(myMPI->max_rho_par, fio->max_rho_par, fio->n_file, MPI::DOUBLE, MPI::SUM);
    }
    if (fio->HeiPar_flag) {
        MPI::COMM_WORLD.Allreduce(myMPI->Hp, fio->Hp, fio->n_file, MPI::DOUBLE, MPI::SUM);
        MPI::COMM_WORLD.Allreduce(myMPI->Hp_in1sigma, fio->Hp_in1sigma, fio->n_file, MPI::DOUBLE, MPI::SUM);
    }
    if (fio->CpuID_flag) {
        for (int i = 0; i != fio->n_file; i++) {
            fio->CpuID_dist[i] = new long[fio->n_cpu];
            MPI::COMM_WORLD.Allreduce(myMPI->cpuid_dist[i], fio->CpuID_dist[i], fio->n_cpu, MPI::LONG, MPI::SUM);
        }
    }
    if (fio->SigmaParY_flag) {
        for (int i = 0; i != fio->n_file; i++) {
            fio->Sigma_par_y[i] = new double[vf->dimensions[0]];
            MPI::COMM_WORLD.Allreduce(myMPI->Sigma_par_y[i], fio->Sigma_par_y[i], vf->dimensions[0], MPI::DOUBLE, MPI::SUM);
        }
    }
    if (fio->VturbGas_flag) {
        for (int i = 0; i != fio->n_file; i++) {
            fio->Vturb_gas[i] = new double[vf->dimensions[2]];
            MPI::COMM_WORLD.Allreduce(myMPI->Vturb_gas[i], fio->Vturb_gas[i], vf->dimensions[2], MPI::DOUBLE, MPI::SUM);
        }
    }
    cout << "Processor " << myMPI->myrank << ": I'm done." << endl;
    myMPI->Barrier();
    if (myMPI->myrank == myMPI->master) {
#endif
        fio->Output_Data();
        /*
        for (int i = 0; i != fio->n_file; i++) {
            cout << "time = " << fio->orbit_time[i] << "; max_rho_par = " << fio->max_rho_par[i] << "; Hp = " << fio->Hp[i] << endl;
        }
         */
        fio->Print_Stars("Master: Finishing Program");
#ifdef ENABLE_MPI
    }
#endif
    
    // have bugs when I delete fio and pl, don't know why
    delete fio;
    delete pl;
    delete vf;
    
    end_t = clock();
    elapsed_secs = double(end_t - begin_t) / CLOCKS_PER_SEC;
#ifdef ENABLE_MPI
    mpi_end_t = MPI::Wtime();
    if (myMPI->myrank == myMPI->master) {
#endif
        cout << "Master: Elapsed time (secs) is "
#ifdef ENABLE_MPI
        << mpi_end_t - mpi_begin_t
#else
        << elapsed_secs
#endif
        << endl;
#ifdef ENABLE_MPI
    }
    myMPI->Finalize();

    delete myMPI;
    
#endif
    return 0;
}
