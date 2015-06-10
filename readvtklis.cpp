//
//  main.cpp
//  readvtklis
//
//  Created by Rixin Li on 1/14/15.
//  Copyright (c) 2015 Rixin Li. All rights reserved.
//

#include "global.h"
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
    std::ios_base::sync_with_stdio(false); // speed up c++
    clock_t begin_t, end_t;
    double elapsed_secs;
    begin_t = clock();
#ifdef ENABLE_MPI
    double mpi_begin_t, mpi_end_t;
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
    myMPI->paras.AllocateMemory(fio->n_file);
    myMPI->Barrier();
    // for debug
    cout << "Processor " << myMPI->myrank << ": " << myMPI->loop_begin << " " << myMPI->loop_end << " " << myMPI->loop_offset << "\n";
#endif
    // for reading V_gas_0 if VpecG_flag is set
    if (fio->VpecG_flag || fio->MeanSigma_flag || fio->VertRho_flag) {
        if (vf->Read_Header_Record_Pos(fio->vtk_filenames[0])) {
            cout << "Having problem reading header..." << "\n";
            exit(1);
        }
        vf->Read_Data(fio->vtk_filenames[0]);
        if (fio->VertRho_flag) {
            fio->paras.V_gas_0 = allocate3d_vector_array<float>(vf->dimensions);
            for (int i = 0; i != vf->dimensions[2]; i++) {
                for (int j = 0; j != vf->dimensions[1]; j++) {
                    for (int k = 0; k != vf->dimensions[0]; k++) {
                        for (int l = 0; l != 3; l++) {
                            fio->paras.V_gas_0[i][j][k][l] = vf->cd_vector[0].data[i][j][k][l]/vf->cd_scalar[0].data[i][j][k];;
                        }
                    }
                }
            }
        }
        if (fio->MeanSigma_flag) {
            vf->Sigma_gas_0_inbox = 0;
            for (int ix = 0; ix != vf->dimensions[0]; ix++) {
                double temp_sigma_gas_0_in_box = 0;
                for (int iy = 0; iy != vf->dimensions[1]; iy++) {
                    for (int iz = 0; iz != vf->dimensions[2]; iz++) {
                        temp_sigma_gas_0_in_box += vf->cd_scalar[0].data[iz][iy][ix];
                    }
                }
                vf->Sigma_gas_0_inbox += temp_sigma_gas_0_in_box;
            }
            vf->Sigma_gas_0_inbox /= vf->dimensions[0];
        }
        fio->paras.AllocateSubMemory(fio->n_file, vf->dimensions);
#ifdef ENABLE_MPI
        myMPI->paras.AllocateSubMemory(fio->n_file, vf->dimensions);
        myMPI->Barrier();
#endif
        
    }
#ifdef ENABLE_MPI
    for (int i = myMPI->loop_begin ; i <= myMPI->loop_end; i += myMPI->loop_offset) {
#else
        for (int i = 0; i < fio->n_file; i++) {
#endif
            // vtk part
            cout
#ifdef ENABLE_MPI
            << "Processor " << myMPI->myrank << ": "
#endif
            << "Reading " << fio->iof.data_basename+"." << setw(4) << setfill('0') << i+fio->start_no << "\n"; //" " << fio->ParNum_flag << fio->RhoParMax_flag << fio->HeiPar_flag << "\n";
            
            if (vf->Read_Header_Record_Pos(fio->vtk_filenames[i])) {
                cout << "Having problem reading header..." << "\n";
                exit(1);
            }
#ifdef ENABLE_MPI
            // some setup in first step
            // since some initialization happens inside the loop
            // we may require n_file should larger than numprocs
            if (i == myMPI->loop_begin) {
                ;
            }
#else
            //vf->Print_File_Info();
#endif
            if (fio->RhoParMax_flag || fio->MeanSigma_flag || fio->VpecG_flag || fio->VertRho_flag || fio->dSigma_flag || fio->CorrL_flag) {
                vf->Read_Data(fio->vtk_filenames[i]);
                vf->Calculate_Mass_Find_Max();
                //cout << "m_gas = " << vf->m_gas << "; m_par = " << vf->m_par << "\n";
                
                // I have checked the total gas mass and par mass, which is corresponding to mratio = 0.02
            }
            if (fio->ParNum_flag || fio->HeiPar_flag) {
                // lis part
                pl->ReadLis(fio->lis_filenames[i]);
            }
            
            
            // recording data
#ifdef ENABLE_MPI
            myMPI->paras.Otime[i] = vf->time;
            if (fio->ParNum_flag) {
                myMPI->paras.N_par[i] = pl->n;
            }
            if (fio->RhoParMax_flag) {
                myMPI->paras.Max_Rhop[i] = vf->Max_Rhop;
            }
            if (fio->HeiPar_flag) {
                pl->ScaleHeight(myMPI->paras.Hp[i], myMPI->paras.Hp_in1sigma[i]);
            }
            if (fio->dSigma_flag) {
                myMPI->paras.dSigma[i] = vf->dSigma;
            }
            if (fio->MeanSigma_flag) {
                vf->MeanSigma(myMPI->paras.MeanSigma[i]);
            }
            if (fio->VpecG_flag) {
                vf->VpecG(myMPI->paras.VpecG[i]);
            }
            if (fio->VertRho_flag) {
                vf->VertRho(myMPI->paras.VertRho[i]);
            }
            if (fio->CorrL_flag) {
                vf->CorrLen(myMPI->paras.CorrL[i]);
            }
            
#else
            fio->paras.Otime[i] = vf->time;
            if (fio->ParNum_flag) {
                fio->paras.N_par[i] = pl->n;
            }
            if (fio->RhoParMax_flag) {
                fio->paras.Max_Rhop[i] = vf->Max_Rhop;
            }
            if (fio->HeiPar_flag) {
                pl->ScaleHeight(fio->paras.Hp[i], fio->paras.Hp_in1sigma[i]);
            }
            if (fio->dSigma_flag) {
                fio->paras.dSigma[i] = vf->dSigma;
            }
            if (fio->MeanSigma_flag) {
                vf->MeanSigma(fio->paras.MeanSigma[i]);
            }
            if (fio->VpecG_flag) {
                vf->VpecG(fio->paras.VpecG[i]);
            }
            if (fio->VertRho_flag) {
                vf->VertRho(fio->paras.VertRho[i]);
            }
            if (fio->CorrL_flag) {
                vf->CorrLen(fio->paras.CorrL[i]);
            }
            
            // for temp output test
            /*
             ofstream file_ccpos;
             string file_ccpos_name = fio->iof.output_path_name;
             file_ccpos.open((file_ccpos_name.append(to_string(i))).c_str(), ofstream::out);
             if (!file_ccpos.is_open()) {
             cout << "Failed to open " << (char)i+fio->iof.output_path_name.c_str() << "\n";
             return 1;
             }
             int temp_z = 32;
             for (int j = 0; j != fio->paras.dimensions[1]; j++) {
             for (int k = 0; k != fio->paras.dimensions[0]; k++) {
             file_ccpos << setw(15) << scientific << vf->cd_vector[0].data[temp_z][j][k][2]/vf->cd_scalar[0].data[temp_z][j][k];
             }
             file_ccpos << "\n";
             }
             file_ccpos.close();
             //*/
            
            
            
#endif
#ifndef RESIZE_LIST
            pl->InitializeList();
#endif
        }
        
#ifdef ENABLE_MPI
        myMPI->Barrier();
        MPI::COMM_WORLD.Allreduce(myMPI->paras.Otime, fio->paras.Otime, fio->n_file, MPI::DOUBLE, MPI::SUM);
        if (fio->ParNum_flag) {
            MPI::COMM_WORLD.Allreduce(myMPI->paras.N_par, fio->paras.N_par, fio->n_file, MPI::LONG, MPI::SUM);
        }
        if (fio->RhoParMax_flag) {
            MPI::COMM_WORLD.Allreduce(myMPI->paras.Max_Rhop, fio->paras.Max_Rhop, fio->n_file, MPI::DOUBLE, MPI::SUM);
        }
        if (fio->HeiPar_flag) {
            MPI::COMM_WORLD.Allreduce(myMPI->paras.Hp, fio->paras.Hp, fio->n_file, MPI::DOUBLE, MPI::SUM);
            MPI::COMM_WORLD.Allreduce(myMPI->paras.Hp_in1sigma, fio->paras.Hp_in1sigma, fio->n_file, MPI::DOUBLE, MPI::SUM);
        }
        if (fio->dSigma_flag) {
            MPI::COMM_WORLD.Allreduce(myMPI->paras.dSigma, fio->paras.dSigma, fio->n_file, MPI::DOUBLE, MPI::SUM);
        }
        if (fio->MeanSigma_flag) {
            for (int i = 0; i != fio->n_file; i++) {
                MPI::COMM_WORLD.Allreduce(myMPI->paras.MeanSigma[i], fio->paras.MeanSigma[i], 2*vf->dimensions[0], MPI::DOUBLE, MPI::SUM);
            }
        }
        if (fio->VpecG_flag) {
            for (int i = 0; i != fio->n_file; i++) {
                MPI::COMM_WORLD.Allreduce(myMPI->paras.VpecG[i], fio->paras.VpecG[i], 3*vf->dimensions[2], MPI::DOUBLE, MPI::SUM);
            }
        }
        if (fio->VertRho_flag) {
            for (int i = 0; i != fio->n_file; i++) {
                MPI::COMM_WORLD.Allreduce(myMPI->paras.VertRho[i], fio->paras.VertRho[i], 2*vf->dimensions[2], MPI::DOUBLE, MPI::SUM);
            }
        }
        if (fio->CorrL_flag) {
            for (int i = 0; i != fio->n_file; i++) {
                MPI::COMM_WORLD.Allreduce(myMPI->paras.CorrL[i], fio->paras.CorrL, vf->dimensions[2], MPI::DOUBLE, MPI::SUM);
            }
        }
        cout << "Processor " << myMPI->myrank << ": I'm done." << "\n";
        myMPI->Barrier();
        if (myMPI->myrank == myMPI->master) {
#endif
            fio->Output_Data();
            /*
             for (int i = 0; i != fio->n_file; i++) {
             cout << "time = " << fio->Otime[i] << "; Max_Rhop = " << fio->Max_Rhop[i] << "; Hp = " << fio->Hp[i] << "\n";
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
            << "\n";
#ifdef ENABLE_MPI
        }
        myMPI->Finalize();
        
        delete myMPI;
        
#endif
        return 0;
    }
