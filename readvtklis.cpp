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
#include "octree.h"
#include <unistd.h>

using namespace std;

/********** global variables **********/
#ifdef ENABLE_MPI
MPI_info *myMPI = new MPI_info;
#endif
FileIO *fio = new FileIO;

/********** sub-routines **********/
int RecordData(int i, VtkFile *vf, ParticleList *pl, Octree * ot);
int IntegrateData(VtkFile *vf, Octree *ot);

/*! \fn int main(int argc, const char * argv[])
 *  \brief main function. */
int main(int argc, const char * argv[]) {
    std::ios_base::sync_with_stdio(false); // speed up c++
    
    /****************************Initialization*******************************/
#ifdef ENABLE_MPI
    double mpi_begin_t, mpi_end_t;
    myMPI->Initialize(argc, argv);
    mpi_begin_t = MPI::Wtime();
#else
    clock_t begin_t = clock();
#endif
    //sleep(15);
    fio->Initialize(argc, argv);
    fio->Generate_Filename();
    ParticleList *pl = new ParticleList;
    VtkFile *vf = new VtkFile;
    Octree *ot = new Octree;
    
#ifdef ENABLE_MPI
    if (myMPI->myrank == myMPI->master) {
#endif
        fio->Check_Input_Path_Filename();
        fio->Print_Stars("Master: Begin to Process Data");
#ifdef ENABLE_MPI
    }
    myMPI->Determine_Loop(fio->n_file);
    if (fio->RhopMaxPerLevel_flag) {
        myMPI->loop_end = -1;
    }
    myMPI->paras.AllocateMemory(fio->n_file);
    myMPI->Barrier();
    // for debug
    cout << "Processor " << myMPI->myrank << ": " << myMPI->loop_begin << " " << myMPI->loop_end << " " << myMPI->loop_offset << endl;
#endif
    
    // reading initial gas properties if needed
    if (fio->MeanSigma_flag) {
        // don't worry, letting multi-cpus read the same address is efficient enough
        if (vf->Read_Header_Record_Pos(fio->vtk_filenames[0])) {
            cout << "Having problem reading header..." << endl;
            exit(1);
        }
        vf->Read_Data(fio->vtk_filenames[0]);
        
        if (fio->MeanSigma_flag) {
            vf->Sigma_gas_0_inbox = 0;
            for (int ix = 0; ix != vf->dimensions[0]; ix++) {
                float temp_sigma_gas_0_in_box = 0;
                for (int iy = 0; iy != vf->dimensions[1]; iy++) {
                    for (int iz = 0; iz != vf->dimensions[2]; iz++) {
                        temp_sigma_gas_0_in_box += vf->cd_scalar[0].data[iz][iy][ix];
                    }
                }
                vf->Sigma_gas_0_inbox += temp_sigma_gas_0_in_box;
            }
            vf->Sigma_gas_0_inbox /= vf->dimensions[0];
        }
    }
    if (fio->MeanSigma_flag || fio->VpecG_flag || fio->VertRho_flag || fio->CorrL_flag) {
        fio->paras.AllocateSubMemory(fio->n_file, vf->dimensions);
#ifdef ENABLE_MPI
        myMPI->paras.AllocateSubMemory(fio->n_file, vf->dimensions);
        myMPI->Barrier();
#endif
    }
    
    /****************************Loop Begin*******************************/
    
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
            << "Reading " << fio->iof.data_basename+"." << setw(4) << setfill('0') << i*fio->interval+fio->start_no << endl;
            if (vf->Read_Header_Record_Pos(fio->vtk_filenames[i])) {
                cout << "Having problem reading header..." << endl;
                exit(1);
            }

            if (fio->RhoParMax_flag || fio->MeanSigma_flag || fio->VpecG_flag || fio->VertRho_flag || fio->dSigma_flag || fio->CorrL_flag) {
                vf->Read_Data(fio->vtk_filenames[i]);
                vf->Calculate_Mass_Find_Max();
                // I have checked the total gas mass and par mass, which is corresponding to mratio = 0.02
            }
            
            // lis part
            if (fio->RhoParMax_flag || fio->ParNum_flag || fio->HeiPar_flag) {
                pl->ReadLis(fio->lis_filenames[i]);
            }
            
            // Then something based on both vtk and lis
            if (fio->RhoParMax_flag) {
                ot->BuildTree(vf, pl);
            }
            
            // record_data
            RecordData(i, vf, pl, ot);
#ifndef RESIZE_LIST
            pl->InitializeList();
#endif
            
#ifdef ENABLE_MPI
        }
#else
    }
#endif
    
    if (fio->RhopMaxPerLevel_flag) {
#ifdef ENABLE_MPI
        for (int i = 0; i < fio->n_file; i++) {
            // read data
            if (vf->Read_Header_Record_Pos(fio->vtk_filenames[i])) {
                cout << "Having problem reading header..." << endl;
                exit(1);
            }
            vf->Read_Data(fio->vtk_filenames[i]);
            pl->ReadLis(fio->lis_filenames[i]);
            // build tree
            ot->BuildTree(vf, pl);
            if (fio->paras.RMPL == NULL) {
                fio->paras.RMPL = new float[ot->level+1];
                for (int i = 0; i <= ot->level; i++) {
                    fio->paras.RMPL[i] = 0;
                }
            }
            double temp_t = myMPI->T();
            myMPI->Barrier();
            myMPI->wait_time += myMPI->T() - temp_t;
            
            ot->RhopMaxPerLevel();
        }
        for (int i = 0; i <= ot->level; i++) {
            ot->RMPL[i] = ot->m1par * ot->RMPL[i] / (ot->foPio3 * ot->Radius[i] * ot->Radius[i] * ot->Radius[i]);
        }
        
#else /* ENABLE_MPI */
        cout << "Such a heavy computation will take forever in single CPU. Please consider using parallel computing." << endl;
        exit(1);
#endif /* ENABLE_MPI */
    }
    
    /************************** Data Output **************************/
    
#ifdef ENABLE_MPI
    myMPI->Barrier();
    IntegrateData(vf, ot);
    cout << "Processor " << myMPI->myrank << ": I'm done. (waiting " << myMPI->wait_time << "sec" << endl;
    myMPI->Barrier();
    
    if (myMPI->myrank == myMPI->master) {
#endif
        fio->Output_Data();
        /*
         for (int i = 0; i != fio->n_file; i++) {
         cout << "time = " << fio->Otime[i] << "; Max_Rhop = " << fio->Max_Rhop[i] << "; Hp = " << fio->Hp[i] << endl;
         }
         */
        fio->Print_Stars("Master: Finishing Program");
#ifdef ENABLE_MPI
    }
#endif
    
    /**************************** Timing ****************************/
#ifdef ENABLE_MPI
    mpi_end_t = MPI::Wtime();
    if (myMPI->myrank == myMPI->master) {
#endif
        cout << "Master: Elapsed time (secs) is "
#ifdef ENABLE_MPI
        << mpi_end_t - mpi_begin_t
#else
        << double(clock() - begin_t) / CLOCKS_PER_SEC
#endif
        << endl;
#ifdef ENABLE_MPI
    }
    myMPI->Finalize();
    //delete myMPI;
#endif
    // have bugs when I delete fio and pl, don't know why
    //delete fio;
    //delete pl;
    //delete vf;

    return 0;
}


/*! \fn int record_data(int i, VtkFile *vf, ParticleList *pl, Octree * ot)
 *  \brief record data of each loop into arrays that sotre output data */
int RecordData(int i, VtkFile *vf, ParticleList *pl, Octree * ot) {
    
    // recording data
#ifdef ENABLE_MPI
    myMPI->paras.Otime[i] = vf->time;
    if (fio->ParNum_flag) {
        myMPI->paras.N_par[i] = pl->n;
    }
    if (fio->RhoParMax_flag) {
        myMPI->paras.Max_Rhop[i] = ot->Max_Rhop;
        myMPI->paras.RpAV[i] = ot->RpAV;
        myMPI->paras.RpSQ[i] = ot->RpSQ;
        myMPI->paras.RpQU[i] = ot->RpQU;
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
        vf->CorrLen(myMPI->paras.CorrL[i]
#ifdef CorrValue
                    , myMPI->paras.CorrV[i]
#endif
                    );
    }
    
#else
    fio->paras.Otime[i] = vf->time;
    if (fio->ParNum_flag) {
        fio->paras.N_par[i] = pl->n;
    }
    if (fio->RhoParMax_flag) {
        fio->paras.Max_Rhop[i] = ot->Max_Rhop;
        fio->paras.RpAV[i] = ot->RpAV;
        fio->paras.RpSQ[i] = ot->RpSQ;
        fio->paras.RpQU[i] = ot->RpQU;
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
        vf->CorrLen(fio->paras.CorrL[i]
#ifdef CorrValue
                    , fio->paras.CorrV[i]
#endif
                    );
    }
    
#endif
    
    return 0;
}


/*! \fn int IntegrateData(VtkFile *vf, Octree *ot)
 *  \brief integrate data in each processor */
int IntegrateData(VtkFile *vf, Octree *ot)
{
    MPI::COMM_WORLD.Allreduce(myMPI->paras.Otime, fio->paras.Otime, fio->n_file, MPI::FLOAT, MPI::SUM);
    if (fio->ParNum_flag) {
        MPI::COMM_WORLD.Allreduce(myMPI->paras.N_par, fio->paras.N_par, fio->n_file, MPI::LONG, MPI::SUM);
    }
    if (fio->RhoParMax_flag) {
        MPI::COMM_WORLD.Allreduce(myMPI->paras.Max_Rhop, fio->paras.Max_Rhop, fio->n_file, MPI::FLOAT, MPI::SUM);
        MPI::COMM_WORLD.Allreduce(myMPI->paras.RpAV, fio->paras.RpAV, fio->n_file, MPI::FLOAT, MPI::SUM);
        MPI::COMM_WORLD.Allreduce(myMPI->paras.RpSQ, fio->paras.RpSQ, fio->n_file, MPI::FLOAT, MPI::SUM);
        MPI::COMM_WORLD.Allreduce(myMPI->paras.RpQU, fio->paras.RpQU, fio->n_file, MPI::FLOAT, MPI::SUM);
    }
    if (fio->HeiPar_flag) {
        MPI::COMM_WORLD.Allreduce(myMPI->paras.Hp, fio->paras.Hp, fio->n_file, MPI::FLOAT, MPI::SUM);
        MPI::COMM_WORLD.Allreduce(myMPI->paras.Hp_in1sigma, fio->paras.Hp_in1sigma, fio->n_file, MPI::FLOAT, MPI::SUM);
    }
    if (fio->dSigma_flag) {
        MPI::COMM_WORLD.Allreduce(myMPI->paras.dSigma, fio->paras.dSigma, fio->n_file, MPI::FLOAT, MPI::SUM);
    }
    if (fio->MeanSigma_flag) {
        for (int i = 0; i != fio->n_file; i++) {
            MPI::COMM_WORLD.Allreduce(myMPI->paras.MeanSigma[i], fio->paras.MeanSigma[i], 2*vf->dimensions[0], MPI::FLOAT, MPI::SUM);
        }
    }
    if (fio->VpecG_flag) {
        for (int i = 0; i != fio->n_file; i++) {
            MPI::COMM_WORLD.Allreduce(myMPI->paras.VpecG[i], fio->paras.VpecG[i], 3*vf->dimensions[2], MPI::FLOAT, MPI::SUM);
        }
    }
    if (fio->VertRho_flag) {
        for (int i = 0; i != fio->n_file; i++) {
            MPI::COMM_WORLD.Allreduce(myMPI->paras.VertRho[i], fio->paras.VertRho[i], 2*vf->dimensions[2], MPI::FLOAT, MPI::SUM);
        }
    }
    if (fio->CorrL_flag) {
        for (int i = 0; i != fio->n_file; i++) {
            MPI::COMM_WORLD.Allreduce(myMPI->paras.CorrL[i], fio->paras.CorrL[i], 3*vf->dimensions[2], MPI::FLOAT, MPI::SUM);
            
#ifdef CorrValue
            MPI::COMM_WORLD.Allreduce(myMPI->paras.CorrV[i], fio->paras.CorrV[i], 3*vf->dimensions[2]*(vf->dimensions[0]/2+1), MPI::FLOAT, MPI::SUM);
#endif
        }
    }
    if (fio->RhopMaxPerLevel_flag) {
        MPI::COMM_WORLD.Allreduce(ot->RMPL, fio->paras.RMPL, ot->level+1, MPI::FLOAT, MPI::MAX);
    }

    return 0;
}




