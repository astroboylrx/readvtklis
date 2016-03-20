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

int GasParDynamic();
int wangbadan(int exitcode);

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
    Octree *ot = new Octree(fio->n_file);
    Quadtree *qt = new Quadtree(fio->n_file);
    
    myMPI->Barrier();
    // for debug
    cout << "Processor " << myMPI->myrank << ": " << myMPI->loop_begin << " " << myMPI->loop_end << " " << myMPI->loop_offset << endl;
#endif
    
    // reading initial gas properties if needed since the dimensions are needed
    if (fio->MeanSigma_flag || fio->VpecG_flag || fio->VertRho_flag || fio->CorrL_flag || fio->GasPar_flag) {
        // don't worry, letting multi-cpus read the same address is efficient enough for small tasks
        if (vf->Read_Header_Record_Pos(fio->vtk_filenames[0])) {
            cout << "Having problem reading header..." << endl;
            exit(1);
        }
        vf->Read_Data(fio->vtk_filenames[0]);
        fio->paras.AllocateSubMemory(fio->n_file, vf->dimensions);
#ifdef ENABLE_MPI
        myMPI->paras.AllocateSubMemory(fio->n_file, vf->dimensions);
        myMPI->Barrier();
#endif
        // MeanSigma need the initial average gas surface density
        if (fio->MeanSigma_flag) {
            vf->Sigma_gas_0_inbox = 0;
            vf->Sigma_gas_0_mid = 0;
            vf->Sigma_gas_0_in2Hp = 0;
            int mid[2] = {vf->dimensions[2]/2-1, vf->dimensions[2]/2};
            int Hp[2] = {vf->dimensions[2]/2-8, vf->dimensions[2]/2+7};
            
            for (int ix = 0; ix != vf->dimensions[0]; ix++) {
                double temp_Sigma_gas_0 = 0;
                for (int iy = 0; iy != vf->dimensions[1]; iy++) {
                    for (int iz = 0; iz != vf->dimensions[2]; iz++) {
                        // letting vf->Sigma_gas_0_inbox to do calculations will induce large error (at the forth effective digit, which will affect the final results)
                        // since there are Nx*Ny*Nz times of adding, so float cannot handle it very well
                        //vf->Sigma_gas_0_inbox += vf->cd_scalar[0].data[iz][iy][ix];
                        temp_Sigma_gas_0 += vf->cd_scalar[0].data[iz][iy][ix];
                    }
                    vf->Sigma_gas_0_mid += (vf->cd_scalar[0].data[mid[0]][iy][ix]+vf->cd_scalar[0].data[mid[1]][iy][ix]);
                    for (int iz = Hp[0]; iz != Hp[1]; iz++) {
                        vf->Sigma_gas_0_in2Hp += vf->cd_scalar[0].data[iz][iy][ix];
                    }
                }
                vf->Sigma_gas_0_inbox += temp_Sigma_gas_0;
            }
            vf->Sigma_gas_0_inbox /= vf->dimensions[0];
            vf->Sigma_gas_0_mid /= vf->dimensions[0];
            vf->Sigma_gas_0_in2Hp /= vf->dimensions[0];            
        }
    }

    /****************************Loop Begin*******************************/
    
#ifdef ENABLE_MPI
    for (int i = myMPI->loop_begin ; i <= myMPI->loop_end; i += myMPI->loop_offset) {
        cout << myMPI->prank()
#else
        for (int i = 0; i < fio->n_file; i++) {
            cout
#endif
            << "Reading " << fio->iof.data_basename+"." << setw(4) << setfill('0') << i*fio->interval+fio->start_no << endl;
            // vtk part
            if (vf->Read_Header_Record_Pos(fio->vtk_filenames[i])) {
                cout << "Having problem reading header..." << endl;
                exit(1);
            }
            
            if (fio->RhoParMax_flag || fio->MeanSigma_flag || fio->VpecG_flag || fio->VertRho_flag || fio->dSigma_flag || fio->CorrL_flag || fio->PointCloud_flag || fio->GasPar_flag) {
                vf->Read_Data(fio->vtk_filenames[i]);
                vf->Calculate_Mass_Find_Max();
                // I have checked the total gas mass and par mass, which is corresponding to mratio = 0.02
            }
            
            // lis part
            if (fio->RhoParMax_flag || fio->ParNum_flag || fio->HeiPar_flag || fio->PointCloud_flag || fio->GasPar_flag) {
                pl->ReadLis(fio->lis_filenames[i]);
            }
            
            // Then something based on both vtk and lis
            if (fio->RhoParMax_flag) {
                ot->BuildTree(vf, pl, i);
            }
            if (fio->PointCloud_flag) {
                pl->Lis2Vtk(fio->lis2vtk_filenames[i], vf->Header);
            }
            
            // calculate and record_data
            RecordData(i, vf, pl, ot);
#ifndef RESIZE_LIST
            pl->InitializeList();
#endif
            
#ifdef ENABLE_MPI
        }
#else
    }
#endif
    /* produce <rho_g>_xy(z)
    float rhog[64], sigma_rhog[64];
    for (int iz = 0; iz != vf->dimensions[2]; iz++) {
        rhog[iz] = 0;
        for (int iy = 0; iy != vf->dimensions[1]; iy++) {
            for (int ix = 0; ix != vf->dimensions[0]; ix++) {
                rhog[iz] += vf->cd_scalar[0].data[iz][iy][ix];
            }
        }
        rhog[iz] /= vf->dimensions[0]*vf->dimensions[1];
        sigma_rhog[iz] = 0;
        for (int iy = 0; iy != vf->dimensions[1]; iy++) {
            for (int ix = 0; ix != vf->dimensions[0]; ix++) {
                sigma_rhog[iz] += (vf->cd_scalar[0].data[iz][iy][ix] - rhog[iz])*(vf->cd_scalar[0].data[iz][iy][ix] - rhog[iz]);
            }
        }
        sigma_rhog[iz] /= vf->dimensions[0]*vf->dimensions[1];
        sigma_rhog[iz] = sqrt(sigma_rhog[iz]);
        cout << vf->cell_center[iz][0][0][2] << " " << rhog[iz] << " " << sigma_rhog[iz] << endl;
    } // */
    
    if (fio->RhopMaxPerLevel_flag) {
#ifdef ENABLE_MPI
        float *temp_RMPL = NULL;
        for (int i = 0; i < fio->n_file; i++) {
            // read data
            if (myMPI->myrank == 0) {
                cout << myMPI->prank() << "Reading " << fio->vtk_filenames[i] << endl;
            }
            if (vf->Read_Header_Record_Pos(fio->vtk_filenames[i])) {
                cout << "Having problem reading header..." << endl;
                exit(1);
            }
            vf->Read_Data(fio->vtk_filenames[i]);
            pl->ReadLis(fio->lis_filenames[i]);
#ifdef OCTREE
            // build tree
            ot->BuildTree(vf, pl, i);
            if (fio->paras.RMPL == NULL) {
                fio->paras.RMPL = new float[ot->level+1];
                for (int i = 0; i <= ot->level; i++) {
                    fio->paras.RMPL[i] = 0;
                }
            }
            if (temp_RMPL == NULL) {
                temp_RMPL = new float[ot->level+1];
                for (int i = 0; i <= ot->level; i++) {
                    temp_RMPL[i] = 0;
                }
            }
            double temp_t = myMPI->T();
            myMPI->Barrier();
            myMPI->wait_time += myMPI->T() - temp_t;
            
            ot->RhopMaxPerLevel(i);
            myMPI->Barrier();
            MPI::COMM_WORLD.Allreduce(ot->RMPL[i], temp_RMPL, ot->level+1, MPI::FLOAT, MPI::MAX);
            for (int level = 0; level <= ot->level; level++) {
                fio->paras.RMPL[level] += temp_RMPL[level];
            }
            
#endif /* OCTREE */
#ifdef QUADTREE
            // build tree
            qt->BuildTree(vf, pl, i);
            if (fio->paras.RMPL == NULL) {
                fio->paras.RMPL = new float[qt->level+1];
                for (int i = 0; i <= qt->level; i++) {
                    fio->paras.RMPL[i] = 0;
                }
            }
            if (temp_RMPL == NULL) {
                temp_RMPL = new float[qt->level+1];
                for (int i = 0; i <= qt->level; i++) {
                    temp_RMPL[i] = 0;
                }
            }
            double temp_t = myMPI->T();
            myMPI->Barrier();
            myMPI->wait_time += myMPI->T() - temp_t;
            
            qt->SigmapMaxPerLevel(i);
            myMPI->Barrier();
            MPI::COMM_WORLD.Allreduce(qt->SMPL[i], temp_RMPL, qt->level+1, MPI::FLOAT, MPI::MAX);
            for (int level = 0; level <= qt->level; level++) {
                fio->paras.RMPL[level] += temp_RMPL[level];
            }
#endif /* QUADTREE */
        }
        
#ifdef OCTREE
        for (int level = 0; level <= ot->level; level++) {
            fio->paras.RMPL[level] = ot->m1par * fio->paras.RMPL[level] / (ot->foPio3 * ot->Radius[level] * ot->Radius[level] * ot->Radius[level]) / fio->n_file;
        }
#endif /* OCTREE */
#ifdef QUADTREE
        for (int level = 0; level <= qt->level; level++) {
            fio->paras.RMPL[level] = qt->m1par * fio->paras.RMPL[level] / (PI * qt->Radius[level] * qt->Radius[level]) / fio->n_file;
        }
#endif /* QUADTREE */

        
#else /* ENABLE_MPI */
        cout << "Such a heavy computation will take forever in single CPU. Please consider using parallel computing." << endl;
        exit(1);
#endif /* ENABLE_MPI */
    }
    
#ifdef ENABLE_MPI
    myMPI->Barrier();
    IntegrateData(vf, ot);
    cout << myMPI->prank() << "I'm done. (waiting " << myMPI->wait_time << "sec" << endl;
    myMPI->Barrier();
#endif
    /************************** Multi-Read **************************/
    // deal with the calculations that need read data many times
    if (fio->GasPar_flag) {
        //GasParDynamic();
    }
    
    /************************** Data Output **************************/
#ifdef ENABLE_MPI
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
    Paras2probe *paras;
#ifdef ENABLE_MPI
    paras = &(myMPI->paras);
#else
    paras = &(fio->paras);
#endif
    paras->Otime[i] = vf->time;
    if (fio->ParNum_flag) {
        paras->N_par[i] = pl->n;
    }
    if (fio->RhoParMax_flag) {
        paras->Max_Rhop[i] = vf->Max_Rhop;
        paras->RpAV[i] = vf->RpAV;
        paras->RpSQ[i] = vf->RpSQ;
        paras->RpQU[i] = vf->RpQU;
    }
    if (fio->HeiPar_flag) {
        pl->ScaleHeight(paras->Hp[i], paras->Hp_in1sigma[i]);
    }
    if (fio->dSigma_flag) {
        paras->dSigma[i] = vf->dSigma;
    }
    if (fio->MeanSigma_flag) {
        vf->MeanSigma(paras->MeanSigma[i]);
        if (i == fio->n_file-1) {
            vf->Sigma(paras->Sigma);
        }
    }
    if (fio->VpecG_flag) {
        vf->VpecG(paras->VpecG[i]);
    }
    if (fio->VertRho_flag) {
        vf->VertRho(paras->VertRho[i]);
    }
    if (fio->CorrL_flag) {
        vf->CorrLen(paras->CorrL[i]
#ifdef CorrValue
                    , paras->CorrV[i]
#endif
                    );
    }
    if (fio->GasPar_flag) {
        vf->GasPar();
        pl->GasPar(vf->dimensions, vf->spacing);
        std::copy(&(vf->dynscal[0]), &(vf->dynscal[16]), &(paras->GasHst[i][0]));
        std::copy(&(vf->dynscal2[0]), &(vf->dynscal2[16]), &(paras->GasHst2[i][0]));
        std::copy(&(vf->dynscal[16]), &(vf->dynscal[32]), &(paras->ParHst[i][0]));
        std::copy(&(pl->dynscal[0]), &(pl->dynscal[20]), &(paras->ParLis[i][0]));
        std::copy(&(vf->GPME[0]), &(vf->GPME[4*vf->dimensions[2]]), &(paras->GPME[i][0]));
        std::copy(&(vf->GPMEPar[0]), &(vf->GPMEPar[4*vf->dimensions[2]]), &(paras->GPMEPar[i][0]));
    }
 
    return 0;
}

#ifdef ENABLE_MPI
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
        MPI::COMM_WORLD.Allreduce(myMPI->paras.RpAV, fio->paras.RpAV, fio->n_file, MPI::DOUBLE, MPI::SUM);
        MPI::COMM_WORLD.Allreduce(myMPI->paras.RpSQ, fio->paras.RpSQ, fio->n_file, MPI::DOUBLE, MPI::SUM);
        MPI::COMM_WORLD.Allreduce(myMPI->paras.RpQU, fio->paras.RpQU, fio->n_file, MPI::DOUBLE, MPI::SUM);
    }
    if (fio->HeiPar_flag) {
        MPI::COMM_WORLD.Allreduce(myMPI->paras.Hp, fio->paras.Hp, fio->n_file, MPI::DOUBLE, MPI::SUM);
        MPI::COMM_WORLD.Allreduce(myMPI->paras.Hp_in1sigma, fio->paras.Hp_in1sigma, fio->n_file, MPI::FLOAT, MPI::SUM);
    }
    if (fio->dSigma_flag) {
        MPI::COMM_WORLD.Allreduce(myMPI->paras.dSigma, fio->paras.dSigma, fio->n_file, MPI::DOUBLE, MPI::SUM);
    }
    if (fio->MeanSigma_flag) {
        for (int i = 0; i != fio->n_file; i++) {
            MPI::COMM_WORLD.Allreduce(myMPI->paras.MeanSigma[i], fio->paras.MeanSigma[i], 4*vf->dimensions[0], MPI::DOUBLE, MPI::SUM);
        }
        for (int i = 0; i != 4*vf->dimensions[0]; i++) {
            MPI::COMM_WORLD.Allreduce(myMPI->paras.Sigma[i], fio->paras.Sigma[i], vf->dimensions[1], MPI::DOUBLE, MPI::SUM);
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
            MPI::COMM_WORLD.Allreduce(myMPI->paras.CorrL[i], fio->paras.CorrL[i], 3*vf->dimensions[2], MPI::FLOAT, MPI::SUM);
            
#ifdef CorrValue
            MPI::COMM_WORLD.Allreduce(myMPI->paras.CorrV[i], fio->paras.CorrV[i], 3*vf->dimensions[2]*(vf->dimensions[0]/2+1), MPI::FLOAT, MPI::SUM);
#endif
        }
    }
#ifdef OCTREE
    //if (fio->RhopMaxPerLevel_flag) {
    //    MPI::COMM_WORLD.Allreduce(ot->RMPL[0], fio->paras.RMPL, ot->level+1, MPI::FLOAT, MPI::MAX);
    //}
#endif
    if (fio->GasPar_flag) {
        for (int i = 0; i != fio->n_file; i++) {
            MPI::COMM_WORLD.Allreduce(myMPI->paras.GasHst[i], fio->paras.GasHst[i], 16, MPI::DOUBLE, MPI::SUM);
            MPI::COMM_WORLD.Allreduce(myMPI->paras.GasHst2[i], fio->paras.GasHst2[i], 16, MPI::DOUBLE, MPI::SUM);
            MPI::COMM_WORLD.Allreduce(myMPI->paras.ParHst[i], fio->paras.ParHst[i], 16, MPI::DOUBLE, MPI::SUM);
            MPI::COMM_WORLD.Allreduce(myMPI->paras.ParLis[i], fio->paras.ParLis[i], 20, MPI::DOUBLE, MPI::SUM);
            MPI::COMM_WORLD.Allreduce(myMPI->paras.GPME[i], fio->paras.GPME[i], 4*(vf->dimensions[2]+1)*9, MPI::DOUBLE, MPI::SUM);
            MPI::COMM_WORLD.Allreduce(myMPI->paras.GPMEPar[i], fio->paras.GPMEPar[i], 4*(vf->dimensions[2]+1)*9, MPI::DOUBLE, MPI::SUM);
        }
    }
    
    return 0;
}

#endif

/*! \fn int GasParDynamic()
 *  \brief calculate dynamic info of gas/par */
int GasParDynamic()
{
    VtkFile *vf = new VtkFile;
    Paras2probe *paras = &fio->paras;
    int Nz = fio->paras.dimensions[2], Nz1 = Nz+1;
    int fourNz = 4*Nz, fourNz1=fourNz+4;
    double toe = 1e-15; // tolerance of error for rho_p
    
#ifdef ENABLE_MPI
    paras = &myMPI->paras;
#endif
    
    /************************** GPME **************************/
    // !!!!!!!!!!! first, copy data to temporary place and caluclate M_0 (at each processor)
    
    double **M_0 = new double*[fio->n_file];
    double **MP_0 = new double*[fio->n_file];
    for (int i = 0; i != fio->n_file; i++) {
        M_0[i] = new double[fourNz];
        MP_0[i] = new double[fourNz];
        for (int j = 0; j != fourNz; j++) {
            M_0[i][j] = fio->paras.GPME[i][j];
            MP_0[i][j] = fio->paras.GPMEPar[i][j];
            fio->paras.GPME[i][j] = 0;
            fio->paras.GPMEPar[i][j] = 0;
#ifdef ENABLE_MPI
            // we need read again and addreduce again
            myMPI->paras.GPME[i][j] = 0;
            myMPI->paras.GPMEPar[i][j] = 0;
#endif
        }
    }
    
    int Norb = 0, Ncount = 0, overhead, delimiter;
    int jd[4] = {0, Nz, 2*Nz, 3*Nz}; // j*dimension[2]
    int jd1[4] = {0, Nz+1, 2*(Nz+1), 3*(Nz+1)}; // j*(dim[2]+1)
    
    for (int i = 0; i != fio->n_file; i++) {
        Norb = int(floor(fio->paras.Otime[i]/(2*3.141592653)));
        Ncount++;
        for (int iz = 0; iz != 3*Nz; iz++) {
            overhead = 1 + iz / Nz; // leave space for entire-volume-averaged M_0
            fio->paras.GPME[0][iz+overhead] += M_0[i][iz];
            fio->paras.GPME[Norb+1][iz+overhead] += M_0[i][iz];
            fio->paras.GPMEPar[0][iz+overhead] += MP_0[i][iz];
            fio->paras.GPMEPar[Norb+1][iz+overhead] += MP_0[i][iz];
        }
        if (i != fio->n_file-1) {
            if (int(floor(fio->paras.Otime[i+1]/(2*3.141592653)))>Norb) {
                for (int iz = 0; iz != 3*Nz; iz++) {
                    overhead = 1 + iz / Nz; // leave space for entire-volume-averaged M_0
                    delimiter = (overhead-1) * (Nz+1); // for entire-volume-averaged M_0
                    fio->paras.GPME[Norb+1][iz+overhead] /= Ncount;
                    fio->paras.GPME[Norb+1][delimiter] += fio->paras.GPME[Norb+1][iz+overhead];
                    fio->paras.GPMEPar[Norb+1][iz+overhead] /= Ncount;
                    fio->paras.GPMEPar[Norb+1][delimiter] += fio->paras.GPMEPar[Norb+1][iz+overhead];
                }
                for (int j = 0; j != 3; j++) {
                    delimiter = jd1[j];
                    fio->paras.GPME[Norb+1][delimiter] /= Nz;
                    fio->paras.GPME[Norb+1][jd1[3]] += fio->paras.GPME[Norb+1][delimiter]*fio->paras.GPME[Norb+1][delimiter];
                    fio->paras.GPMEPar[Norb+1][delimiter] /= Nz;
                    fio->paras.GPMEPar[Norb+1][jd1[3]] += fio->paras.GPMEPar[Norb+1][delimiter]*fio->paras.GPMEPar[Norb+1][delimiter];
                }
                fio->paras.GPME[Norb+1][3*(Nz+1)] = sqrt(fio->paras.GPME[Norb+1][3*(Nz+1)]);
                fio->paras.GPMEPar[Norb+1][3*(Nz+1)] = sqrt(fio->paras.GPMEPar[Norb+1][3*(Nz+1)]);
                // the total momentum should = sqrt(xx+yy+zz) in the last step
                for (int iz = 0; iz != Nz; iz++) {
                    fio->paras.GPME[Norb+1][Nz1*3+1+iz] = 0;
                    fio->paras.GPMEPar[Norb+1][Nz1*3+1+iz] = 0;
                }
                for (int iz = 0; iz != Nz; iz++) {
                    fio->paras.GPME[Norb+1][Nz1*3+1+iz] = sqrt(fio->paras.GPME[Norb+1][Nz1*2+1+iz]*fio->paras.GPME[Norb+1][Nz1*2+1+iz]+fio->paras.GPME[Norb+1][Nz1*1+1+iz]*fio->paras.GPME[Norb+1][Nz1*1+1+iz]+fio->paras.GPME[Norb+1][1+iz]*fio->paras.GPME[Norb+1][1+iz]);
                    fio->paras.GPMEPar[Norb+1][Nz1*3+1+iz] = sqrt(fio->paras.GPMEPar[Norb+1][Nz1*2+1+iz]*fio->paras.GPMEPar[Norb+1][Nz1*2+1+iz]+fio->paras.GPMEPar[Norb+1][Nz1*1+1+iz]*fio->paras.GPMEPar[Norb+1][Nz1*1+1+iz]+fio->paras.GPMEPar[Norb+1][1+iz]*fio->paras.GPMEPar[Norb+1][1+iz]);
                }
                Ncount = 0;
            }
        } else {
            for (int iz = 0; iz != 3*Nz; iz++) {
                overhead = 1 + iz / Nz;
                delimiter = (overhead-1) * (Nz+1);
                fio->paras.GPME[Norb+1][iz+overhead] /= Ncount;
                fio->paras.GPME[Norb+1][delimiter] += fio->paras.GPME[Norb+1][iz+overhead];
                fio->paras.GPMEPar[Norb+1][iz+overhead] /= Ncount;
                fio->paras.GPMEPar[Norb+1][delimiter] += fio->paras.GPMEPar[Norb+1][iz+overhead];
                fio->paras.GPME[0][iz+overhead] /= fio->n_file;
                fio->paras.GPME[0][delimiter] += fio->paras.GPME[0][iz+overhead];
                fio->paras.GPMEPar[0][iz+overhead] /= fio->n_file;
                fio->paras.GPMEPar[0][delimiter] += fio->paras.GPMEPar[0][iz+overhead];
            }
            for (int j = 0; j != 3; j++) {
                delimiter = jd1[j];
                fio->paras.GPME[Norb+1][delimiter] /= Nz;
                fio->paras.GPME[Norb+1][3*(Nz+1)] += fio->paras.GPME[Norb+1][delimiter]*fio->paras.GPME[Norb+1][delimiter];
                fio->paras.GPMEPar[Norb+1][delimiter] /= Nz;
                fio->paras.GPMEPar[Norb+1][3*(Nz+1)] += fio->paras.GPMEPar[Norb+1][delimiter]*fio->paras.GPMEPar[Norb+1][delimiter];
                fio->paras.GPME[0][delimiter] /= Nz;
                fio->paras.GPME[0][3*(Nz+1)] += fio->paras.GPME[0][delimiter]*fio->paras.GPME[0][delimiter];
                fio->paras.GPMEPar[0][delimiter] /= Nz;
                fio->paras.GPMEPar[0][3*(Nz+1)] += fio->paras.GPMEPar[0][delimiter]*fio->paras.GPMEPar[0][delimiter];
            }
            fio->paras.GPME[Norb+1][3*(Nz+1)] = sqrt(fio->paras.GPME[Norb+1][3*(Nz+1)]);
            fio->paras.GPMEPar[Norb+1][3*(Nz+1)] = sqrt(fio->paras.GPMEPar[Norb+1][3*(Nz+1)]);
            fio->paras.GPME[0][3*(Nz+1)] = sqrt(fio->paras.GPME[0][3*(Nz+1)]);
            fio->paras.GPMEPar[0][3*(Nz+1)] = sqrt(fio->paras.GPMEPar[0][3*(Nz+1)]);
            // the total momentum should = sqrt(xx+yy+zz) in the last step
            for (int iz = 0; iz != Nz; iz++) {
                fio->paras.GPME[Norb+1][Nz1*3+1+iz] = 0;
                fio->paras.GPMEPar[Norb+1][Nz1*3+1+iz] = 0;
                fio->paras.GPME[0][Nz1*3+1+iz] = 0;
                fio->paras.GPMEPar[0][Nz1*3+1+iz] = 0;
            }
            for (int iz = 0; iz != Nz; iz++) {
                fio->paras.GPME[Norb+1][Nz1*3+1+iz] = sqrt(fio->paras.GPME[Norb+1][Nz1*2+1+iz]*fio->paras.GPME[Norb+1][Nz1*2+1+iz]+fio->paras.GPME[Norb+1][Nz1*1+1+iz]*fio->paras.GPME[Norb+1][Nz1*1+1+iz]+fio->paras.GPME[Norb+1][1+iz]*fio->paras.GPME[Norb+1][1+iz]);
                fio->paras.GPMEPar[Norb+1][Nz1*3+1+iz] = sqrt(fio->paras.GPMEPar[Norb+1][Nz1*2+1+iz]*fio->paras.GPMEPar[Norb+1][Nz1*2+1+iz]+fio->paras.GPMEPar[Norb+1][Nz1*1+1+iz]*fio->paras.GPMEPar[Norb+1][Nz1*1+1+iz]+fio->paras.GPMEPar[Norb+1][1+iz]*fio->paras.GPMEPar[Norb+1][1+iz]);
                fio->paras.GPME[0][Nz1*3+1+iz] = sqrt(fio->paras.GPME[0][Nz1*2+1+iz]*fio->paras.GPME[0][Nz1*2+1+iz]+fio->paras.GPME[0][Nz1*1+1+iz]*fio->paras.GPME[0][Nz1*1+1+iz]+fio->paras.GPME[0][1+iz]*fio->paras.GPME[0][1+iz]);
                fio->paras.GPMEPar[0][Nz1*3+1+iz] = sqrt(fio->paras.GPMEPar[0][Nz1*2+1+iz]*fio->paras.GPMEPar[0][Nz1*2+1+iz]+fio->paras.GPMEPar[0][Nz1*1+1+iz]*fio->paras.GPMEPar[0][Nz1*1+1+iz]+fio->paras.GPMEPar[0][1+iz]*fio->paras.GPMEPar[0][1+iz]);
            }
        }
    }
    
    Norb = 0; Ncount = 0;
    
    // !!!!!!!!!!! Second, calcualte the factor of e_0 and M_1 (we need distribute tasks to each processors)
    float temp_p_gas_i, temp_p_par_i;
    double Area =  (fio->paras.dimensions[1]*fio->paras.dimensions[0]);
    double Volume = (fio->paras.dimensions[2]*fio->paras.dimensions[1]*fio->paras.dimensions[0]);
    //cout << "Area = " << Area << ", Volume = " << Volume << endl;
    
#ifdef ENABLE_MPI
    for (int i = myMPI->loop_begin ; i <= myMPI->loop_end; i += myMPI->loop_offset) {
        cout << myMPI->prank()
#else
        for (int i = 0; i < fio->n_file; i++) {
            cout
#endif
            << "Reading for GPME " << fio->iof.data_basename+"." << setw(4) << setfill('0') << i*fio->interval+fio->start_no << endl;
            if (vf->Read_Header_Record_Pos(fio->vtk_filenames[i])) {
                cout << "Having problem reading header..." << endl;
                exit(1);
            }
            vf->Read_Data(fio->vtk_filenames[i]);
            
            // now computing M_1 and a factor of e_0
            int e0_s = 3*fourNz1;  // start index of e_0
            for (int iz = 0; iz != vf->dimensions[2]; iz++) {
                for (int iy = 0; iy != vf->dimensions[1]; iy++) {
                    for (int ix = 0; ix != vf->dimensions[0]; ix++) {
                        for (int j = 0; j != 3; j++) {
                            paras->GPME[i][e0_s+jd[j]+iz] += 1/vf->cd_scalar[0].data[iz][iy][ix]; // e_0 per dz
                            if (vf->cd_scalar[1].data[iz][iy][ix] > toe || vf->cd_scalar[1].data[iz][iy][ix] < -toe) {
                                paras->GPMEPar[i][e0_s+jd[j]+iz] += 1/vf->cd_scalar[1].data[iz][iy][ix]; // e_0 per dz
                            } // toe is to avoid inf
                            /////////////////////////////////////////////////////////////////////////////
                            temp_p_gas_i = vf->cd_vector[0].data[iz][iy][ix][j];
                            paras->GPME[i][fourNz1+jd1[j]] += (temp_p_gas_i - fio->paras.GPME[0][jd1[j]]); // <M_1(t)>_V[x,y,z,tot]
                            paras->GPME[i][fourNz1+jd1[j]+1+iz] += (temp_p_gas_i - fio->paras.GPME[0][jd1[j]+1+iz]); // <M_1(z,t)>_V[x,y,z,tot] per dz
                            temp_p_par_i = vf->cd_vector[1].data[iz][iy][ix][j];
                            paras->GPMEPar[i][fourNz1+jd1[j]] += (temp_p_par_i - fio->paras.GPMEPar[0][jd1[j]]); // <M_1(t)>_V[x,y,z,tot]
                            paras->GPMEPar[i][fourNz1+jd1[j]+1+iz] += (temp_p_par_i - fio->paras.GPMEPar[0][jd1[j]+1+iz]); // <M_1(z,t)>_V[x,y,z,tot] per dz
                        }
                    }
                }
                // <M_1(z,t)>_V[x,y,z,tot] per dz
                for (int j = 0; j != 3; j++) {
                    paras->GPME[i][e0_s+jd[j]+iz] /= Area; // e_0 per dz
                    paras->GPMEPar[i][e0_s+jd[j]+iz] /= Area; // e_0 per dz
                    /////////////////////////////////////////////////////////////////////////////
                    paras->GPME[i][fourNz1+jd1[j]+1+iz] /= Area; // <M_1(z,t)>_V[x,y,z,tot] per dz
                    paras->GPME[i][fourNz1+jd1[3]+1+iz] = paras->GPME[i][fourNz1+jd1[j]+1+iz]*paras->GPME[i][fourNz1+jd1[j]+1+iz];
                    paras->GPMEPar[i][fourNz1+jd1[j]+1+iz] /= Area;
                    paras->GPMEPar[i][fourNz1+jd1[3]+1+iz] = paras->GPMEPar[i][fourNz1+jd1[j]+1+iz]*paras->GPMEPar[i][fourNz1+jd1[j]+1+iz];
                }
                paras->GPME[i][fourNz1+jd1[3]+1+iz] = sqrt(paras->GPME[i][fourNz1+jd1[3]+1+iz]);
                paras->GPMEPar[i][fourNz1+jd1[3]+1+iz] = sqrt(paras->GPMEPar[i][fourNz1+jd1[3]+1+iz]);
            }
            
            // <M_1(t)>_V[x,y,z,tot]
            for (int j = 0; j != 3; j++) {
                paras->GPME[i][fourNz1+jd1[j]] /= Volume; // <M_1(t)>_V[x,y,z,tot]
                paras->GPME[i][fourNz1+jd1[3]] += paras->GPME[i][fourNz1+jd1[j]]*paras->GPME[i][fourNz1+jd1[j]];
                paras->GPMEPar[i][fourNz1+jd1[j]] /= Volume;
                paras->GPMEPar[i][fourNz1+jd1[3]] += paras->GPMEPar[i][fourNz1+jd1[j]]*paras->GPMEPar[i][fourNz1+jd1[j]];
            }
            paras->GPME[i][fourNz1+jd1[3]] = sqrt(paras->GPME[i][fourNz1+jd1[3]]);
            paras->GPMEPar[i][fourNz1+jd1[3]] = sqrt(paras->GPMEPar[i][fourNz1+jd1[3]]);

            // now computing M_2 and every other thing
            int e1_s = 4*fourNz1; // start index of e_1
            int e2_s = 5*fourNz1; // start index of e_2
            int M2_s = 2*fourNz1; // start index of M_2
            int m01_s = 6*fourNz1; // start index of m0*m1
            int m02_s = 7*fourNz1; // start index of m0*m2
            int m12_s = 8*fourNz1; // start index of m1*m2
            double tempM2;

            for (int iz = 0; iz != vf->dimensions[2]; iz++) {
                for (int iy = 0; iy != vf->dimensions[1]; iy++) {
                    for (int ix = 0; ix != vf->dimensions[0]; ix++) {
                        for (int j = 0; j != 3; j++) {
                            paras->GPME[i][e1_s+jd1[j]] += 0.5*paras->GPME[i][fourNz1+jd1[j]]*paras->GPME[i][fourNz1+jd1[j]]/vf->cd_scalar[0].data[iz][iy][ix]; // <e_1(t)>_V[x,y,z,tot]
                            paras->GPME[i][e1_s+jd1[j]+1+iz] += 0.5*paras->GPME[i][fourNz1+jd1[j]+1+iz]*paras->GPME[i][fourNz1+jd1[j]+1+iz]/vf->cd_scalar[0].data[iz][iy][ix]; // <e_1(z,t)>_V[x,y,z,tot]
                            if (vf->cd_scalar[1].data[iz][iy][ix] > toe || vf->cd_scalar[1].data[iz][iy][ix] < -toe) {
                                paras->GPMEPar[i][e1_s+jd1[j]] += 0.5*paras->GPMEPar[i][fourNz1+jd1[j]]*paras->GPMEPar[i][fourNz1+jd1[j]]/vf->cd_scalar[1].data[iz][iy][ix]; // <e_1(t)>_V[x,y,z,tot]
                                paras->GPMEPar[i][e1_s+jd1[j]+1+iz] += 0.5*paras->GPMEPar[i][fourNz1+jd1[j]+1+iz]*paras->GPMEPar[i][fourNz1+jd1[j]+1+iz]/vf->cd_scalar[1].data[iz][iy][ix]; // <e_1(z,t)>_V[x,y,z,tot]
                            } // toe is to avoid inf
                            /////////////////////////////////////////////////////////////////////////////
                            paras->GPME[i][m01_s+jd1[j]] += 0.5*fio->paras.GPME[0][jd1[j]]*paras->GPME[i][fourNz1+jd1[j]]/vf->cd_scalar[0].data[iz][iy][ix]; // <M_0*M_1/rho>_V(t)[x,y,z,tot]
                            paras->GPME[i][m01_s+jd1[j]+1+iz] += 0.5*fio->paras.GPME[0][jd1[j]+1+iz]*paras->GPME[i][fourNz1+jd1[j]+1+iz]/vf->cd_scalar[0].data[iz][iy][ix]; // <M_0*M_1/rho>_V(z,t)[x,y,z,tot] per dz
                            if (vf->cd_scalar[1].data[iz][iy][ix] > toe || vf->cd_scalar[1].data[iz][iy][ix] < -toe) {
                                paras->GPMEPar[i][m01_s+jd1[j]] += 0.5*fio->paras.GPMEPar[0][jd1[j]]*paras->GPMEPar[i][fourNz1+jd1[j]]/vf->cd_scalar[1].data[iz][iy][ix]; // <M_0*M_1/rho>_V(t)[x,y,z,tot]
                                paras->GPMEPar[i][m01_s+jd1[j]+1+iz] += 0.5*fio->paras.GPMEPar[0][jd1[j]+1+iz]*paras->GPMEPar[i][fourNz1+jd1[j]+1+iz]/vf->cd_scalar[1].data[iz][iy][ix]; // <M_0*M_1/rho>_V(z,t)[x,y,z,tot] per dz
                            } // toe is to avoid inf
                            /////////////////////////////////////////////////////////////////////////////
                            temp_p_gas_i = vf->cd_vector[0].data[iz][iy][ix][j];
                            tempM2 = (temp_p_gas_i - fio->paras.GPME[0][jd1[j]] - paras->GPME[i][fourNz1+jd1[j]]);
                            paras->GPME[i][M2_s+jd1[j]] += tempM2; // M_2(t)[x,y,z,tot]/V
                            paras->GPME[i][m02_s+jd1[j]] += 0.5*fio->paras.GPME[0][jd1[j]]*tempM2/vf->cd_scalar[0].data[iz][iy][ix]; // M_0*M_2/rho(t)[x,y,z,tot]/V
                            paras->GPME[i][m12_s+jd1[j]] += 0.5*paras->GPME[i][fourNz1+jd1[j]]*tempM2/vf->cd_scalar[0].data[iz][iy][ix]; // M_1*M_2/rho(t)[x,y,z,tot]/V
                            paras->GPME[i][e2_s+jd1[j]] += 0.5*tempM2*tempM2/vf->cd_scalar[0].data[iz][iy][ix]; // e_2(t)[x,y,z,tot]/V
                            temp_p_par_i = vf->cd_vector[1].data[iz][iy][ix][j];
                            tempM2 = (temp_p_par_i - fio->paras.GPMEPar[0][jd1[j]] - paras->GPMEPar[i][fourNz1+jd1[j]]);
                            paras->GPMEPar[i][M2_s+jd1[j]] += tempM2; // M_2(t)[x,y,z,tot]/V
                            if (vf->cd_scalar[1].data[iz][iy][ix] > toe || vf->cd_scalar[1].data[iz][iy][ix] < -toe) {
                                paras->GPMEPar[i][m02_s+jd1[j]] += 0.5*fio->paras.GPMEPar[0][jd1[j]]*tempM2/vf->cd_scalar[1].data[iz][iy][ix]; // M_0*M_2/rho(t)[x,y,z,tot]/V
                                paras->GPMEPar[i][m12_s+jd1[j]] += 0.5*paras->GPMEPar[i][fourNz1+jd1[j]]*tempM2/vf->cd_scalar[1].data[iz][iy][ix]; // M_1*M_2/rho(t)[x,y,z,tot]/V
                                paras->GPMEPar[i][e2_s+jd1[j]] += 0.5*tempM2*tempM2/vf->cd_scalar[1].data[iz][iy][ix]; // e_2(t)[x,y,z,tot]/V
                            } // toe is to avoid inf
                            /////////////////////////////////////////////////////////////////////////////
                            tempM2 =  (temp_p_gas_i - fio->paras.GPME[0][jd1[j]+1+iz] - paras->GPME[i][fourNz1+jd1[j]+1+iz]);
                            paras->GPME[i][M2_s+jd1[j]+1+iz] += tempM2; // M_2(z,t)[x,y,z,tot]/A
                            paras->GPME[i][e2_s+jd1[j]+1+iz] += 0.5*tempM2*tempM2/vf->cd_scalar[0].data[iz][iy][ix]; // e_2(z,t)[x,y,z,tot]/A,
                            paras->GPME[i][m02_s+jd1[j]+1+iz] += 0.5*fio->paras.GPME[0][jd[j]+1+iz]*tempM2/vf->cd_scalar[0].data[iz][iy][ix]; // M_0*M_2/rho(z,t)[x,y,z,tot]/A
                            paras->GPME[i][m12_s+jd1[j]+1+iz] += 0.5*paras->GPME[i][fourNz1+jd1[1]+1+iz]*tempM2/vf->cd_scalar[0].data[iz][iy][ix]; // M_1*M_2/rho(z,t)[x,y,z,tot]/A
                            tempM2 =  (temp_p_par_i - fio->paras.GPMEPar[0][jd1[j]+1+iz] - paras->GPMEPar[i][fourNz1+jd1[j]+1+iz]);
                            paras->GPMEPar[i][M2_s+jd1[j]+1+iz] += tempM2; // M_2(z,t)[x,y,z,tot]/A
                            if (vf->cd_scalar[1].data[iz][iy][ix] > toe || vf->cd_scalar[1].data[iz][iy][ix] < -toe) {
                                paras->GPMEPar[i][e2_s+jd1[j]+1+iz] += 0.5*tempM2*tempM2/vf->cd_scalar[1].data[iz][iy][ix]; // e_2(z,t)[x,y,z,tot]/A,
                                paras->GPMEPar[i][m02_s+jd1[j]+1+iz] += 0.5*fio->paras.GPMEPar[0][jd[j]+1+iz]*tempM2/vf->cd_scalar[1].data[iz][iy][ix]; // M_0*M_2/rho(z,t)[x,y,z,tot]/A
                                paras->GPMEPar[i][m12_s+jd1[j]+1+iz] += 0.5*paras->GPMEPar[i][fourNz1+jd1[1]+1+iz]*tempM2/vf->cd_scalar[1].data[iz][iy][ix]; // M_1*M_2/rho(z,t)[x,y,z,tot]/A
                            } // toe is to avoid inf
                        }
                    }
                }
                
                for (int j = 0; j != 3; j++) {
                    paras->GPME[i][e1_s+jd1[j]+1+iz] /= Area;
                    paras->GPME[i][e1_s+jd1[3]+1+iz] += paras->GPME[i][e1_s+jd1[j]+1+iz];
                    paras->GPME[i][M2_s+jd1[j]+1+iz] /= Area;
                    paras->GPME[i][M2_s+jd1[3]+1+iz] += paras->GPME[i][M2_s+jd1[j]+1+iz]*paras->GPME[i][M2_s+jd1[j]+1+iz];
                    paras->GPME[i][e2_s+jd1[j]+1+iz] /= Area;
                    paras->GPME[i][e2_s+jd1[3]+1+iz] += paras->GPME[i][e2_s+jd1[j]+1+iz];
                    paras->GPME[i][m01_s+jd1[j]+1+iz] /= Area;
                    paras->GPME[i][m01_s+jd1[3]+1+iz] += paras->GPME[i][m01_s+jd1[j]+1+iz]*paras->GPME[i][m01_s+jd1[j]+1+iz];
                    paras->GPME[i][m02_s+jd1[j]+1+iz] /= Area;
                    paras->GPME[i][m02_s+jd1[3]+1+iz] += paras->GPME[i][m02_s+jd1[j]+1+iz]*paras->GPME[i][m02_s+jd1[j]+1+iz];
                    paras->GPME[i][m12_s+jd1[j]+1+iz] /= Area;
                    paras->GPME[i][m12_s+jd1[3]+1+iz] += paras->GPME[i][m12_s+jd1[j]+1+iz]*paras->GPME[i][m12_s+jd1[j]+1+iz];
                    paras->GPMEPar[i][e1_s+jd1[j]+1+iz] /= Area;
                    paras->GPMEPar[i][e1_s+jd1[3]+1+iz] += paras->GPMEPar[i][e1_s+jd1[j]+1+iz];
                    paras->GPMEPar[i][M2_s+jd1[j]+1+iz] /= Area;
                    paras->GPMEPar[i][M2_s+jd1[3]+1+iz] += paras->GPMEPar[i][M2_s+jd1[j]+1+iz]*paras->GPMEPar[i][M2_s+jd1[j]+1+iz];
                    paras->GPMEPar[i][e2_s+jd1[j]+1+iz] /= Area;
                    paras->GPMEPar[i][e2_s+jd1[3]+1+iz] += paras->GPMEPar[i][e2_s+jd1[j]+1+iz];
                    paras->GPMEPar[i][m01_s+jd1[j]+1+iz] /= Area;
                    paras->GPMEPar[i][m01_s+jd1[3]+1+iz] += paras->GPMEPar[i][m01_s+jd1[j]+1+iz]*paras->GPMEPar[i][m01_s+jd1[j]+1+iz];
                    paras->GPMEPar[i][m02_s+jd1[j]+1+iz] /= Area;
                    paras->GPMEPar[i][m02_s+jd1[3]+1+iz] += paras->GPMEPar[i][m02_s+jd1[j]+1+iz]*paras->GPMEPar[i][m02_s+jd1[j]+1+iz];
                    paras->GPMEPar[i][m12_s+jd1[j]+1+iz] /= Area;
                    paras->GPMEPar[i][m12_s+jd1[3]+1+iz] += paras->GPMEPar[i][m12_s+jd1[j]+1+iz]*paras->GPMEPar[i][m12_s+jd1[j]+1+iz];
                }
                paras->GPME[i][M2_s+jd1[3]+1+iz] = sqrt(paras->GPME[i][M2_s+jd1[3]+1+iz]);
                paras->GPME[i][m01_s+jd1[3]+1+iz] = sqrt(paras->GPME[i][m01_s+jd1[3]+1+iz]);
                paras->GPME[i][m02_s+jd1[3]+1+iz] = sqrt(paras->GPME[i][m02_s+jd1[3]+1+iz]);
                paras->GPME[i][m12_s+jd1[3]+1+iz] = sqrt(paras->GPME[i][m12_s+jd1[3]+1+iz]);
                paras->GPMEPar[i][M2_s+jd1[3]+1+iz] = sqrt(paras->GPMEPar[i][M2_s+jd1[3]+1+iz]);
                paras->GPMEPar[i][m01_s+jd1[3]+1+iz] = sqrt(paras->GPMEPar[i][m01_s+jd1[3]+1+iz]);
                paras->GPMEPar[i][m02_s+jd1[3]+1+iz] = sqrt(paras->GPMEPar[i][m02_s+jd1[3]+1+iz]);
                paras->GPMEPar[i][m12_s+jd1[3]+1+iz] = sqrt(paras->GPMEPar[i][m12_s+jd1[3]+1+iz]);
            }
            
            for (int j = 0; j != 3; j++) {
                paras->GPME[i][M2_s+jd1[j]] /= Volume;
                paras->GPME[i][M2_s+jd1[3]] += paras->GPME[i][M2_s+jd1[j]]*paras->GPME[i][M2_s+jd1[j]];
                paras->GPME[i][e1_s+jd1[j]] /= Volume;
                paras->GPME[i][e1_s+jd1[3]] += paras->GPME[i][e1_s+jd1[j]];
                paras->GPME[i][e2_s+jd1[j]] /= Volume;
                paras->GPME[i][e2_s+jd1[3]] += paras->GPME[i][e2_s+jd1[j]];
                paras->GPME[i][m01_s+jd1[j]] /= Volume;
                paras->GPME[i][m01_s+jd1[3]] += paras->GPME[i][m01_s+jd1[j]]*paras->GPME[i][m01_s+jd1[j]];
                paras->GPME[i][m02_s+jd1[j]] /= Volume;
                paras->GPME[i][m02_s+jd1[3]] += paras->GPME[i][m02_s+jd1[j]]*paras->GPME[i][m02_s+jd1[j]];
                paras->GPME[i][m12_s+jd1[j]] /= Volume;
                paras->GPME[i][m12_s+jd1[3]] += paras->GPME[i][m12_s+jd1[j]]*paras->GPME[i][m12_s+jd1[j]];
                paras->GPMEPar[i][M2_s+jd1[j]] /= Volume;
                paras->GPMEPar[i][M2_s+jd1[3]] += paras->GPMEPar[i][M2_s+jd1[j]]*paras->GPMEPar[i][M2_s+jd1[j]];
                paras->GPMEPar[i][e1_s+jd1[j]] /= Volume;
                paras->GPMEPar[i][e1_s+jd1[3]] += paras->GPMEPar[i][e1_s+jd1[j]];
                paras->GPMEPar[i][e2_s+jd1[j]] /= Volume;
                paras->GPMEPar[i][e2_s+jd1[3]] += paras->GPMEPar[i][e2_s+jd1[j]];
                paras->GPMEPar[i][m01_s+jd1[j]] /= Volume;
                paras->GPMEPar[i][m01_s+jd1[3]] += paras->GPMEPar[i][m01_s+jd1[j]]*paras->GPMEPar[i][m01_s+jd1[j]];
                paras->GPMEPar[i][m02_s+jd1[j]] /= Volume;
                paras->GPMEPar[i][m02_s+jd1[3]] += paras->GPMEPar[i][m02_s+jd1[j]]*paras->GPMEPar[i][m02_s+jd1[j]];
                paras->GPMEPar[i][m12_s+jd1[j]] /= Volume;
                paras->GPMEPar[i][m12_s+jd1[3]] += paras->GPMEPar[i][m12_s+jd1[j]]*paras->GPMEPar[i][m12_s+jd1[j]];
            }
            paras->GPME[i][M2_s+jd1[3]] = sqrt(paras->GPME[i][M2_s+jd1[3]]);
            paras->GPME[i][m01_s+jd1[3]] = sqrt(paras->GPME[i][m01_s+jd1[3]]);
            paras->GPME[i][m02_s+jd1[3]] = sqrt(paras->GPME[i][m02_s+jd1[3]]);
            paras->GPME[i][m12_s+jd1[3]] = sqrt(paras->GPME[i][m12_s+jd1[3]]);
            paras->GPMEPar[i][M2_s+jd1[3]] = sqrt(paras->GPMEPar[i][M2_s+jd1[3]]);
            paras->GPMEPar[i][m01_s+jd1[3]] = sqrt(paras->GPMEPar[i][m01_s+jd1[3]]);
            paras->GPMEPar[i][m02_s+jd1[3]] = sqrt(paras->GPMEPar[i][m02_s+jd1[3]]);
            paras->GPMEPar[i][m12_s+jd1[3]] = sqrt(paras->GPMEPar[i][m12_s+jd1[3]]);
            
#ifdef ENABLE_MPI
        }
#else
    }
#endif
    
    // !!!!!!!!!!! Third, if MPI collect data, then calculate e_0
#ifdef ENABLE_MPI
    for (int i = 0; i != fio->n_file; i++) {
        MPI::COMM_WORLD.Allreduce(&(myMPI->paras.GPME[i][fourNz1]), &(fio->paras.GPME[i][fourNz1]), fourNz1*8, MPI::DOUBLE, MPI::SUM);
        MPI::COMM_WORLD.Allreduce(&(myMPI->paras.GPMEPar[i][fourNz1]), &(fio->paras.GPMEPar[i][fourNz1]), fourNz1*8, MPI::DOUBLE, MPI::SUM);
    }
#endif
    int e0_s = 3*fourNz1;  // start index of e_0
    double **e_0 = new double*[fio->n_file];
    double **eP_0 = new double*[fio->n_file];
    for (int i = 0; i != fio->n_file; i++) {
        e_0[i] = new double[fourNz];
        eP_0[i] = new double[fourNz];
        for (int j = 0; j != fourNz; j++) {
            e_0[i][j] = fio->paras.GPME[i][e0_s+j];
            eP_0[i][j] = fio->paras.GPMEPar[i][e0_s+j];
            fio->paras.GPME[i][e0_s+j] = 0;
            fio->paras.GPMEPar[i][e0_s+j] = 0;
#ifdef ENABLE_MPI
            // we need read again and addreduce again
            myMPI->paras.GPME[i][e0_s+j] = 0;
            myMPI->paras.GPMEPar[i][e0_s+j] = 0;
#endif
        }
    }

    for (int i = 0; i != fio->n_file; i++) {
        Norb = int(floor(fio->paras.Otime[i]/(2*3.141592653)));
        Ncount++;
        for (int iz = 0; iz != 3*Nz; iz++) {
            overhead = 1 + iz / Nz; // leave space for entire-volume-averaged e_0
            fio->paras.GPME[0][e0_s+iz+overhead] += 0.5*fio->paras.GPME[0][iz+overhead]*fio->paras.GPME[0][iz+overhead]*e_0[i][iz];
            fio->paras.GPME[Norb+1][e0_s+iz+overhead] += 0.5*fio->paras.GPME[Norb+1][iz+overhead]*fio->paras.GPME[Norb+1][iz+overhead]*e_0[i][iz];
            fio->paras.GPMEPar[0][e0_s+iz+overhead] += 0.5*fio->paras.GPMEPar[0][iz+overhead]*fio->paras.GPMEPar[0][iz+overhead]*eP_0[i][iz];
            fio->paras.GPMEPar[Norb+1][e0_s+iz+overhead] += 0.5*fio->paras.GPMEPar[Norb+1][iz+overhead]*fio->paras.GPMEPar[Norb+1][iz+overhead]*eP_0[i][iz];
        }
        if (i != fio->n_file-1) {
            if (int(floor(fio->paras.Otime[i+1]/(2*3.141592653)))>Norb) {
                for (int iz = 0; iz != 3*Nz; iz++) {
                    overhead = 1 + iz / Nz; // leave space for entire-volume-averaged e_0
                    delimiter = (overhead-1) * (Nz+1); // for entire-volume-averaged e_0
                    fio->paras.GPME[Norb+1][e0_s+iz+overhead] /= Ncount;
                    fio->paras.GPME[Norb+1][e0_s+delimiter] += fio->paras.GPME[Norb+1][e0_s+iz+overhead];
                    fio->paras.GPMEPar[Norb+1][e0_s+iz+overhead] /= Ncount;
                    fio->paras.GPMEPar[Norb+1][e0_s+delimiter] += fio->paras.GPMEPar[Norb+1][e0_s+iz+overhead];
                }
                for (int j = 0; j != 3; j++) {
                    delimiter = jd1[j];
                    fio->paras.GPME[Norb+1][e0_s+delimiter] /= Nz;
                    fio->paras.GPME[Norb+1][e0_s+jd1[3]] += fio->paras.GPME[Norb+1][e0_s+delimiter];
                    fio->paras.GPMEPar[Norb+1][e0_s+delimiter] /= Nz;
                    fio->paras.GPMEPar[Norb+1][e0_s+jd1[3]] += fio->paras.GPMEPar[Norb+1][e0_s+delimiter];
                }
                
                // the total energy = x+y+z in the last step
                for (int iz = 0; iz != Nz; iz++) {
                    fio->paras.GPME[Norb+1][e0_s+jd1[3]+1+iz] = 0;
                    fio->paras.GPMEPar[Norb+1][e0_s+jd1[3]+1+iz] = 0;
                }
                for (int iz = 0; iz != Nz; iz++) {
                    fio->paras.GPME[Norb+1][e0_s+jd1[3]+1+iz] = fio->paras.GPME[Norb+1][e0_s+jd1[2]+1+iz]+fio->paras.GPME[Norb+1][e0_s+jd1[1]+1+iz]+fio->paras.GPME[Norb+1][e0_s+1+iz];
                    fio->paras.GPMEPar[Norb+1][e0_s+jd1[3]+1+iz] = fio->paras.GPMEPar[Norb+1][e0_s+jd1[2]+1+iz]+fio->paras.GPMEPar[Norb+1][e0_s+jd1[1]+1+iz]+fio->paras.GPMEPar[Norb+1][e0_s+1+iz];
                }
                Ncount = 0;
            }
        } else {
            for (int iz = 0; iz != 3*Nz; iz++) {
                overhead = 1 + iz / Nz;
                delimiter = (overhead-1) * (Nz+1);
                fio->paras.GPME[Norb+1][e0_s+iz+overhead] /= Ncount;
                fio->paras.GPME[Norb+1][e0_s+delimiter] += fio->paras.GPME[Norb+1][e0_s+iz+overhead];
                fio->paras.GPMEPar[Norb+1][e0_s+iz+overhead] /= Ncount;
                fio->paras.GPMEPar[Norb+1][e0_s+delimiter] += fio->paras.GPMEPar[Norb+1][e0_s+iz+overhead];
                fio->paras.GPME[0][e0_s+iz+overhead] /= fio->n_file;
                fio->paras.GPME[0][e0_s+delimiter] += fio->paras.GPME[0][e0_s+iz+overhead];
                fio->paras.GPMEPar[0][e0_s+iz+overhead] /= fio->n_file;
                fio->paras.GPMEPar[0][e0_s+delimiter] += fio->paras.GPMEPar[0][e0_s+iz+overhead];
            }
            for (int j = 0; j != 3; j++) {
                delimiter = jd1[j];
                fio->paras.GPME[Norb+1][e0_s+delimiter] /= Nz;
                fio->paras.GPME[Norb+1][e0_s+jd1[3]] += fio->paras.GPME[Norb+1][e0_s+delimiter];
                fio->paras.GPMEPar[Norb+1][e0_s+delimiter] /= Nz;
                fio->paras.GPMEPar[Norb+1][e0_s+jd1[3]] += fio->paras.GPMEPar[Norb+1][e0_s+delimiter];
                fio->paras.GPME[0][e0_s+delimiter] /= Nz;
                fio->paras.GPME[0][e0_s+jd1[3]] += fio->paras.GPME[0][e0_s+delimiter];
                fio->paras.GPMEPar[0][e0_s+delimiter] /= Nz;
                fio->paras.GPMEPar[0][e0_s+jd1[3]] += fio->paras.GPMEPar[0][e0_s+delimiter];
            }
            
            // the total energy = x+y+z in the last step
            for (int iz = 0; iz != Nz; iz++) {
                fio->paras.GPME[Norb+1][e0_s+Nz1*3+1+iz] = 0;
                fio->paras.GPMEPar[Norb+1][e0_s+Nz1*3+1+iz] = 0;
                fio->paras.GPME[0][e0_s+Nz1*3+1+iz] = 0;
                fio->paras.GPMEPar[0][e0_s+Nz1*3+1+iz] = 0;
            }
            for (int iz = 0; iz != Nz; iz++) {
                fio->paras.GPME[Norb+1][e0_s+Nz1*3+1+iz] = fio->paras.GPME[Norb+1][e0_s+Nz1*2+1+iz]+fio->paras.GPME[Norb+1][e0_s+Nz1*1+1+iz]+fio->paras.GPME[Norb+1][e0_s+1+iz];
                fio->paras.GPMEPar[Norb+1][e0_s+Nz1*3+1+iz] = fio->paras.GPMEPar[Norb+1][e0_s+Nz1*2+1+iz]+fio->paras.GPMEPar[Norb+1][e0_s+Nz1*1+1+iz]+fio->paras.GPMEPar[Norb+1][e0_s+1+iz];
                fio->paras.GPME[0][e0_s+Nz1*3+1+iz] = fio->paras.GPME[0][e0_s+Nz1*2+1+iz]+fio->paras.GPME[0][e0_s+Nz1*1+1+iz]+fio->paras.GPME[0][e0_s+1+iz];
                fio->paras.GPMEPar[0][e0_s+Nz1*3+1+iz] = fio->paras.GPMEPar[0][e0_s+Nz1*2+1+iz]+fio->paras.GPMEPar[0][e0_s+Nz1*1+1+iz]+fio->paras.GPMEPar[0][e0_s+1+iz];
            }
        }

    }

    /************************** GPME **************************/
    return 0;
}

int wangbadan(int exitcode)
{
    for (int i = 0; i != fio->n_file; i++) {
        for (int j = 0; j != 4*9*(fio->paras.dimensions[2]+1); j++) {
            if (std::isnan(fio->paras.GPME[i][j])) {
                cout << "Found nan GPME[" << i << "][" << j << "]\n";
                exit(exitcode);
            }
            if (std::isinf(fio->paras.GPME[i][j])) {
                cout << "Found inf GPME[" << i << "][" << j << "]\n";
                exit(exitcode);
            }
            if (std::isnan(fio->paras.GPMEPar[i][j])) {
                cout << "Found nan GPMEPar[" << i << "][" << j << "]\n";
                exit(exitcode);
            }
            if (std::isinf(fio->paras.GPMEPar[i][j])) {
                cout << "Found inf GPMEPar[" << i << "][" << j << "]\n";
                exit(exitcode);
            }
        }
    }
    return 0;
}
