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
    
    // we need local variables to record data and then MPI_Allreduce
    double *orbit_time = new double[fio->n_file];
    double *max_rho_par = new double[fio->n_file];
    double *Hp = new double[fio->n_file];
    long *n_par = new long[fio->n_file];
    long **cpuid_dist = new long*[fio->n_file];
    for (int i = 0; i != fio->n_file; i++) {
        // initialize if you don't assign all of them values but use them for calculation
        orbit_time[i] = 0;
        max_rho_par[i] = 0;
        Hp[i] = 0;
        n_par[i] = 0;
        cpuid_dist[i] = new long[fio->n_cpu];
        for (int j = 0; j != fio->n_cpu; j++) {
            cpuid_dist[i][j] = 0;
        }
    }
    myMPI->Barrier();
    // for debug
    cout << "Processor " << myMPI->myrank << ": " << myMPI->loop_begin << " " << myMPI->loop_end << " " << myMPI->loop_offset << endl;

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
#ifndef ENABLE_MPI
        //vf->Print_File_Info();
#endif
        if (fio->RhoParMax_flag) {
            vf->Read_Data(fio->vtk_filenames[i]);
            vf->Calculate_Mass_Find_Max();
        }
        if (fio->ParNum_flag || fio->HeiPar_flag || fio->CpuID_flag) {
            // lis part
            pl->ReadLis(fio->lis_filenames[i]);
        }
                
        // recording data
#ifdef ENABLE_MPI
        orbit_time[i] = vf->time;
        if (fio->ParNum_flag) {
            n_par[i] = pl->n;
        }
        if (fio->RhoParMax_flag) {
            max_rho_par[i] = vf->max_rho_par;
        }
        if (fio->HeiPar_flag) {
            Hp[i] = pl->ScaleHeight();
        }
        if (fio->CpuID_flag) {
            pl->CpuID();
            delete [] cpuid_dist[i];
            cpuid_dist[i] = pl->CpuID_dist;
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
            fio->Hp[i] = pl->ScaleHeight();
        }
        if (fio->CpuID_flag) {
            pl->CpuID();
            fio->CpuID_dist[i] = pl->CpuID_dist;
        }

#endif
#ifndef RESIZE_LIST
        pl->InitializeList();
#endif
    }
    
#ifdef ENABLE_MPI
    myMPI->Barrier();
    MPI::COMM_WORLD.Allreduce(orbit_time, fio->orbit_time, fio->n_file, MPI::DOUBLE, MPI::SUM);
    if (fio->ParNum_flag) {
        MPI::COMM_WORLD.Allreduce(n_par, fio->n_par, fio->n_file, MPI::LONG, MPI::SUM);
    }
    if (fio->RhoParMax_flag) {
        MPI::COMM_WORLD.Allreduce(max_rho_par, fio->max_rho_par, fio->n_file, MPI::DOUBLE, MPI::SUM);
    }
    if (fio->HeiPar_flag) {
        MPI::COMM_WORLD.Allreduce(Hp, fio->Hp, fio->n_file, MPI::DOUBLE, MPI::SUM);
    }
    if (fio->CpuID_flag) {
        for (int i = 0; i != fio->n_file; i++) {
            fio->CpuID_dist[i] = new long[fio->n_cpu];
            for (int j = 0; j != fio->n_cpu; j++) {
                fio->CpuID_dist[i][j] = 0;
            }
            MPI::COMM_WORLD.Allreduce(cpuid_dist[i], fio->CpuID_dist[i], fio->n_cpu, MPI::LONG, MPI::SUM);
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
    
    delete [] orbit_time;
    delete [] max_rho_par;
    delete [] Hp;
    delete [] cpuid_dist;
    delete myMPI;
    
#endif
    return 0;
}
