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

int main(int argc, const char * argv[]) {
    
    clock_t begin_t, end_t;
    double elapsed_secs, mpi_begin_t, mpi_end_t;
    begin_t = clock();
#ifdef ENABLE_MPI
    MPI_info *myMPI = new MPI_info;
    myMPI->Initialize(argc, argv);
    mpi_begin_t = MPI::Wtime();
#endif
    FileIO *fio = new FileIO(argc, argv);
    fio->Generate_Filename();
    ParticleList *pl = new ParticleList;
    VtkFile *vf = new VtkFile;
#ifdef ENABLE_MPI
    if (myMPI->myrank == myMPI->master) {
#endif
        fio->Check_Path_Filename();
        fio->print_stars("Master: Begin to Process Data");
#ifdef ENABLE_MPI
    }
    myMPI->Determine_Loop(fio->n_file);
    double *orbit_time = new double[fio->n_file];
    double *max_rho_par = new double[fio->n_file];
    double *Hp = new double[fio->n_file];
    myMPI->Barrier();
    // for debug
    //cout << "Processor " << myMPI->myrank << ": " << myMPI->loop_begin << " " << myMPI->loop_end << " " << myMPI->loop_offset << endl;

    for (int i = myMPI->loop_begin ; i <= myMPI->loop_end; i += myMPI->loop_offset) {
#else
    for (int i = 0; i < fio->n_file; i++) {
#endif
        // vtk part
        cout
#ifdef ENABLE_MPI
        << "Processor " << myMPI->myrank << ": "
#endif
        << "Reading " << fio->data_basename+"." << setw(4) << setfill('0') << i+fio->start_no << endl;
        
        if (vf->Read_Header_Record_Pos(fio->vtk_filenames[i])) {
            cout << "Having problem reading header..." << endl;
            exit(1);
        }
#ifndef ENABLE_MPI
        //vf->Print_File_Info();
#endif
        vf->Read_Data(fio->vtk_filenames[i]);
        vf->Calculate_Mass_Find_Max();

        // lis part
        pl->ReadLis(fio->lis_filenames[i]);
        
        // recording data
#ifdef ENABLE_MPI
        orbit_time[i] = vf->time;
        max_rho_par[i] = vf->max_rho_par*fio->mratio*vf->m_gas/vf->m_par;
        Hp[i] = pl->ScaleHeight();
#else
        fio->orbit_time[i] = vf->time;
        // rescale the total mass of particles to 0.02 gas mass, so this density is relative to gas density 1
        fio->max_rho_par[i] = vf->max_rho_par*fio->mratio*vf->m_gas/vf->m_par;
        fio->Hp[i] = pl->ScaleHeight();
#endif
#ifndef RESIZE_LIST
        pl->InitializeList();
#endif
    }
    
#ifdef ENABLE_MPI
    myMPI->Barrier();
    MPI::COMM_WORLD.Allreduce(orbit_time, fio->orbit_time, fio->n_file, MPI::DOUBLE, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(max_rho_par, fio->max_rho_par, fio->n_file, MPI::DOUBLE, MPI::SUM);
    MPI::COMM_WORLD.Allreduce(Hp, fio->Hp, fio->n_file, MPI::DOUBLE, MPI::SUM);
    cout << "Processor " << myMPI->myrank << ": I'm done." << endl;
    myMPI->Barrier();
    if (myMPI->myrank == myMPI->master) {
#endif
        fio->output_data();
        /*
        for (int i = 0; i != fio->n_file; i++) {
            cout << "time = " << fio->orbit_time[i] << "; max_rho_par = " << fio->max_rho_par[i] << "; Hp = " << fio->Hp[i] << endl;
        }
         */
        fio->print_stars("Master: Finishing Program");
#ifdef ENABLE_MPI
    }
#endif
    
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
    delete myMPI;
#endif
    return 0;
}
