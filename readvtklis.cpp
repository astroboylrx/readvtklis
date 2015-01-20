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
    double elapsed_secs;
    begin_t = clock();
    
    FileIO *fio = new FileIO(argc, argv);
    fio->Generate_Filename();
    ParticleList *pl = new ParticleList;
    VtkFile *vf = new VtkFile;

    fio->print_stars("Begin to Process Data");
    for (int i = 0 ; i <= (fio->end_no - fio->start_no); i++) {
        // vtk part
        cout << "\nReading " << fio->vtk_filenames[i].substr(fio->vtk_filenames[i].find_last_of("/\\")+1) << endl;
        
        if (vf->Read_Header_Record_Pos(fio->vtk_filenames[i])) {
            cout << "Having problem reading header..." << endl;
            exit(1);
        }
        vf->Print_File_Info();
        vf->Construct_Coor();
        vf->Read_Data(fio->vtk_filenames[i]);
        vf->Calculate_Mass_Find_Max();
        
        cout << "\nReading " << fio->lis_filenames[i].substr(fio->lis_filenames[i].find_last_of("/\\")+1) << endl;
        // lis part
        pl->ReadLis(fio->lis_filenames[i]);
        
        // recording data
        fio->orbit_time[i] = vf->time;
        // rescale the total mass of particles to 0.02 gas mass, so this density is relative to gas density 1
        fio->max_rho_par[i] = vf->max_rho_par*fio->mratio*vf->m_gas/vf->m_par;
        fio->Hp[i] = pl->ScaleHeight();
#ifndef RESIZE_LIST
        pl->InitializeList();
#endif
    }
    fio->output_data();
    
    fio->print_stars("Finishing Program...");
    delete fio;
    delete pl;
    delete vf;
    
    end_t = clock();
    elapsed_secs = double(end_t - begin_t) / CLOCKS_PER_SEC;
    cout << "Elapsed time (secs) is " << elapsed_secs << endl;
    
    return 0;
}
