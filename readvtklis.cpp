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
        cout << "file number: " << fio->start_no+i << endl;
        pl->ReadLis(fio->lis_filenames[i]);
        cout << pl->ScaleHeight() << endl;
#ifndef RESIZE_LIST
        pl->InitializeList();
#endif
        if (vf->Read_Header_Record_Pos(fio->vtk_filenames[i])) {
            cout << "Having problem reading header..." << endl;
            exit(1);
        }
        vf->Read_Data(fio->vtk_filenames[i]);
    }
    fio->print_stars("Finishing Program...");
    delete fio;
    delete pl;
    delete vf;
    
    end_t = clock();
    elapsed_secs = double(end_t - begin_t) / CLOCKS_PER_SEC;
    cout << "Elapsed time (secs) is " << elapsed_secs << endl;
    
    return 0;
}
