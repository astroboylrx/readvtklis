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
        
        pl->ReadLis(fio->lis_filenames[i]);
        cout << pl->ScaleHeight() << endl;
#ifndef RESIZE_LIST
        pl->InitializeList();
#endif
        if (vf->Read_Header_Record_Pos(fio->vtk_filenames[i])) {
            cout << "Having problem reading header..." << endl;
            exit(1);
        }
        vf->cd_scaler[0].Read_Scaler_Data(fio->vtk_filenames[i]);
        vf->cd_scaler[1].Read_Scaler_Data(fio->vtk_filenames[i]);
        
        float temp = 0;
        for (int i = 0; i != vf->dimensions[2]; i++) {
            for (int j = 0; j != vf->dimensions[1]; j++) {
                for (int k = 0; k != vf->dimensions[0]; k++) {
                    if (vf->cd_scaler[1].data[i][j][k] > temp) {
                        temp = vf->cd_scaler[1].data[i][j][k];
                    }
                }
            }
        }
        cout << "max of p density: " << temp << endl;
        
    }
    
    delete fio;
    delete pl;
    delete vf;
    
    end_t = clock();
    elapsed_secs = double(end_t - begin_t) / CLOCKS_PER_SEC;
    cout << "Elapsed time (secs) is " << elapsed_secs << endl;
    
    return 0;
}
