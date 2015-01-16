//
//  main.cpp
//  readvtklis
//
//  Created by Rixin Li on 1/14/15.
//  Copyright (c) 2015 Rixin Li. All rights reserved.
//

#include <iostream>
#include "fop.h"
#include "readlis.h"
#include "readvtk.h"

using namespace std;

int main(int argc, const char * argv[]) {
    
    FileIO *fio = new FileIO(argc, argv);
    fio->Generate_Filename();
    ParticleList *pl = new ParticleList;

    for (int i = fio->start_no; i <= fio->end_no; i++) {
        pl->ReadLis(fio->lis_filenames[i]);
        cout << pl->ScaleHeight() << endl;
        pl->InitializeList();
    }
    
    
    delete fio;
    delete pl;

    return 0;
}
