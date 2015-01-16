//
//  main.cpp
//  readvtklis
//
//  Created by Rixin Li on 1/14/15.
//  Copyright (c) 2015 Rixin Li. All rights reserved.
//

#include <iostream>
#include "readlis.h"
#include "readvtk.h"

using namespace std;

int main(int argc, const char * argv[]) {
    
    if (argc < 3) {
        printf("USAGE: %s -i <data_path> -b <data_basename> -s <post-name> -f <# (range(f1:f2)) -o <output_path_name>\n", argv[0]);
        exit(1);
    } else {
        switch (argc)
        {
            case 3:
                if (strcmp(argv[1], "-c") == 0 || strcmp(argv[1], "c") == 0) {
                    fop->conti = 1;
                    strcpy(fop->outputpath, argv[2]);
                    break;
                }
                strcpy(fop->outputpath, argv[1]);
                strcpy(fop->cont, argv[2]);
                if (strcmp(fop->cont, "-c") == 0 || strcmp(fop->cont, "c") == 0) {
                    fop->conti = 1;
                } else {
                    fop->conti = 0;
                }
                //printf("fop->cont = %s, fop->conti = %d, argv[2] = %s, argv[1] = %s\n", fop->cont, fop->conti, argv[2], argv[1]);
                break;
                /********** waiting for implementation **********/
                
            default:
                printf("Only one arguments: USAGE: %s, <data_output_path> [-c]\n", argv[0]);
                exit(1);
        }
    }
    
    ParticleList *pl;
    pl = (ParticleList *)malloc(sizeof(ParticleList));

    
    
    
    cout << "Hello, World!\n";
    return 0;
}
