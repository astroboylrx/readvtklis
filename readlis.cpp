//
//  readlis.cpp
//  readvtklis
//
//  Created by Rixin Li on 1/14/15.
//  Copyright (c) 2015 Rixin Li. All rights reserved.
//

#include "readlis.h"

/********** Read particle list from file **********/
int ParticleList::ReadLis(string filename)
{
    ifstream file (filename.c_str(), ios::binary);
    if (file.is_open()) {
        file.read((char *)(coorlim), 12*sizeof(float));
        // you can read an array in a bundle, like above
        // The values in coorlim are the coordinate limits of this grid and domain:
        // x1l, x1u, x2l, x2u, x3l, x3u, x1dl, x1du, x2dl, x2du, x3dl, x3du
        // here l means lower limit, u means upper limit, d means domain
        
        file.read((char *)(&ntype), sizeof(float));
        typeinfo = (float *)malloc(ntype*sizeof(float));
        for (int i = 0; i < ntype; i++) {
            file.read((char *)(&typeinfo[i]), sizeof(float));
        }
        file.read((char *)(&time), sizeof(float));
        file.read((char *)(&dt), sizeof(float));
        file.read((char *)(&n), sizeof(long));
        List = (Particle *)malloc(n*sizeof(Particle));
        for (long i = 0; i < n; i++) {
            for (int j = 0; j < 3; j++) {
                file.read((char *)(&List[i].x[j]), sizeof(float));
            }
            for (int j = 0; j < 3; j++) {
                file.read((char *)(&List[i].v[j]), sizeof(float));
            }
            file.read((char *)(&List[i].rad), sizeof(float));
            file.read((char *)(&List[i].mass), sizeof(float));
            file.read((char *)(&List[i].pid), sizeof(long));
            file.read((char *)(&List[i].cpuid), sizeof(int));
        }
    } else {
        cout << "Failed to open " << filename << endl;
    }
    file.close();
    return 0;
}

float ParticleList::ScaleHeight()
{
    float Hp = 0;
    for (long i = 0; i < n; i++) {
        Hp += List[i].x[2]*List[i].x[2];
    }
    Hp = sqrt(Hp/n);
    return Hp;
}

ParticleList::~ParticleList()
{
    delete List;
}

