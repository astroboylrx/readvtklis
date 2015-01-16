//
//  readlis.h
//  readvtklis
//
//  Created by Rixin Li on 1/14/15.
//  Copyright (c) 2015 Rixin Li. All rights reserved.
//

#ifndef __readvtklis__readlis__
#define __readvtklis__readlis__

#include "fop.h"
/********** Particle class **********/
class Particle {
private:
    
public:
    int cpuid;
    long pid; // particle ID in list
    float x[3], v[3], rad, mass; // partile infomation
};

/********** ParticleList class **********/
class ParticleList {
private:
    
public:
    // the number of types
    int ntype;
    // the number of particles
    long n;
    // the coordinate limit, info of partile type
    // the time, and time step
    float coorlim[12], *typeinfo, time, dt;
    // particle list vector
    vector <Particle> List;
    // constructor and destructor
    ParticleList();
    ~ParticleList();
    
    // Read particle list from file
    int ReadLis(string filename);
    // calculate the scale height of partiles
    float ScaleHeight();
    // free List to control memory
    int InitializeList();
};


#endif /* defined(__readvtklis__readlis__) */
