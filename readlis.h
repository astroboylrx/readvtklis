//
//  readlis.h
//  readvtklis
//
//  Created by Rixin Li on 1/14/15.
//  Copyright (c) 2015 Rixin Li. All rights reserved.
//

#ifndef __readvtklis__readlis__
#define __readvtklis__readlis__

#include "global.h"
#include "fop.h"

/*! \class Particle
 *  \brief Information about one single particle
 */
class Particle {
private:
    
public:
    int cpuid;                      /*!< the id of cpu (which contains it) */
    long pid;                       /*!< particle ID in list */
    float x[3], v[3], rad, mass;    /*!< partile infomation */
};

/*! \class ParticleList
 *  \brief Information about the entire particle list
 */
class ParticleList {
private:
    
public:
    int ntype;                      /*!< number of par types */
    long n, npar_ghost;             /*!< number of particles */
    float coorlim[12];              /*!< the coordinate limit */
    float *typeinfo;                /*!< info of par type */
    float time;                     /*!< current time */
    float dt;                       /*!< the length of time step */
    float x1l, x1u;                 /*!< lower/upper x1 range for data in this grid */
    float x2l, x2u;                 /*!< lower/upper x2 range for data in this grid  */
    float x3l, x3u;                 /*!< lower/upper x3 range for data in this grid  */
    float x1dl, x1du;               /*!< lower/upper x1 range for data in this domain */
    float x2dl, x2du;               /*!< lower/upper x2 range for data in this domain  */
    float x3dl, x3du;               /*!< lower/upper x3 range for data in this domain  */
    vector <Particle> List;         /*!< particle list vector */
    
    ParticleList();                 /*!< constructor */
    ~ParticleList();                /*!< destructor */
    
    /*! \fn int ReadLis(string filename)
     *  \brief Read particle list from file */
    int ReadLis(string filename);
    
    /*! \fn int ScaleHeight(double &Hp, double &Hp_in1sigma);
     *  \brief calculate the scale height of partiles */
    int ScaleHeight(double &Hp, double &Hp_in1sigma);
    
    /*! \fn int InitializeList()
     *  \brief free List to control memory */
    int InitializeList();
    
    /*! \fn int PrintInfo()
     *  \brief print basic info */
    int PrintInfo();
    
    /*! \fn int GetNumprocs()
     *  \brief get the number of processors */
    int GetNumprocs();

};


#endif /* defined(__readvtklis__readlis__) */
