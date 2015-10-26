//
//  octree.hpp
//  readvtklis
//
//  Created by Rixin Li on 10/12/15.
//  Copyright © 2015 Rixin Li. All rights reserved.
//

#ifndef octree_h
#define octree_h

#include "global.h"
#include "fop.h"
#include "readvtk.h"
#include "readlis.h"
const float PI = 3.1415927;

/*! \class OctreeNode
 *  \brief the node class of Octree
 */
class OctreeNode {
private:
    
public:
    int level;                                      /*!< the level of this node */
    OctreeNode *Father;                             /*!< father node pointer of this node */
    vector<OctreeNode *> Daughter;                  /*!< eight daughters in the increasing order of x, y and z */
    float rhop;                                     /*!< the sum of particle density for all the cells in this node */
    int Nx;                                         /*!< the number of cells in 1 direction,  forced to be the same in each directions */
    float Lx;                                       /*!< length of this node */
    float center[3];                                /*!< the coordinates of the center of this node */
    long np;                                        /*!< number of particles in this node */
    
    /*! \fn OctreeNode()
     *  \brief Constructor of tree node */
    OctreeNode();
    
    /*! \fn ~OctreeNode()
     *  \brief Destructor */
    ~OctreeNode();
    
};

/*! \class Octree
 *  \brief establish octree to calculate weighted particle density
 */
class Octree {
private:
    
public:
    int level;                                      /*!< total levels of the tree */
    long *NxCubic;                                  /*!< Nx^3 at each level */
    OctreeNode *root;                               /*!< the root of the whole tree */
    VtkFile *vf;                                    /*!< the VtkFile pointer used to build this tree */
    ParticleList *pl;                               /*!< the ParticleList pointer used to build this tree */
    
    float Max_Rhop;                                 /*!< maximum density of particle */
    float RpAV;                                     /*!< <rho_p> */
    float RpSQ;                                     /*!< <rho_p^2>^0.5 */
    float RpQU;                                     /*!< <rho_p^4>^0.25 */
    
    enum { etar = 16 };                             /*!< etar in # of cells */
    float s3o2;                                     /*!< sqrt(3)/2 */
    float foPio3;                                   /*!< 4*pi/3 */
    float Radius;                                   /*!< half etar */
    float **MaxD;                                   /*!< max allowed distance that we do not consider periodic boundary */
    
    float *RMPL;                                    /*!< Rhop_Max Per Level, in reverse order */
    
    /*! \fn Octree();
     *  \brief Constructor of tree sturcture */
    Octree();
    
    /*! \fn Initialize();
     *  \brief refresh the class */
    int Initialize();
    
    /*! \fn ~Octree()
     *  \brief Destructor */
    ~Octree();
    
    /*! \fn void CleanMem(OctreeNode *p);
     *  \brief clean all the memory */
    void CleanMem(OctreeNode *p);
    
    /*! \fn int BuildTree(VtkFile *VF, ParticleList *PL)
     *  \brief Build the whole tree */
    int BuildTree(VtkFile *VF, ParticleList *PL);
    
    /*! \fn OctreeNode *AddCell(float cc[3], float rhop)
     *  \brief using cell center to find the tree node, create nodes if needed */
    OctreeNode *AddCell(float cc[3], float rhop);
    
    /*! \fn int GetOctant(float center[3], T pos[3]);
     *  \brief return the index of the daughter */
    template<typename T>
    int GetOctant(float center[3], T pos[3]);
    
    /*! \fn int AddParticle(float x[3])
     *  \brief assign particle to one node */
    int AddParticle(float x[3]);
    
    /*! \fn void EvaluateOnePoint(OctreeNode *p, T x[3], float *rp, long *cells)
     *  \brief calculate rhop for a sphere centered at x with r=Radius */
    template<typename T>
    void EvaluateOnePoint(OctreeNode *p, T x[3], float *rp, long *cells);
    
    /*! \fn float Distance(OctreeNode *p, float x[3])
     *  \brief calculate distance between cell center and particle */
    template<typename T>
    float Distance(OctreeNode *p, T x[3]);
    
    /*! \fn void RhopMaxPerLevel()
     *  \brief Find the max rhop within a sphere with radius of N*dx */
    void RhopMaxPerLevel();
};




#endif /* octree_h */
