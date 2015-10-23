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

/*! \class OctreeNode
 *  \brief the node class of Octree
 */
class OctreeNode {
private:
    
public:
    int level;                                      /*!< the level of this node */
    OctreeNode *Father;                             /*!< father node pointer of this node */
    vector<OctreeNode *> Daughter;                  /*!< eight daughters in the increasing order of x, y and z */
    double rhop;                                    /*!< the sum of particle density for all the cells in this node */
    int Nx;                                         /*!< the number of cells,  forced to be the same in each directions */
    double Lx;                                      /*!< length of this node */
    double center[3];                               /*!< the coordinates of the center of this node */
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
    OctreeNode *root;                               /*!< the root of the whole tree */
    VtkFile *vf;                                    /*!< the VtkFile pointer used to build this tree */
    ParticleList *pl;                               /*!< the ParticleList pointer used to build this tree */
    
    double Max_Rhop;                                /*!< maximum density of particle */
    double RpAV;                                    /*!< <rho_p> */
    double RpSQ;                                    /*!< <rho_p^2>^0.5 */
    double RpQU;                                    /*!< <rho_p^4>^0.25 */
    
    enum { etar = 16 };                             /*!< etar in # of cells */
    double s3o2;                                    /*!< sqrt(3)/2 */
    double Radius;                                  /*!< half etar */
    double MaxD[3];                                 /*!< max allowed distance that we do not consider periodic boundary */
    double RpEtar;                                  /*!< weighted particle density computed by summing all the related cells around each particle, using eta*r as the radius; define "related" by measuing the distance between cell center and particle */
    
    double *maxrhop;
    
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
    
    /*! \fn OctreeNode *AddCell(double cc[3], double rhop)
     *  \brief using cell center to find the tree node, create nodes if needed */
    OctreeNode *AddCell(double cc[3], double rhop);
    
    /*! \fn int GetOctant(double center[3], T pos[3]);
     *  \brief return the index of the daughter */
    template<typename T>
    int GetOctant(double center[3], T pos[3]);
    
    /*! \fn int AddParticle(float x[3])
     *  \brief assign particle to one node */
    int AddParticle(float x[3]);
    
    /*! \fn void EvaluateOnePoint(OctreeNode *p, T x[3], double *rp)
     *  \brief calculate rhop for a sphere centered at x with r=Radius */
    template<typename T>
    void EvaluateOnePoint(OctreeNode *p, T x[3], double *rp);
    
    /*! \fn double Distance(OctreeNode *p, float x[3])
     *  \brief calculate distance between cell center and particle */
    template<typename T>
    double Distance(OctreeNode *p, T x[3]);
};




#endif /* octree_h */
