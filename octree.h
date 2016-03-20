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
const double PI = 3.141592653589793;

/*! \class QuadtreeNode
 *  \brief the node class of Quadtree */
class QuadtreeNode {
private:
    
public:
    int level;                          /*!< level of this node */
    QuadtreeNode *Father;               /*!< father node pointer */
    vector<QuadtreeNode *> Daughter;    /*!< four daughters in the increasing order of x, y */
    
    float Lx;                           /*!< length of this node */
    float center[2];                    /*!< the coordinates of the center of this node */
    long pid;                           /*!< particle index in this node */
    long np;                            /*!< how many particles in this node */
    vector<long> deep_pids;             /*!< keep more paritcles in this node if exceeding a certain level (particles do overlap since we only account for x and y) */
    
    /*! \fn QuadtreeNode()
     *  \brief constructor */
    QuadtreeNode();
    
    /*! \fn QuadtreeNode()
     *  \brief destructor */
    ~QuadtreeNode();
    
};

/*! \class Quadtree
 *  \brief establish quadtree to calculate particle surface density per level */
class Quadtree {
private:
    
public:
    int level;                          /*!< total levels of the tree down to cell size */
    int level_limit;                    /*!< the deeper limit of tree levels (to put particle in deep_pids) */
    
    QuadtreeNode *root;                 /*!< the root of the whole tree */
    VtkFile *vf;                        /*!< the VtkFile pointer used to build this tree */
    ParticleList *pl;                   /*!< the ParticleList pointer used to build this tree */
    
    float s2o2;                         /*!< sqrt(2)/2 */
    float *Radius;                      /*!< radius of the circle for surface density calculation */
    float **MaxD;                       /*!< max allowed distance that we do not consider periodic boundary */
    
    float m1par;                        /*!< msaa of 1 particle */
    float **SMPL;                       /*!< Sigmap_Max Per Level, in reverse order for all files */
    float shortestL;                    /*!< the shortest side length of numerical domain */
    float tempcenter[2];                /*!< temp center for advanced search on small radius */
    
    /*! \fn Quadtree()
     *  \brief constructor of tree structure */
    Quadtree(int n_file);
    
    /*! \fn Initialize()
     *  \brief refresh the class */
    int Initialize();
    
    /*! \fn ~Quadtree()
     *  \brief destructor */
    ~Quadtree();
    
    /*! \fn void CleanMem(QuadtreeNode *p)
     *  \brief clean all the memory */
    void CleanMem(QuadtreeNode *p);
    
    /*! \fn int BuildTree(VtkFile *VF, ParticleList *PL, int file_i)
     *  \brief build the whole tree */
    int BuildTree(VtkFile *VF, ParticleList *PL, int file_i);

    
    /*! \fn int GetQuadrant(float center[], T pos[]);
     *  \brief return the index of the daughter */
    template<typename T>
    int GetQuadrant(float center[], T pos[]);
    
    /*! \fn QuadtreeNode *CreateNode(QuadtreeNode *p, int &quadrant)
     *  \brief create new node by the father node and quadrant */
    QuadtreeNode *CreateNode(QuadtreeNode *p, int &quadrant);
    
    /*! \fn int AddParticle(QuadtreeNode *p, Particle &it)
     *  \brief assign particle's index to one node */
    int AddParticle(QuadtreeNode *p, Particle &it);
    
    /*! \fn void EvaluateAccurateSphere(QuadtreeNode *p, T x[], long &npar, float &R, int &i_MaxD)
     *  \brief calculate npar for an accurate sphere centered at x with r=Radius */
    template<typename T>
    void EvaluateAccurateSphere(QuadtreeNode *p, T x[], long &npar, float &R, int &i_MaxD);
    
    /*! \fn float Distance(T x[], T y[], int &i_MaxD)
     *  \brief calculate distance between two locations, i_MaxD is index of MaxD */
    template<typename T>
    float Distance(T x[], T y[], int &i_MaxD);
    
    /*! \fn void SigmapMaxPerLevel(int file_i)
     *  \brief Find the max rhop within a sphere with radius of N*dx */
    void SigmapMaxPerLevel(int file_i);
};

/*! \class OctreeNode
 *  \brief the node class of Octree
 */
class OctreeNode {
private:
    
public:
    int level;                          /*!< the level of this node */
    OctreeNode *Father;                 /*!< father node pointer of this node */
    vector<OctreeNode *> Daughter;      /*!< eight daughters in the increasing order of x, y and z */

    float Lx;                           /*!< length of this node */
    float center[3];                    /*!< the coordinates of the center of this node */
    long pid;                           /*!< particle index in this node */
    long np;                            /*!< number of particles in this node */
    vector<long> deep_pids;             /*!< keep more particles in this node if exceeding a certain level */
    
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
    int level;                          /*!< total levels of the tree down to cell size */
    int level_limit;                    /*!< the deeper limit of tree levels (to put particle in deep_pids */
    
    OctreeNode *root;                   /*!< the root of the whole tree */
    VtkFile *vf;                        /*!< the VtkFile pointer used to build this tree */
    ParticleList *pl;                   /*!< the ParticleList pointer used to build this tree */
    
    float Max_Rhop;                     /*!< maximum density of particle */
    double RpAV;                        /*!< <rho_p> */
    double RpSQ;                        /*!< <rho_p^2>^0.5 */
    double RpQU;                        /*!< <rho_p^4>^0.25 */
    
    enum { etar = 16 };                 /*!< etar in # of cells */
    float s3o2;                         /*!< sqrt(3)/2 */
    float foPio3;                       /*!< 4*pi/3 */
    float *Radius;                      /*!< half etar */
    float **MaxD;                       /*!< max allowed distance that we do not consider periodic boundary */
    float m1par;                        /*!< mass of 1 particle */
    float **RMPL;                       /*!< Rhop_Max Per Level, in reverse order for all files */
    float shortestL;                    /*!< The shortest side length of the numerical domain */
    float tempcenter[3];
    
    /*! \fn Octree(int n_file);
     *  \brief Constructor of tree sturcture */
    Octree(int n_file);
    
    /*! \fn Initialize();
     *  \brief refresh the class */
    int Initialize();
    
    /*! \fn ~Octree()
     *  \brief Destructor */
    ~Octree();
    
    /*! \fn void CleanMem(OctreeNode *p);
     *  \brief clean all the memory */
    void CleanMem(OctreeNode *p);
    
    /*! \fn int BuildTree(VtkFile *VF, ParticleList *PL, int file_i)
     *  \brief Build the whole tree */
    int BuildTree(VtkFile *VF, ParticleList *PL, int file_i);
    
    /*! \fn OctreeNode *CreateNode(OctreeNode *p, int &octant)
     *  \brief create new node depend on father node and octant */
    OctreeNode *CreateNode(OctreeNode *p, int &octant);
    
    /*! \fn int GetOctant(float center[3], T pos[3]);
     *  \brief return the index of the daughter */
    template<typename T>
    int GetOctant(float center[3], T pos[3]);
    
    /*! \fn int AddParticle(OctreeNode *p, Particle &it)
     *  \brief assign particle to one node */
    int AddParticle(OctreeNode *p, Particle &it);
    
    /*! \fn void EvaluateOctreeSphere(OctreeNode *p, T x[3], double &rp, long &cells, float &R)
     *  \brief calculate rhop for an octree-approximate sphere centered at x with r=Radius */
    template<typename T>
    void EvaluateOctreeSphere(OctreeNode *p, T x[3], double &rp, long &cells, float &R);
    
    /*! \fn void EvaluateAccurateSphere(OctreeNode *p, T x[3], long &npar, float &R, int &i_MaxD)
     *  \brief calculate rhop for an accurate sphere centered at x with r=Radius */
    template<typename T>
    void EvaluateAccurateSphere(OctreeNode *p, T x[3], long &npar, float &R, int &i_MaxD);
    
    /*! \fn float Distance(T x[3], T y[3], int &i_MaxD)
     *  \brief calculate distance between two locations, i_MaxD is index of MaxD */
    template<typename T>
    float Distance(T x[3], T y[3], int &i_MaxD);
    
    /*! \fn void RhopMaxPerLevel(int file_i)
     *  \brief Find the max rhop within a sphere with radius of N*dx */
    void RhopMaxPerLevel(int file_i);
    
};




#endif /* octree_h */
