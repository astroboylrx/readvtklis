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

/*! \class TreeNode
 *  \brief class for tree node */
template<int D>
class TreeNode {
private:
    
public:
    int level;
    TreeNode<D> *father;
    vector<TreeNode<D> *> daughter;
    
    float L;
    float center[D];
    int pid;
    int np;
    vector<int> deep_pids;
    
    /*! \fn TreeNode()
     *  \brief constructor */
    TreeNode() {
        level = 0;
        father = nullptr;
        daughter.resize(2<<D);
        for (typename vector<TreeNode<D> *>::iterator it = daughter.begin(); it != daughter.end(); it++) {
            *it = nullptr;
        }
        L = 0.0;
        for (int i = 0; i != D; i++) {
            center[i] = 0.0;
        }
        pid = -1;
        np = 0;
    }
    
    /*! \fn ~TreeNode()
     *  \brief destructor */
    ~TreeNode() {
        if (daughter.size() > 0) {
            vector<TreeNode<D> *> temp;
            daughter.swap(temp);
        }
        if (deep_pids.size() > 0) {
            vector<int> temp;
            deep_pids.swap(temp);
        }
    }
};

/*! \class Tree
 *  \brief tree class */
template <int D>
class Tree {
private:
    
public:
    int level;
    int level_limit;
    
    TreeNode<D> *root;
    VtkFile *vf;
    ParticleList *pl;
    
    float size_to_diagonal;
    float *half_diagonal;
    float *Radius;
    float **MaxD;
    
    float m1par;
    float **RMPL;
    float shortestL;
    float tempcenter[D];
    
    /*! \fn Tree(int n_file)
     *  \brief constructor */
    Tree(int n_file) {
        size_to_diagonal = sqrt(D)/2.0;
        MaxD = nullptr;
        RMPL = new float*[n_file];
        for (int i = 0; i != n_file; i++) {
            RMPL[i] = nullptr;
        }
        Radius = nullptr;
        half_diagonal = nullptr;
        Initialize();
    }
    
    /*! \fn void Initialize()
     *  \brief initialization */
    void Initialize() {
        root = new TreeNode<D>;
        shortestL = 1e5;
    }
    
    /*! \fn void CleanMem(TreeNode<D> *p)
     *  \brief clean memory */
    void CleanMem(TreeNode<D> *p) {
        for (typename vector<TreeNode<D> *>::iterator it = p->daughter.begin(); it != p->daughter.end(); it++) {
            if (*it != nullptr) {
                CleanMem(*it);
            }
            delete (*it);
            (*it) = nullptr;
        }
    }
    
    /*! \fn ~Tree()
     *  \brief destructor */
    ~Tree() {
        CleanMem(root);
        delete root; root = nullptr;
        for (int i = 0; i != fio->n_file; i++) {
            if (RMPL[i] != nullptr) {
                delete [] RMPL[i];
            }
        }
        if (RMPL != nullptr) {
            delete [] RMPL;
        }
        if (Radius != nullptr) {
            delete [] Radius;
        }
        if (half_diagonal != nullptr) {
            delete [] half_diagonal;
        }
    }
    
    /*! \fn int GetDivision(float center[], float pos[])
     *  \brief calculate which daughter the particle belongs to */
    int GetDivision(float center[], float pos[]) {
        int i = 0;
        for (int xyz = 0; xyz != D; xyz++) {
            i += (pos[xyz] > center[xyz])<<xyz;
        }
        return i;
    }
    
    /*! \fn TreeNode<D> *CreateNode(TreeNode<D> *p, int &division)
     *  \brief Create a tree node and return the pointer */
    TreeNode<D> *CreateNode(TreeNode<D> *p, int &division) {
        TreeNode<D> *q = new TreeNode<D>;
        q->level = p->level + 1;
        q->father = p;
        q->L = p->L / 2.0;
        for (int xyz = 0; xyz != D; xyz++) {
            q->center[xyz] = p->center[xyz] + q->L * (((division>>xyz) & 1) - 0.5);
        }
        p->daughter[division] = q;
        return q;
    }
    
    /*! \fn void AddParticle(TreeNode<D> *p, Particle &it)
     *  \brief add a particle into the tree structure */
    void AddParticle(TreeNode<D> *p, Particle &it) {
        p->np++; // anyway we need to do this
        
        // if it is an empty leaf previously, just put particle id int
        if (p->np == 1) {
            p->pid = static_cast<int>(it.pid);
            return;
        }
        // if it touch the level limit, just put particle id into vector
        if (p->level == level_limit) {
            p->deep_pids.push_back(static_cast<int>(it.pid));
            return;
        }
        
        int di;
        // Or else it has one already and level < level limit, move it to daughter first
        if (p->np == 2) {
            di = GetDivision(p->center, pl->List[p->pid].x);
            AddParticle(CreateNode(p, di), pl->List[p->pid]);
            p->pid = -1; // means we have particles in daughter nodes
        }
        
        // Deal with new particle now
        di = GetDivision(p->center, it.x);
        if (p->daughter[di] == nullptr) {
            CreateNode(p, di);
        }
        AddParticle(p->daughter[di], it);
    }
    
    /*! \fn void BuildTree(VtkFile *VF, ParticleList *PL, int file_i)
     *  \brief build the tree by add particles into it */
    void BuildTree(VtkFile *VF, ParticleList *PL, int file_i) {
        
        CleanMem(root);
        delete root; root = nullptr;
        Initialize();
        vf = VF;
        pl = PL;
        
        // first assign root node and calculate m1par
        int rootNx = 0;
        for (int i = 0; i != D; i++) {
            if (rootNx < vf->dimensions[i]) {
                rootNx = vf->dimensions[i];
            }
            if (shortestL > vf->L[i]) {
                shortestL = vf->L[i];
            }
        }
        
        root->L = rootNx * vf->spacing[0]; // assume spacing is the same
        m1par = 0.02*sqrt(2*PI)*vf->L[0]*vf->L[1]/pl->n; // calculate mass first, then surface density
        
        /***** The following need the "level" to allocate space *****/
        level = round(log10(rootNx)/log10(2.0));
        level_limit = round(log10(root->L/0.0001)/log10(2.0)); // no point to go further down to <5e-4
        
        if (MaxD == nullptr) {
            MaxD = new float*[level+1];
            if (Radius == nullptr) {
                Radius = new float[level+1];
            }
            for (int i = 0; i <= level; i++) {
                MaxD[i] = new float[D];
                Radius[i] = vf->spacing[0] * pow(2.0, i-1); // used in later func
                for (int j = 0; j != D; j++) {
                    MaxD[i][j] = vf->L[j] - Radius[i];
                }
#ifdef OutflowRate
                if (D == 3) {
                    // let it bigger than any distance so no periodic computation happen in x3 direction
                    // in fact all particles are in the middle plane, and the RMPL[file_i][root_level] is calculated by another way. This is not main concern
                    MaxD[i][2] = vf->L[2] * 2;
                }
#endif
                //if (myMPI->myrank == 0) {
                //    cout << myMPI->prank() << "D = " << D << ", MaxD[" << i << "]=" << pvector<D, float>(MaxD[i]) << ", R[i]=" << Radius[i] << endl;
                //}
            }
        }
        if (RMPL[file_i] == nullptr) {
            RMPL[file_i] = new float[level+1];
            for (int i = 0; i <= level; i++) {
                RMPL[file_i][i] = 0;
            }
        }
        if (half_diagonal == nullptr) {
            half_diagonal = new float[level_limit+1];
            half_diagonal[0] = root->L * size_to_diagonal;
            for (int i = 1; i <= level_limit; i++) {
                half_diagonal[i] = half_diagonal[i-1] / 2.0;
            }
        }
        
        /***** compute particle numbers for each node *****/
        for (vector<Particle>::iterator it = pl->List.begin(); it != pl->List.end(); ++it) {
            AddParticle(root, *it);
        } //*/

    }
    
    /*! \fn Distance(float x[], float y[], int &i_MaxD)
     *  \brief calculate distance between two points, considering the periodic BCs */
    float Distance(float x[], float y[], int &i_MaxD) {
        float d = 0, dr[D];
        for (int i = 0; i != D; i++) {
            dr[i] = fabs(x[i] - y[i]);
            if (dr[i] > MaxD[i_MaxD][i]) {
                dr[i] = vf->L[i] - dr[i];
            }
            d += dr[i] * dr[i];
        }
        return sqrt(d);
    }
    
    /*! \fn void EvaluateAccurateRadius(TreeNode<D> *p, float x[], int &npar, float &R, int &i_MaxD)
     *  \brief now calculate npar for an accurate general sphere in D dimensions at x with r=Radius */
    void EvaluateAccurateRadius(TreeNode<D> *p, float x[], int &npar, float &R, int &i_MaxD) {
        float d = Distance(p->center, x, i_MaxD);
        if (d > R + size_to_diagonal * p->L) { // half_diagonal[p->level]
            // completely outside radius
            ;
        } else if (d < R - size_to_diagonal * p->L) { // half_diagonal[p->level]
            // completely inside radius
            npar += p->np;
        } else {
            // intersect
            // if it is leaf node
            if (p->np == 1 || p->level == level_limit) {
                if (Distance(x, pl->List[p->pid].x, i_MaxD) <= R) {
                    npar++;
                }
                for (int i = 0; i < p->deep_pids.size(); i++) {
                    if (Distance(x, pl->List[p->deep_pids[i]].x, i_MaxD) <= R) {
                        npar++;
                    }
                }
            } else {
                for (int i = 0; i != 2<<D; i++) {
                    if (p->daughter[i] != nullptr) {
                        EvaluateAccurateRadius(p->daughter[i], x, npar, R, i_MaxD);
                    }
                }
            }
        }
    }
    
    /*! \fn void RhopMaxPerLevel(int file_i)
     *  \brief find the max rhop/sigamp/... within a sphere with radius of N*dx */
    void RhopMaxPerLevel(int file_i) {
#ifdef ENABLE_MPI
        int Npoints, **indices, index;
        
        if (D == 1) {
            Npoints = vf->dimensions[0];
            indices = new int*[Npoints];
            for (int npo = 0; npo != Npoints; npo++) {
                indices[npo] = new int[1];
                indices[npo][0] = npo % vf->dimensions[0];
            }
        }
        if (D == 2) {
            Npoints = vf->dimensions[0] * vf->dimensions[1];
            indices = new int*[Npoints];
            
            // list all the indices of points we need to check
            for (int npo = 0; npo != Npoints; npo++) {
                indices[npo] = new int[2];
                indices[npo][1] = npo / vf->dimensions[0];
                indices[npo][0] = npo % vf->dimensions[0];
            }
        }
        if (D == 3) {
            // certainly those outside the center zone (with rhop = 0) can be omitted
            int ks = -1, ke = -1, iz, iy, ix;
            for (iz = vf->dimensions[2]/2-1; iz != 0; iz--) {
                ks = -1;
                for (iy = 0; iy != vf->dimensions[1]; iy++) {
                    for (ix = 0; ix != vf->dimensions[0]; ix++) {
                        if (vf->cd_scalar[1].data[iz][iy][ix] > 0) {
                            ks = -2; break;
                        }
                    }
                    if (ks == -2) {
                        break;
                    }
                }
                if (ks == -1) {
                    ks = iz; break;
                }
            }
            if (ks == -2) {
                ks = 0;
            }
            ks += 2; // save some time
            for (iz = vf->dimensions[2]/2; iz != vf->dimensions[2]; iz++) {
                ke = -1;
                for (iy = 0; iy != vf->dimensions[1]; iy++) {
                    for (ix = 0; ix != vf->dimensions[0]; ix++) {
                        if (vf->cd_scalar[1].data[iz][iy][ix] > 0) {
                            ke = -2; break;
                        }
                    }
                    if (ke == -2) {
                        break;
                    }
                }
                if (ke == -1) {
                    ke = iz; break;
                }
            }
            if (ke == -2) {
                ke = vf->dimensions[2]-1;
            }
            ke -= 2; // save some time
            
            // list all the points we need to check
            Npoints = vf->dimensions[0] * vf->dimensions[1] * (ke-ks+1);
            indices = new int*[Npoints];
            int temp = vf->dimensions[0] * vf->dimensions[1];
            for (int npo = 0; npo != Npoints; npo++) {
                indices[npo] = new int[3];
                indices[npo][2] = npo / temp + ks;
                indices[npo][1] = (npo % temp) / vf->dimensions[0];
                indices[npo][0] = (npo % temp) % vf->dimensions[0];
            }
        }
        
        // Begin to process data
        int mine_contribution = 0;
        for (index = myMPI->myrank; index < Npoints; index += myMPI->numprocs) {
            int npar;

            // for increasing sample when it comes to smallest sphere
            float *temp_cellcenter;
            short selection, xyz, yesorno;
            
            if (D == 3) {
                temp_cellcenter = vf->cell_center[indices[index][2]][indices[index][1]][indices[index][0]];
            } else if (D == 2) {
                temp_cellcenter = vf->cell_center[0][indices[index][1]][indices[index][0]];
            } else if (D == 1) {
                temp_cellcenter = vf->cell_center[0][0][indices[index][0]];
            }
            
            for (int i = 0; i <= level; i++) {
                if (Radius[i] >= shortestL) {
                    break; // further step make no sense, so save some computation time
                }
                npar = 0;
                
                // Count real particle numbers in accurate sphere
                EvaluateAccurateRadius(root, temp_cellcenter, npar, Radius[i], i);
                
                if (npar > RMPL[file_i][i]) {
                    RMPL[file_i][i] = npar;
                }
                
                // some small circles need more sampling
                if (i < 3) {
                    for (selection = 1; selection != 2<<D; selection++) {
                        npar = 0;
                        for (xyz = 0; xyz != D; xyz++) {
                            yesorno = (selection >> xyz) & 1;
                            tempcenter[xyz] = temp_cellcenter[xyz] - yesorno * Radius[0];
                        }
                        EvaluateAccurateRadius(root, tempcenter, npar, Radius[i], i);
                        if (npar > RMPL[file_i][i]) {
                            RMPL[file_i][i] = npar;
                        }
                        if (i == 0) {
                            float tempcenter2[3];
                            short selection2, xyz2, yesorno2;
                            for (selection2 = 1; selection2 != 2<<D; selection2++) {
                                npar = 0;
                                for (xyz2 = 0; xyz2 != D; xyz2++) {
                                    yesorno2 = (selection2 >> xyz2) & 1;
                                    tempcenter2[xyz2] = tempcenter[xyz2] - yesorno2 * Radius[0]/2.0;
                                }
                                EvaluateAccurateRadius(root, tempcenter2, npar, Radius[i], i);
                                if (npar > RMPL[file_i][i]) {
                                    RMPL[file_i][i] = npar;
                                }
                            }
                        }
                    }
                }
            }
            mine_contribution++;
        }
        cout << myMPI->prank() << mine_contribution << " jobs are done." << endl;

    }
#endif /* ENABLE_MPI */
    
};


////////////// The following are deprecated!


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
