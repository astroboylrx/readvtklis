//
//  octree.cpp
//  readvtklis
//
//  Created by Rixin Li on 10/12/15.
//  Copyright © 2015 Rixin Li. All rights reserved.
//

#include "octree.h"

/**********************************/
/********** QuadTreeNode **********/
/**********************************/

/********** Constructor **********/
QuadtreeNode::QuadtreeNode()
{
    level = 0;
    Father = NULL;
    Daughter.resize(4);
    for (vector<QuadtreeNode *>::iterator it = Daughter.begin(); it != Daughter.end(); it++) {
        (*it) = NULL;
    }
    Lx = 0.0;
    center[0] = 0.0;
    center[0] = 0.0;
    pid = -1;
    np = 0;
}

/********** Destructor **********/
QuadtreeNode::~QuadtreeNode()
{
    if (Daughter.size() > 0) {
        vector<QuadtreeNode *> temp;
        Daughter.swap(temp);
    }
}


/******************************/
/********** Quadtree **********/
/******************************/

/********** Constructor **********/
Quadtree::Quadtree(int n_file)
{
    s2o2 = sqrt(2.0)/2.0;
    MaxD = NULL;
    SMPL = new float*[n_file];
    for (int i = 0; i != n_file; i++) {
        SMPL[i] = NULL;
    }
    Radius = NULL;
    Initialize();
}

/********** Initialize **********/
int Quadtree::Initialize()
{
    root = new QuadtreeNode;
    shortestL = 1e10;
    return 0;
}

/********** Destructor **********/
Quadtree::~Quadtree()
{
    CleanMem(root);
    delete root; root = NULL;
    for (int i = 0; i != fio->n_file; i++) {
        if (SMPL[i] != NULL) {
            delete [] SMPL[i];
        }
    }
    if (SMPL != NULL) delete [] SMPL;
    if (Radius != NULL) delete [] Radius;
}

/********** CleanMem **********/
void Quadtree::CleanMem(QuadtreeNode *p)
{
    for (vector<QuadtreeNode *>::iterator it = p->Daughter.begin(); it != p->Daughter.end(); it++) {
        if ((*it) != NULL) {
            CleanMem((*it));
            delete (*it);
            (*it) = NULL;
        }
    }
}


/********** BuildTree **********/
/*! \fn BuildTree(VtkFile *VF, ParticleList *PL, int file_i)
 *  \brief Build the whole tree */
int Quadtree::BuildTree(VtkFile *VF, ParticleList *PL, int file_i)
{
    CleanMem(root);
    delete root; root = NULL;
    Initialize();
    vf = VF;
    pl = PL;
    
    /***** first assign root node and calculate m1par *****/
    int rootNx = 0;
    for (int i = 0; i != 2; i++) {
        if (rootNx < vf->dimensions[i]) {
            rootNx = vf->dimensions[i];
        }
        if (shortestL > vf->L[i]) {
            shortestL = vf->L[i];
        }
    } // Nx and Ny should be the same
    
    for (int i = 0; i != 2; i++) {
        root->center[i] = 0.0; // this code only considers special cases
    }
    root->Lx = rootNx * vf->spacing[0]; // assuming spacing are all the same
    m1par = 0.02*sqrt(2*PI)*vf->L[0]*vf->L[1]/pl->n; // calculate mass first, then surface density
    
    
    /***** The following need the "level" to allocate space *****/
    level = round(log10(rootNx)/log10(2.0));
    level_limit = round(log10(root->Lx/0.0005)/log10(2.0)); // no point to go further down to <5e-4
    
    if (MaxD == NULL) {
        MaxD = new float*[level+1];
        if (Radius == NULL) {
            Radius = new float[level+1];
        }
        for (int i = 0; i <= level; i++) {
            MaxD[i] = new float[2];
            Radius[i] = vf->spacing[0] * pow(2.0, i-1); // used in later func
            MaxD[i][0] = vf->L[0] - Radius[i];
            MaxD[i][1] = vf->L[1] - Radius[i];
            //cout << myMPI->prank() << "MaxD[" << i << "]=" << pvector<float>(MaxD[i]) << ", R[i]=" << Radius[i] << endl;
        }
    }
    if (SMPL[file_i] == NULL) {
        SMPL[file_i] = new float[level+1];
        for (int i = 0; i <= level; i++) {
            SMPL[file_i][i] = 0;
        }
    }
    
    /***** compute particle numbers for each node *****/
    for (vector<Particle>::iterator it = pl->List.begin(); it != pl->List.end(); ++it) {
        AddParticle(root, *it);
    } //*/
    
    return 0;
}

/********** GetOctant **********/
/*! \fn int GetQuadrant(float center[], T pos[]);
 *  \brief return the index of the daughter */
template<typename T>
int Quadtree::GetQuadrant(float center[], T pos[])
{
    int i = 0;
    i += (pos[0]>center[0]);     // x
    i += (pos[1]>center[1])<<1;  // y
    return i;
}

/********** CreateNode **********/
/*! \fn QuadtreeNode *CreateNode(QuadtreeNode *p, int &quadrant)
 *  \brief create new node depend on father node and quadrant */
QuadtreeNode * Quadtree::CreateNode(QuadtreeNode *p, int &quadrant)
{
    QuadtreeNode *q = new QuadtreeNode;
    q->level = p->level + 1;
    q->Father = p;
    q->Lx = p->Lx / 2;
    q->center[0] = p->center[0] + q->Lx * ((quadrant & 1) - 0.5);
    q->center[1] = p->center[1] + q->Lx * ((quadrant > 1) - 0.5);
    p->Daughter[quadrant] = q;
    return q;
}

/********** AddParticle **********/
/*! \fn int AddParticle(QuadtreeNode *p, Particle &it)
 *  \brief assign particle to one node */
int Quadtree::AddParticle(QuadtreeNode *p, Particle &it)
{
    p->np++; // anyway we need to do this
    
    // if it is an empty leaf previously, just put particle id int
    if (p->np == 1) {
        p->pid = it.pid;
        return 0;
    }
    // if it touch the level limit, just put particle id into vector
    if (p->level == level_limit) {
        p->deep_pids.push_back(it.pid);
        return 0;
    }
    
    int di;
    // Or else it has one already and level < level limit, move it to daughter first
    if (p->np == 2) {
        di = GetQuadrant<float>(p->center, pl->List[p->pid].x);
        AddParticle(CreateNode(p, di), pl->List[p->pid]);
        p->pid = -1; // means we have particles in daughter nodes
    }
    
    // Deal with new particle now
    di = GetQuadrant<float>(p->center, it.x);
    if (p->Daughter[di] == NULL) {
        CreateNode(p, di);
    }
    AddParticle(p->Daughter[di], it);
    
    return 0;
}


/********** Distance **********/
/*! \fn float Distance(T x[], T y[], int &i_MaxD)
 *  \brief calculate distance between two locations, i_MaxD is index of MaxD */
template<typename T>
float Quadtree::Distance(T x[], T y[], int &i_MaxD)
{
    float d = 0, dr[2];
    for (int i = 0; i != 2; i++) {
        dr[i] = fabs(x[i] - y[i]);
        if (dr[i] > MaxD[i_MaxD][i]) {
            dr[i] = vf->L[i] - dr[i];
        }
        d += dr[i] * dr[i];
    }
    return sqrt(d);
}

/********** EvaluateAccurateSphere **********/
/*! \fn void EvaluateAccurateSphere(QuadtreeNode *p, T x[], long &npar, float &R, int &i_MaxD)
 *  \brief calculate npar for an accurate sphere centered at x with r=Radius */
template<typename T>
void Quadtree::EvaluateAccurateSphere(QuadtreeNode *p, T x[], long &npar, float &R, int &i_MaxD)
{
    float d = Distance<T>(p->center, x, i_MaxD);
    if (d > R + s2o2 * p->Lx) {
        // completely outside radius
        ;
    } else if (d < R - s2o2 * p->Lx) {
        // completely inside radius
        npar += p->np;
    } else {
        // intersect
        // if it is lead node
        if (p->np == 1 || p->level == level_limit) {
            if (Distance<float>(x, pl->List[p->pid].x, i_MaxD) <= R) {
                npar++;
            }
            for (int i = 0; i < p->deep_pids.size(); i++) {
                if (Distance<float>(x, pl->List[p->deep_pids[i]].x, i_MaxD) <= R) {
                    npar++;
                }
            }
        } else {
            for (int i = 0; i != 4; i++) {
                if (p->Daughter[i] != NULL) {
                    EvaluateAccurateSphere(p->Daughter[i], x, npar, R, i_MaxD);
                }
            }
        }
    }
}

/********** MaxSigmapPerLevel **********/
/*! \fn void SigmapMaxPerLevel(int file_i)
 *  \brief Find the max rhop within a sphere with radius of N*dx */
void Quadtree::SigmapMaxPerLevel(int file_i)
{
#ifdef ENABLE_MPI

    int Npoints = vf->dimensions[0]*vf->dimensions[1];
    int **indices, index = 0, goodboy;
    indices = new int*[Npoints];
    
    // list all the indices of points we need to check
    for (int npo = 0; npo != Npoints; npo++) {
        indices[npo] = new int[2];
        indices[npo][1] = npo / vf->dimensions[0];
        indices[npo][0] = npo % vf->dimensions[0];
    }
    
    if (myMPI->myrank == 0) {
        vector<int> working (myMPI->numprocs-1, 1); // used for ending the jobs
        for (int i = 1; i < min(Npoints, myMPI->numprocs); i++) {
            MPI::COMM_WORLD.Send(&index, 1, MPI::INT, i, i);
            //cout << "Master: sent point " << index << " to Processor " << i << endl;
            index++;
        }
        
        while (index < Npoints) {
            MPI::COMM_WORLD.Recv(&goodboy, 1, MPI::INT, MPI::ANY_SOURCE, MPI::ANY_TAG, myMPI->status);
            // in case we need something else
            //sender = myMPI->status.Get_source();
            //tag = myMPI->status.Get_tag();
            MPI::COMM_WORLD.Send(&index, 1, MPI::INT, goodboy, goodboy);
            if (index % 4000 == 0) {
                cout << "Master: sent point " << index << " to Processor " << goodboy << endl;
            }
            index++;
        }
        index = -1;
        while (working.size() > 0) {
            MPI::COMM_WORLD.Recv(&goodboy, 1, MPI::INT, MPI::ANY_SOURCE, MPI::ANY_TAG, myMPI->status);
            MPI::COMM_WORLD.Send(&index, 1, MPI::INT, goodboy, goodboy);
            working.pop_back();
        }
        cout << "Master: points distribution is done." << endl;
    } else {
        long npar;
        long mine_contribution = 0;
        // for increasing sample when it comes to smallest sphere
        float *temp_cellcenter;
        short used_for_and_op = 1, selection, xy, yesorno;
        
        //cout << "Worker " << myMPI->myrank << ": I'm in. " << endl;
        MPI::COMM_WORLD.Recv(&index, 1, MPI::INT, 0, myMPI->myrank, myMPI->status);
        while (index != -1) {
            // workers do the calculations
            temp_cellcenter = vf->cell_center[0][indices[index][1]][indices[index][0]];
            for (int i = 0; i <= level; i++) {
                if (Radius[i] >= shortestL) {
                    break; // further step make no sense, so save some computation time
                }
                npar = 0;
                
                // Count real particle numbers in accurate sphere
                EvaluateAccurateSphere<float>(root, temp_cellcenter, npar, Radius[i], i);
                
                if (npar > SMPL[file_i][i]) {
                    SMPL[file_i][i] = npar;
                }
                
                // some small circles need more sampling
                if (i < 4) {
                    for (selection = 1; selection != 4; selection++) {
                        npar = 0;
                        for (xy = 0; xy != 2; xy++) {
                            yesorno = (selection >> xy) & used_for_and_op;
                            tempcenter[xy] = temp_cellcenter[xy] - yesorno * Radius[0];
                        }
                        EvaluateAccurateSphere(root, tempcenter, npar, Radius[i], i);
                        if (npar > SMPL[file_i][i]) {
                            SMPL[file_i][i] = npar;
                        }
                        if (i == 0) {
                            float tempcenter2[3];
                            short selection2, xy2, yesorno2;
                            for (selection2 = 1; selection2 != 4; selection2++) {
                                npar = 0;
                                for (xy2 = 0; xy2 != 2; xy2++) {
                                    yesorno2 = (selection2 >> xy2) & used_for_and_op;
                                    tempcenter2[xy2] = tempcenter[xy2] - yesorno2 * Radius[0]/2.0;
                                }
                                EvaluateAccurateSphere(root, tempcenter2, npar, Radius[i], i);
                                if (npar > SMPL[file_i][i]) {
                                    SMPL[file_i][i] = npar;
                                }
                            }
                        }
                    }
                }
                /*
                 if (myMPI->myrank == 2 && mine_contribution == 6750) {
                 cout << "z=" << temp_cellcenter[2] << ", R=" << Radius[i] << ", npar=" << npar << endl;
                 }
                 */
            }
            mine_contribution++;
            MPI::COMM_WORLD.Send(&myMPI->myrank, 1, MPI::INT, 0, 0);
            MPI::COMM_WORLD.Recv(&index, 1, MPI::INT, 0, myMPI->myrank, myMPI->status);
        }
        cout << "Worker " << myMPI->myrank << ": " << mine_contribution << " jobs are done." << endl;
    }
    
    // memory control
    for (int npo = 0; npo != Npoints; npo++) {
        delete [] indices[npo];
    }
    delete [] indices;
    
    
#else /* ENABLE_MPI */
    cout << "Such a heavy computation will take forever in single CPU. Please consider using parallel computing." << endl;
    exit(1);
#endif /* ENABLE_MPI */
    
}






/********************************/
/********** OctreeNode **********/
/********************************/

/********** Constructor **********/
OctreeNode::OctreeNode()
{
    level = 0;
    np = 0;
    Father = NULL;
    Daughter.resize(8);
    for (vector<OctreeNode *>::iterator it = Daughter.begin(); it != Daughter.end(); it++) {
        (*it) = NULL;
    }
}

/********** Destructor **********/
OctreeNode::~OctreeNode()
{
    if (Daughter.size() > 0) {
        vector<OctreeNode *> temp;
        Daughter.swap(temp);
    }
}


/****************************/
/********** Octree **********/
/****************************/

/********** Constructor **********/
Octree::Octree(int n_file)
{
    s3o2 = sqrt(3.0)/2.0;
    foPio3 = 4.0*PI/3.0;
    MaxD = NULL;
    RMPL = new float*[n_file];
    for (int i = 0; i != n_file; i++) {
        RMPL[i] = NULL;
    }
    Radius = NULL;
    Initialize();
}

/********** Initialize **********/
int Octree::Initialize()
{
    root = new OctreeNode;
    Max_Rhop = 0;
    RpAV = 0;
    RpSQ = 0;
    RpQU = 0;
    shortestL = 1e10;
    return 0;
}

/********** Destructor **********/
Octree::~Octree()
{
    CleanMem(root);
    delete root; root = NULL;
    for (int i = 0; i != fio->n_file; i++) {
        if (RMPL[i] != NULL) {
            delete [] RMPL[i];
        }
    }
    if (RMPL != NULL) delete [] RMPL;
    if (Radius != NULL) delete [] Radius;
}

/********** CleanMem **********/
void Octree::CleanMem(OctreeNode *p)
{
    for (vector<OctreeNode *>::iterator it = p->Daughter.begin(); it != p->Daughter.end(); it++) {
        if ((*it) != NULL) {
            CleanMem((*it));
            delete (*it);
            (*it) = NULL;
        }
    }
}

/********** BuildTree **********/
/*! \fn BuildTree(VtkFile *VF, ParticleList *PL, int file_i)
 *  \brief Build the whole tree */
int Octree::BuildTree(VtkFile *VF, ParticleList *PL, int file_i)
{
    CleanMem(root);
    delete root; root = NULL;
    Initialize();
    vf = VF;
    pl = PL;
    
    /***** first assign root node and calculate m1par *****/
    int rootNx = 0;
    for (int i = 0; i != 3; i++) {
        if (rootNx < vf->dimensions[i]) {
            rootNx = vf->dimensions[i];
        }
        if (shortestL > vf->L[i]) {
            shortestL = vf->L[i];
        }
    }
    for (int i = 0; i != 3; i++) {
        root->center[i] = 0.0; // this code only considers special cases
    }
    root->Lx = rootNx * vf->spacing[0]; // assuming spacing are all the same
    m1par = 0.02*sqrt(2*PI)*vf->L[0]*vf->L[1]/pl->n;
    //if (vf->L[0] != vf->L[1] || vf->L[1] != vf->L[2]) {
    //    cout << "WARNING: this code should be used only for simualations in cubic boxes. " << endl;
    //    exit(1);
    //}
    
    /***** The following need the "level" to allocate space *****/
    level = round(log10(rootNx)/log10(2.0)); // count root level
    level_limit = round(log10(root->Lx/0.0001)/log10(2.0)); // stop at 1e-4
    
    if (MaxD == NULL) {
        MaxD = new float*[level+1];
        if (Radius == NULL) {
            Radius = new float[level+1];
        }
        for (int i = 0; i <= level; i++) {
            MaxD[i] = new float[3];
            Radius[i] = vf->spacing[0] * pow(2.0, i-1); // used in later func
            MaxD[i][0] = vf->L[0] - Radius[i];
            MaxD[i][1] = vf->L[1] - Radius[i];
            MaxD[i][2] = vf->L[2] - Radius[i];
#ifdef OutflowRate
            MaxD[i][2] = vf->L[2] * 2;
            // let it bigger than any distance so no periodic computation happen in x3 direction
            // in fact all particles are in the middle plane, and the RMPL[file_i][root_level] is calculated by another way. This is not main concern
#endif
            //cout << myMPI->prank() << "MaxD[" << i << "]=" << pvector<float>(MaxD[i]) << ", R[i]=" << Radius[i] << endl;
        }
    }
    if (RMPL[file_i] == NULL) {
        RMPL[file_i] = new float[level+1];
        for (int i = 0; i <= level; i++) {
            RMPL[file_i][i] = 0;
        }
    }

    /***** compute particle numbers for each node *****/
    for (vector<Particle>::iterator it = pl->List.begin(); it != pl->List.end(); ++it) {
        AddParticle(root, *it);
    } //*/
    
    return 0;
}

/********** GetOctant **********/
/*! \fn int GetOctant(float center[3], T pos[3]);
 *  \brief return the index of the daughter */
template<typename T>
int Octree::GetOctant(float center[3], T pos[3])
{
    int i = 0;
    i += (pos[0]>center[0]);     // x
    i += (pos[1]>center[1])<<1;  // y
    i += (pos[2]>center[2])<<2;  // z
    return i;
}

/********** CreateNode **********/
/*! \fn OctreeNode *CreateNode(OctreeNode *p, int &octant)
 *  \brief create new node depend on father node and octant */
OctreeNode * Octree::CreateNode(OctreeNode *p, int &octant)
{
    OctreeNode *q = new OctreeNode;
    q->level = p->level + 1;
    q->Father = p;
    q->Lx = p->Lx / 2;
    q->center[0] = p->center[0] + q->Lx * ((octant & 1) - 0.5);
    q->center[1] = p->center[1] + q->Lx * (bool(octant & 2) - 0.5);
    q->center[2] = p->center[2] + q->Lx * ((octant > 3) - 0.5);
    p->Daughter[octant] = q;
    return q;
}


/********** AddParticle **********/
/*! \fn int AddParticle(OctreeNode *p, Particle &it)
 *  \brief assign particle to one node */
int Octree::AddParticle(OctreeNode *p, Particle &it)
{
    p->np++; // anyway we need to do this
    
    // if it is an empty leaf previously, just put particle id int
    if (p->np == 1) {
        p->pid = it.pid;
        return 0;
    }
    // if it touch the level limit, just put particle id into vector
    if (p->level == level_limit) {
        p->deep_pids.push_back(it.pid);
        return 0;
    }
    
    int di;
    // Or else it has one already and level < level limit, move it to daughter first
    if (p->np == 2) {
        di = GetOctant<float>(p->center, pl->List[p->pid].x);
        AddParticle(CreateNode(p, di), pl->List[p->pid]);
        p->pid = -1; // means we have particles in daughter nodes
    }
    
    // Deal with new particle now
    di = GetOctant<float>(p->center, it.x);
    if (p->Daughter[di] == NULL) {
        CreateNode(p, di);
    }
    AddParticle(p->Daughter[di], it);
    
    return 0;
}

/********** Distance **********/
/*! \fn float Distance(T x[3], T y[3], int &i_MaxD)
 *  \brief calculate distance between two locations, i_MaxD is index of MaxD */
template<typename T>
float Octree::Distance(T x[3], T y[3], int &i_MaxD)
{
    float d = 0, dr[3];
    for (int i = 0; i != 3; i++) {
        dr[i] = fabs(x[i] - y[i]);
        if (dr[i] > MaxD[i_MaxD][i]) {
            dr[i] = vf->L[i] - dr[i];
        }
        d += dr[i] * dr[i];
    }
    return sqrt(d);
}

/*! \fn void EvaluateAccurateSphere(OctreeNode *p, T x[3], long &npar, float &R, int &i_MaxD)
 *  \brief calculate rhop for an accurate sphere centered at x with r=Radius */
template<typename T>
void Octree::EvaluateAccurateSphere(OctreeNode *p, T x[3], long &npar, float &R, int &i_MaxD)
{
    float d = Distance<T>(p->center, x, i_MaxD);
    if (d > R + s3o2 * p->Lx) {
        // completely outside radius
        ;
    } else if (d < R - s3o2 * p->Lx) {
        // completely inside radius
        npar += p->np;
    } else {
        // intersect
        // if it is leaf node
        if (p->np == 1 || p->level == level_limit) {
            if (Distance<float>(x, pl->List[p->pid].x, i_MaxD) <= R) {
                npar++;
            }
            for (int i = 0; i < p->deep_pids.size(); i++) {
                if (Distance<float>(x, pl->List[p->deep_pids[i]].x, i_MaxD) <= R) {
                    npar++;
                }
            }
        } else {
            for (int i = 0; i != 8; i++) {
                if (p->Daughter[i] != NULL) {
                    EvaluateAccurateSphere(p->Daughter[i], x, npar, R, i_MaxD);
                }
            }
        }
    }
}

/********** MaxRhopPerLevel **********/
/*! \fn void RhopMaxPerLevel(int file_i)
 *  \brief Find the max rhop within a sphere with radius of N*dx */
void Octree::RhopMaxPerLevel(int file_i)
{
    int ks = -1, ke = -1, iz, iy, ix;
    
    // in case needed
    //RMPL[file_i][0] = Max_Rhop;
    //RMPL[file_i][level] = RpAV * (vf->kpe-vf->kps) / vf->dimensions[2];
    
    // certainly those outside the center zone (with rhop = 0) can be omitted
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
    
#ifdef ENABLE_MPI
    if (myMPI->myrank == 0) {
        cout << "Pr" << myMPI->myrank << ": ks = " << ks << ", ke = " << ke << ", m1par = " << m1par << endl;
    }
    
    int Npoints = vf->dimensions[0]*vf->dimensions[1]*(ke-ks+1), temp = vf->dimensions[0]*vf->dimensions[1];
    int **indices, index = 0, goodboy;
    indices = new int*[Npoints];
    
    // list all the indices of points we need to check
    for (int npo = 0; npo != Npoints; npo++) {
        indices[npo] = new int[3];
        indices[npo][2] = npo / temp + ks;
        indices[npo][1] = (npo % temp) / vf->dimensions[0];
        indices[npo][0] = (npo % temp) % vf->dimensions[0];
        /*
        if (myMPI->myrank == 0) {
            cout << indices[npo][2] << " " << indices[npo][1] << " " << indices[npo][0] << endl;
        }
         */
    }

    if (myMPI->myrank == 0) {
        vector<int> working (myMPI->numprocs-1, 1); // used for ending the jobs
        for (int i = 1; i < min(Npoints, myMPI->numprocs); i++) {
            MPI::COMM_WORLD.Send(&index, 1, MPI::INT, i, i);
            //cout << "Master: sent point " << index << " to Processor " << i << endl;
            index++;
        }
        
        while (index < Npoints) {
            MPI::COMM_WORLD.Recv(&goodboy, 1, MPI::INT, MPI::ANY_SOURCE, MPI::ANY_TAG, myMPI->status);
            // in case we need something else
            //sender = myMPI->status.Get_source();
            //tag = myMPI->status.Get_tag();
            MPI::COMM_WORLD.Send(&index, 1, MPI::INT, goodboy, goodboy);
            if (index % 50000 == 0) {
                cout << "Master: sent point " << index << " to Processor " << goodboy << endl;
            }
            index++;
        }
        index = -1;
        while (working.size() > 0) {
            MPI::COMM_WORLD.Recv(&goodboy, 1, MPI::INT, MPI::ANY_SOURCE, MPI::ANY_TAG, myMPI->status);
            MPI::COMM_WORLD.Send(&index, 1, MPI::INT, goodboy, goodboy);
            working.pop_back();
        }
        cout << "Master: points distribution is done." << endl;
    } else {
        long npar;
        long mine_contribution = 0;
        // for increasing sample when it comes to smallest sphere
        float *temp_cellcenter;
        short used_for_and_op = 1, selection, xyz, yesorno;
        
        //cout << "Worker " << myMPI->myrank << ": I'm in. " << endl;
        MPI::COMM_WORLD.Recv(&index, 1, MPI::INT, 0, myMPI->myrank, myMPI->status);
        while (index != -1) {
            // workers do the calculations
            temp_cellcenter = vf->cell_center[indices[index][2]][indices[index][1]][indices[index][0]];
            for (int i = 0; i <= level; i++) {
                if (Radius[i] >= shortestL) {
                    break; // further step make no sense, so save some computation time
                }
                npar = 0;
                
                // Count real particle numbers in accurate sphere
                EvaluateAccurateSphere<float>(root, temp_cellcenter, npar, Radius[i], i);
                
                if (npar > RMPL[file_i][i]) {
                    RMPL[file_i][i] = npar;
                }
                
                // some small sphere need more sample
                if (i < 4) {
                    for (selection = 1; selection != 8; selection++) {
                        npar = 0;
                        for (xyz = 0; xyz != 3; xyz++) {
                            yesorno = (selection >> xyz) & used_for_and_op;
                            tempcenter[xyz] = temp_cellcenter[xyz] - yesorno * Radius[0];
                        }
                        EvaluateAccurateSphere(root, tempcenter, npar, Radius[i], i);
                        if (npar > RMPL[file_i][i]) {
                            RMPL[file_i][i] = npar;
                        }
                        if (i == 0) {
                            float tempcenter2[3];
                            short selection2, xyz2, yesorno2;
                            for (selection2 = 1; selection2 != 8; selection2++) {
                                npar = 0;
                                for (xyz2 = 0; xyz2 != 3; xyz2++) {
                                    yesorno2 = (selection2 >> xyz2) & used_for_and_op;
                                    tempcenter2[xyz2] = tempcenter[xyz2] - yesorno2 * Radius[0]/2.0;
                                }
                                EvaluateAccurateSphere(root, tempcenter2, npar, Radius[i], i);
                                if (npar > RMPL[file_i][i]) {
                                    RMPL[file_i][i] = npar;
                                }
                            }
                        }
                    }
                }
                /*
                if (myMPI->myrank == 2 && mine_contribution == 6750) {
                    cout << "z=" << temp_cellcenter[2] << ", R=" << Radius[i] << ", npar=" << npar << endl;
                }
                 */
            }
            mine_contribution++;
            MPI::COMM_WORLD.Send(&myMPI->myrank, 1, MPI::INT, 0, 0);
            MPI::COMM_WORLD.Recv(&index, 1, MPI::INT, 0, myMPI->myrank, myMPI->status);
        }
        cout << "Worker " << myMPI->myrank << ": " << mine_contribution << " jobs are done." << endl;
    }
    
    // memory control
    for (int npo = 0; npo != Npoints; npo++) {
        delete [] indices[npo];
    }
    delete [] indices;

    
#else /* ENABLE_MPI */
    cout << "Such a heavy computation will take forever in single CPU. Please consider using parallel computing." << endl;
    exit(1);
#endif /* ENABLE_MPI */

    
    
}






