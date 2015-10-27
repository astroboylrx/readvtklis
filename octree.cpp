//
//  octree.cpp
//  readvtklis
//
//  Created by Rixin Li on 10/12/15.
//  Copyright © 2015 Rixin Li. All rights reserved.
//

#include "octree.h"

/********************************/
/********** OctreeNode **********/
/********************************/

/********** Constructor **********/
OctreeNode::OctreeNode()
{
    level = 0;
    rhop = 0;
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
Octree::Octree()
{
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
    s3o2 = sqrt(3.0)/2.0;
    foPio3 = 4.0*PI/3.0;
    return 0;
}

/********** Destructor **********/
Octree::~Octree()
{
    CleanMem(root);
    delete root; root = NULL;
    delete [] NxCubic;
    delete [] MaxD;
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
/*! \fn BuildTree(VtkFile *VF, ParticleList *PL)
 *  \brief Build the whole tree */
int Octree::BuildTree(VtkFile *VF, ParticleList *PL)
{
    CleanMem(root);
    delete root; root = NULL;
    Initialize();
    vf = VF;
    pl = PL;
    
    root->Nx = 0;
    for (int i = 0; i != 3; i++) {
        if (root->Nx < vf->dimensions[i]) {
            root->Nx = vf->dimensions[i];
        }
    }
    for (int i = 0; i != 3; i++) {
        root->center[i] = 0.0; // this code only considers special cases
    }
    root->Lx = root->Nx * vf->spacing[0]; // assuming spacing are all the same
    
    level = int(log10(root->Nx)/log10(2.0)); // count root level
    NxCubic = new long[level+1];
    MaxD = new float*[level+1];
    m1par = 0.02*sqrt(2*PI)*root->Lx*root->Lx/pl->n;
    
    if (vf->L[0] != vf->L[1] || vf->L[1] != vf->L[2]) {
        if (myMPI->myrank == 0) {
            cout << "WARNING: this code should be used only for simualations in cubic boxes. " << endl;
        }
    }

    for (int i = 0; i <= level; i++) {
        NxCubic[i] = long(pow(root->Nx/pow(2.0, i), 3));
        MaxD[i] = new float[3];
        
        Radius = root->Lx/pow(2.0, i+1);
        MaxD[i][0] = vf->L[0] - Radius;
        MaxD[i][1] = vf->L[1] - Radius;
        MaxD[i][2] = vf->L[2] - Radius;
        
#ifdef OutflowRate
        MaxD[i][2] = vf->L[2] * 2;
        // let it bigger than any distance so no periodic computation happen in x3 direction
        // in fact all particles are in the middle plane, and the RMPL[root_level] is calculated by another way. This is not main concern
#endif
        if (myMPI->myrank == 0) {
            cout << "NxCubic[" << i << "] = " << NxCubic[i] << ", MaxD[" << i << "] = " << MaxD[i][0] << endl;
        }
    }
    
    // Add cells to activate nodes
    OctreeNode *p;
    float temprhop, tempSQ, tempV;
    for (int iz = 0; iz != vf->dimensions[2]; iz++) {
        for (int iy = 0; iy != vf->dimensions[1]; iy++) {
            for (int ix = 0; ix != vf->dimensions[0]; ix++) {
                temprhop = vf->cd_scalar[1].data[iz][iy][ix];
                if (temprhop != 0) {
                    p = AddCell(vf->cell_center[iz][iy][ix], vf->cd_scalar[1].data[iz][iy][ix]);
                    if (Max_Rhop < p->rhop) {
                        Max_Rhop = p->rhop;
                    }
                    if (iz >= vf->kps and iz < vf->kpe) {
                        RpAV += p->rhop;
                        tempSQ = p->rhop * p->rhop;
                        RpSQ += tempSQ;
                        RpQU += tempSQ * tempSQ;
                    }
                }
            }
        }
    }
    tempV = (vf->kpe-vf->kps)*vf->dimensions[1]*vf->dimensions[0];
    RpAV = RpAV/tempV;
    RpSQ = sqrt(RpSQ/tempV);
    RpQU = sqrt(sqrt(RpQU/tempV));
    
    // compute particle numbers for each node
    for (vector<Particle>::iterator it = pl->List.begin(); it != pl->List.end(); ++it) {
        AddParticle(*it);
    } //*/
    
    
    return 0;
}

/********** AddCell **********/
/*! \fn OctreeNode *AddCell(float cc[], float &rhop)
 *  \brief using cell center to find the tree node, create nodes if needed */
OctreeNode *Octree::AddCell(float cc[], float &rhop)
{
    OctreeNode *p = root, *q;
    for (int l = 0; l != level; l++) {
        int di = GetOctant<float>(p->center, cc); // daughter index
        if (p->Daughter[di] == NULL) {
            q = new OctreeNode;
            q->level = p->level + 1;
            q->Nx = p->Nx / 2;
            q->Father = p;
            q->Lx = p->Lx / 2;
            for (int i = 0; i != 3; i++) {
                q->center[i] = p->center[i] - q->Lx/2.0 + (cc[i]>p->center[i])*q->Lx;
            }
            p->Daughter[di] = q;
        }
        p->rhop += rhop;
        p = p->Daughter[di];
    }
    p->rhop += rhop;
    return p;
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

/********** AddParticle **********/
/*! \fn int AddParticle(Particle &it)
 *  \brief assign particle to one node */
int Octree::AddParticle(Particle &it)
{
    OctreeNode *p = root;
    for (int l = 0; l != level; l++) {
        int di = GetOctant<float>(p->center, it.x); // daughter index
        if (p->Daughter[di] == NULL) {
            cout << "Should not create node while adding particles.\n";
            exit(1);
        }
        p->np++;
        p = p->Daughter[di];
    }
    p->np++;
    p->parlist.push_back(it.pid); // only leaf node need it
    
    return 0;
}

/********** Distance **********/
/*! \fn float Distance(OctreeNode *p, float x[3])
 *  \brief calculate distance between cell center and particle */
template<typename T>
float Octree::Distance(OctreeNode *p, T x[3])
{
    float d = 0, dr[3];
    for (int i = 0; i != 3; i++) {
        dr[i] = fabs(x[i] - p->center[i]);
        if (dr[i] > MaxD[p->level][i]) {
            dr[i] = vf->L[i] - dr[i];
        }
        d += dr[i] * dr[i];
    }
    return sqrt(d);
}

/********** EvaluateOneP **********/
/*! \fn void EvaluateOctreeSphere(OctreeNode *p, T x[3], float *rp, long *cells)
 *  \brief calculate rhop for a sphere centered at x with r=Radius */
template<typename T>
void Octree::EvaluateOctreeSphere(OctreeNode *p, T x[3], float *rp, long *cells)
{
    float d;
    d = Distance<T>(p, x);
    if (d > Radius + s3o2 * p->Lx) {
        // completely outside radius
        ;
    } else if (d < Radius - s3o2 * p->Lx) {
        // completely inside radius
        *rp += p->rhop;
        *cells += NxCubic[p->level];
    } else {
        // mostly intersect with radius
        if (p->level == level) {
            if (d <= Radius) {
                *rp += p->rhop;
                (*cells)++;
            }
        } else {
            for (int i = 0; i != 8; i++) {
                if (p->Daughter[i] != NULL) {
                    EvaluateOctreeSphere(p->Daughter[i], x, rp, cells);
                }
            }
        }
    }
}

/*! \fn void EvaluateAccurateSphere(OctreeNode *p, T x[3], long &npar)
 *  \brief calculate rhop for an accurate sphere centered at x with r=Radius */
template<typename T>
void Octree::EvaluateAccurateSphere(OctreeNode *p, T x[3], long &npar)
{
    float d;
    d = Distance<T>(p, x);
    if (d > Radius + s3o2 * p->Lx) {
        // completely outside radius
        ;
    } else if (d < Radius - s3o2 * p->Lx) {
        // completely inside radius
        npar += p->np;
    } else {
        // intersect
        if (p->level == level) {
            for (vector<long>::iterator it = p->parlist.begin(); it != p->parlist.end(); it++) {
                if (Distance<T>(p, pl->List[*it].x) <= Radius) {
                    npar++;
                }
            }
        } else {
            for (int i = 0; i != 8; i++) {
                if (p->Daughter[i] != NULL) {
                    EvaluateAccurateSphere(p->Daughter[i], x, npar);
                }
            }
        }
    }
}

/********** MaxRhopPerLevel **********/
/*! \fn void RhopMaxPerLevel()
 *  \brief Find the max rhop within a sphere with radius of N*dx */
void Octree::RhopMaxPerLevel()
{
    int ks = -1, ke = -1, iz, iy, ix;
    RMPL = new float[level+1];
    RMPL[0] = Max_Rhop;
    RMPL[level] = RpAV * (vf->kpe-vf->kps) / vf->dimensions[2];
    for (int i = 1; i < level; i++) {
        RMPL[i] = 0;
    }
    
    // certainly those outside the center zone (with rhop = 0) can be omitted
    for (iz = vf->dimensions[2]/2-1; iz != 0; iz--) {
        ks = -1;
        for (iy = 0; iy != vf->dimensions[1]; iy++) {
            for (ix = 0; ix != vf->dimensions[0]; ix++) {
                if (vf->cd_scalar[1].data[iz][iy][iz] != 0) {
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
    for (iz = vf->dimensions[2]/2; iz != vf->dimensions[2]; iz++) {
        ke = -1;
        for (iy = 0; iy != vf->dimensions[1]; iy++) {
            for (ix = 0; ix != vf->dimensions[0]; ix++) {
                if (vf->cd_scalar[1].data[iz][iy][iz] != 0) {
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
    
#ifdef ENABLE_MPI
    if (myMPI->myrank == 0) {
        cout << "Pr" << myMPI->myrank << ": ks = " << ks << ", ke = " << ke << endl;
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
            cout << "Master: sent point " << index << " to Processor " << i << endl;
            index++;
        }
        
        while (index < Npoints) {
            MPI::COMM_WORLD.Recv(&goodboy, 1, MPI::INT, MPI::ANY_SOURCE, MPI::ANY_TAG, myMPI->status);
            // in case we need something else
            //sender = myMPI->status.Get_source();
            //tag = myMPI->status.Get_tag();
            MPI::COMM_WORLD.Send(&index, 1, MPI::INT, goodboy, goodboy);
            if (index % 2500 == 0) {
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
        float rp; long cells; long npar;
        long mine_contribution = 0;
        cout << "Worker " << myMPI->myrank << ": I'm in. " << endl;
        MPI::COMM_WORLD.Recv(&index, 1, MPI::INT, 0, myMPI->myrank, myMPI->status);
        while (index != -1) {
            // workers do the calculations
            for (int i = 0; i <= level; i++) {
                rp = 0; cells = 0; npar = 0;
                Radius = vf->spacing[0] * pow(2.0, i-1); // calculating radius (= diameter/2.), assuming each cell are cubic

                /* Count octree-style approximate sphere
                EvaluateOctreeSphere<float>(root, vf->cell_center[indices[index][2]][indices[index][1]][indices[index][0]], &rp, &cells);
                rp /= cells;
                //*/
                
                // Count real particle numbers in accurate sphere
                EvaluateAccurateSphere<float>(root, vf->cell_center[indices[index][2]][indices[index][1]][indices[index][0]], npar);
                rp = m1par * npar / (foPio3 * Radius * Radius * Radius);
                
                if (rp > RMPL[i]) {
                    RMPL[i] = rp;
                }
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






