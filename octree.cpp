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
    s3o2 = sqrt(3.0)/2.0;
    foPio3 = 4.0*PI/3.0;
    NxCubic = NULL;
    MaxD = NULL;
    RMPL = NULL;
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
    if (NxCubic != NULL) delete [] NxCubic;
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
/*! \fn BuildTree(VtkFile *VF, ParticleList *PL)
 *  \brief Build the whole tree */
int Octree::BuildTree(VtkFile *VF, ParticleList *PL)
{
    CleanMem(root);
    delete root; root = NULL;
    Initialize();
    vf = VF;
    pl = PL;
    
    /***** first assign root node and calculate m1par *****/
    root->Nx = 0;
    for (int i = 0; i != 3; i++) {
        if (root->Nx < vf->dimensions[i]) {
            root->Nx = vf->dimensions[i];
        }
        if (shortestL > vf->L[i]) {
            shortestL = vf->L[i];
        }
    }
    for (int i = 0; i != 3; i++) {
        root->center[i] = 0.0; // this code only considers special cases
    }
    root->Lx = root->Nx * vf->spacing[0]; // assuming spacing are all the same
    m1par = 0.02*sqrt(2*PI)*vf->L[0]*vf->L[1]/pl->n;
    //if (vf->L[0] != vf->L[1] || vf->L[1] != vf->L[2]) {
    //    cout << "WARNING: this code should be used only for simualations in cubic boxes. " << endl;
    //    exit(1);
    //}
    
    /***** The following need the "level" to allocate space *****/
    level = int(log10(root->Nx)/log10(2.0)); // count root level
    if (NxCubic == NULL) { // only used in EvaluateOctreeSphere()
        NxCubic = new long[level+1];
        for (int i = 0; i <= level; i++) {
            NxCubic[i] = long(pow(root->Nx/pow(2.0, i), 3));
        }
    }
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
            // in fact all particles are in the middle plane, and the RMPL[root_level] is calculated by another way. This is not main concern
#endif
            //cout << myMPI->prank() << "MaxD[" << i << "]=" << pvector<float>(MaxD[i]) << ", R[i]=" << Radius[i] << endl;
        }
    }
    if (RMPL == NULL) {
        RMPL = new float[level+1];
        for (int i = 0; i <= level; i++) {
            RMPL[i] = 0;
        }
    }
    
    /***** Add cells to activate nodes *****/
    OctreeNode *p;
    double temprhop, tempSQ, tempV;
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
    
    /***** compute particle numbers for each node *****/
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
            // check if the point is on the domain surface
            int i = 0;
            i += (it.x[0]>=p->center[0]);     // x
            i += (it.x[1]>=p->center[1])<<1;  // y
            i += (it.x[2]>=p->center[2])<<2;  // z
            if (p->Daughter[i] == NULL) {
                cout << "Should not create node while adding particles.\n";
                cout << "par.x=" << pvector(it.x) << " center=" << pvector(p->center) << " di=" << di << " level=" << p->level << endl;
                exit(1);
            } else {
                di = i;
            }
        }
        p->np++;
        p = p->Daughter[di];
    }
    p->np++;
    p->parlist.push_back(it.pid); // only leaf node need it
    
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

/********** EvaluateOneP **********/
/*! \fn void EvaluateOctreeSphere(OctreeNode *p, T x[3], double &rp, long &cells,  float &R)
 *  \brief calculate rhop for a sphere centered at x with r=Radius */
template<typename T>
void Octree::EvaluateOctreeSphere(OctreeNode *p, T x[3], double &rp, long &cells, float &R)
{
    float d = Distance<T>(p, x);
    if (d > R + s3o2 * p->Lx) {
        // completely outside radius
        ;
    } else if (d < R - s3o2 * p->Lx) {
        // completely inside radius
        rp += p->rhop;
        cells += NxCubic[p->level];
    } else {
        // mostly intersect with radius
        if (p->level == level) {
            if (d <= R) {
                rp += p->rhop;
                cells++;
            }
        } else {
            for (int i = 0; i != 8; i++) {
                if (p->Daughter[i] != NULL) {
                    EvaluateOctreeSphere(p->Daughter[i], x, rp, cells, R);
                }
            }
        }
    }
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
        if (p->level == level) {
            for (long it = 0; it < p->parlist.size(); it++) {
                if (Distance<float>(x, pl->List[p->parlist[it]].x, i_MaxD) <= R) {
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
/*! \fn void RhopMaxPerLevel()
 *  \brief Find the max rhop within a sphere with radius of N*dx */
void Octree::RhopMaxPerLevel()
{
    int ks = -1, ke = -1, iz, iy, ix;
    
    // in case needed
    //RMPL[0] = Max_Rhop;
    //RMPL[level] = RpAV * (vf->kpe-vf->kps) / vf->dimensions[2];
    
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
        double rp; long cells; long npar;
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
                rp = 0; cells = 0; npar = 0;

                /* Count octree-style approximate sphere
                EvaluateOctreeSphere<float>(root, vf->cell_center[indices[index][2]][indices[index][1]][indices[index][0]], rp, cells, Radius[i]);
                rp /= cells;
                //*/
                
                // Count real particle numbers in accurate sphere
                EvaluateAccurateSphere<float>(root, temp_cellcenter, npar, Radius[i], i);
                
                if (npar > RMPL[i]) {
                    RMPL[i] = npar;
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
                        if (npar > RMPL[i]) {
                            RMPL[i] = npar;
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
                                if (npar > RMPL[i]) {
                                    RMPL[i] = npar;
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






