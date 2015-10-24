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
    return 0;
}

/********** Destructor **********/
Octree::~Octree()
{
    delete root;
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
    //root->FreeNode();
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
    level = int(log10(root->Nx)/log10(2.0)); // count root level
    maxrhop = (double *)calloc(level+1, sizeof(double));
    
    for (int i = 0; i != 3; i++) {
        root->center[i] = 0.0; // this code only considers special cases
    }
    root->Lx = root->Nx * vf->spacing[0]; // assuming spacing are all the same
    
    
    // Add cells to activate nodes
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
                        tempSQ = p->rhop * p->rhop;
                        RpSQ += tempSQ;
                        RpQU += tempSQ * tempSQ;
                    }
                }
            }
        }
    }
    tempV = (vf->kpe-vf->kps)*vf->dimensions[1]*vf->dimensions[0];
    RpAV = root->rhop/tempV;
    RpSQ = sqrt(RpSQ/tempV);
    RpQU = sqrt(sqrt(RpQU/tempV));
    
    /* compute particle numbers for each node
    for (vector<Particle>::iterator it = pl->List.begin(); it != pl->List.end(); ++it) {
        AddParticle(it->x);
    } //*/
    
    
    free(maxrhop);
    return 0;
}

/********** AddCell **********/
/*! \fn OctreeNode *AddCell(double cc[3], double rhop)
 *  \brief using cell center to find the tree node, create nodes if needed */
OctreeNode *Octree::AddCell(double cc[3], double rhop)
{
    OctreeNode *p = root, *q;
    for (int l = 0; l != level; l++) {
        int di = GetOctant<double>(p->center, cc); // daughter index
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
/*! \fn int GetOctant(double pos[3]);
 *  \brief return the index of the daughter */
template<typename T>
int Octree::GetOctant(double center[3], T pos[3])
{
    int i = 0;
    i += (pos[0]>center[0]);     // x
    i += (pos[1]>center[1])<<1;  // y
    i += (pos[2]>center[2])<<2;  // z
    return i;
}

/********** AddParticle **********/
/*! \fn int AddParticle(float x[3])
 *  \brief assign particle to one node */
int Octree::AddParticle(float x[3])
{
    OctreeNode *p = root;
    for (int l = 0; l != level; l++) {
        int di = GetOctant<float>(p->center, x); // daughter index
        if (p->Daughter[di] == NULL) {
            cout << "Should not create node while adding particles.\n";
            exit(1);
        }
        p->np++;
        p = p->Daughter[di];
    }
    p->np++;
    
    return 0;
}

/********** Distance **********/
/*! \fn double Distance(OctreeNode *p, float x[3])
 *  \brief calculate distance between cell center and particle */
template<typename T>
double Octree::Distance(OctreeNode *p, T x[3])
{
    double d = 0, dr[3];
    for (int i = 0; i != 3; i++) {
        dr[i] = fabs(x[0] - p->center[0]);
        if (dr[i] > MaxD[0]) {
            dr[i] = vf->L[0] - dr[i];
        }
        d += dr[i] * dr[i];
    }
    return d;
}


/********** EvaluateOneP **********/
/*! \fn void EvaluateOnePoint(OctreeNode *p, T x[3], double *rp)
 *  \brief calculate rhop for a sphere centered at x with r=Radius */
template<typename T>
void Octree::EvaluateOnePoint(OctreeNode *p, T x[3], double *rp)
{
    double d;
    d = Distance<T>(p, x);
    if (d > Radius + s3o2 * p->Lx) {
        ;
    } else if (p->level != level && d < Radius - s3o2 * p->Lx) {
        //cout << "hit\n";
        *rp += p->rhop;
    } else if (p->level == level && d < Radius) {
        //cout << "hit2\n";
        *rp += p->rhop;
    } else {
        for (int i = 0; i != 8; i++) {
            if (p->Daughter[i] != NULL) {
                EvaluateOnePoint(p->Daughter[i], x, rp);
            }
        }
    }
}







