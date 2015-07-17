//
//  readlis.cpp
//  readvtklis
//
//  Created by Rixin Li on 1/14/15.
//  Copyright (c) 2015 Rixin Li. All rights reserved.
//

#include "readlis.h"

/********** Constructor **********/
ParticleList::ParticleList()
{
    ;
}

/********** Read particle list from file **********/
/*! \fn int ReadLis(string filename)
 *  \brief Read particle list from file */
int ParticleList::ReadLis(string filename)
{
    ifstream file (filename.c_str(), ios::binary);
    if (file.is_open()) {
        file.read((char *)(coorlim), 12*sizeof(float));
        // you can read an array in a bundle, like above
        // The values in coorlim are the coordinate limits of this grid and domain:
        // x1l, x1u, x2l, x2u, x3l, x3u, x1dl, x1du, x2dl, x2du, x3dl, x3du
        // here l means lower limit, u means upper limit, d means domain
        
        file.read((char *)(&ntype), sizeof(float));
        typeinfo = (float *)malloc(ntype*sizeof(float));
        for (int i = 0; i < ntype; i++) {
            file.read((char *)(&typeinfo[i]), sizeof(float));
        }
        file.read((char *)(&time), sizeof(float));
        file.read((char *)(&dt), sizeof(float));
        file.read((char *)(&n), sizeof(long));
        if (fio->ParNum_flag == 1 && fio->HeiPar_flag == 0) {
            return 0;
        }
#ifdef RESERVE_PUSH_BACK
        /*** reserve and push_back ***/
        List.reserve(n);
        for (long i = 0; i < n; i++) {
            Particle temp;
            for (int j = 0; j < 3; j++) {
                file.read((char *)(&temp.x[j]), sizeof(float));
            }
            for (int j = 0; j < 3; j++) {
                file.read((char *)(&temp.v[j]), sizeof(float));
            }
            
            file.read((char *)(&temp.rad), sizeof(float));
            file.read((char *)(&temp.mass), sizeof(float));
            file.read((char *)(&temp.pid), sizeof(long));
            file.read((char *)(&temp.cpuid), sizeof(int));
            List.push_back(temp);
        }
        /*** reserve and push_back ***/
#endif
        
#ifdef FROM_ARRAY_TO_VECTOR
        /*** construct array and give it to vector ***/
        Particle *temp = new Particle;
        temp = (Particle *)malloc(n*sizeof(Particle));
        for (long i = 0; i < n; i++) {
            for (int j = 0; j < 3; j++) {
                file.read((char *)(&temp[i].x[j]), sizeof(float));
            }
            for (int j = 0; j < 3; j++) {
                file.read((char *)(&temp[i].v[j]), sizeof(float));
            }
            file.read((char *)(&temp[i].rad), sizeof(float));
            file.read((char *)(&temp[i].mass), sizeof(float));
            file.read((char *)(&temp[i].pid), sizeof(long));
            file.read((char *)(&temp[i].cpuid), sizeof(int));
        }
        vector<Particle> tempList(temp, temp+n);
        List.swap(tempList);
        delete temp;
        /*** construct array and give it to vector ***/
#endif
        
#ifdef RESIZE_LIST
        /*** resize list ***/
        if (List.size() != n) {
            List.resize(n);
        }
        for (vector<Particle>::iterator it = List.begin(); it != List.end(); it++) {
            for (int j = 0; j < 3; j++) {
                file.read((char *)(&it->x[j]), sizeof(float));
            }
            for (int j = 0; j < 3; j++) {
                file.read((char *)(&it->v[j]), sizeof(float));
            }
            file.read((char *)(&it->rad), sizeof(float));
            file.read((char *)(&it->mass), sizeof(float));
            file.read((char *)(&it->pid), sizeof(long));
            file.read((char *)(&it->cpuid), sizeof(int));
        }
        /*** resize list ***/
#endif
        
        
        
        
    } else {
        cout << "Failed to open " << filename << "\n";
    }
    file.close();
    return 0;
}

/********** Calculate the scale height of partiles **********/
/*! \fn float ScaleHeight(double &Hp, double &Hp_in1sigma);
 *  \brief calculate the scale height of partiles */
int ParticleList::ScaleHeight(double &Hp, double &Hp_in1sigma)
{
    double one_sigma = 0.682689492137;
    if (List.size() == 0) {
        Hp = 0;
        Hp_in1sigma = 0;
        return 0;
    }
    vector<float> par_z;
    par_z.resize(List.size());
    Hp = 0;
    for (long i = 0; i != List.size(); i++) {
        Hp += List[i].x[2]*List[i].x[2];
        par_z[i] = List[i].x[2];
    }
    Hp = sqrt(Hp/n);
    
    sort(par_z.begin(), par_z.end());
    vector<float> par_z_reverse(par_z);
    sort(par_z_reverse.begin(), par_z_reverse.end(), std::greater<double>());
    while ((par_z.size()+par_z_reverse.size()-List.size())/((double)List.size()) > one_sigma) {
        if (par_z.back() > fabs(par_z_reverse.back())) {
            par_z.pop_back();
        } else {
            par_z_reverse.pop_back();
        }
    }
    Hp_in1sigma = par_z.back()?: par_z.back() > fabs(par_z_reverse.back());
    
    return 0;
}

/********** Free List to control memory **********/
/*! \fn int InitializeList()
 *  \brief free List to control memory */
int ParticleList::InitializeList()
{
    vector<Particle> temp;
    List.swap(temp);
    return 0;
}

/********** Print basic info **********/
/*! \fn int PrintInfo()
 *  \brief print basic info */
int ParticleList::PrintInfo()
{
    cout << "N_par = " << n << "\n";
    cout << "N_par_type = " << ntype << "\n";
    cout << "Mass of particles:" << "\n";
    for (int i = 0; i != ntype; i++) {
        cout << typeinfo[i] << "\n";
    }
    
    return 0;
}

/********** Get numprocs **********/
/*! \fn int GetNumprocs()
 *  \brief get the number of processors */
int ParticleList::GetNumprocs()
{
    int n_cpu = 0;
    // this is assuming each cpu contains particles, wrong idea
    for (vector<Particle>::iterator it = List.begin(); it != List.end(); ++it) {
        if (it->cpuid > n_cpu) {
            n_cpu = it->cpuid;
        }
    }
    return n_cpu+1;
}

/********** Destructor **********/
ParticleList::~ParticleList()
{
    if (List.size() > 0) {
        vector<Particle> temp;
        List.swap(temp);
    }
}

