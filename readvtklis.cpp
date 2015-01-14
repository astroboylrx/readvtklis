#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <cstdlib>
#include <cmath>

using namespace std;

/********** Particle class **********/
class Particle{
private:

public:
	int cpuid;
	long pid; // particle ID in list
	float x[3], v[3], rad, mass; // partile infomation
};

/********** ParticleList class **********/
class ParticleList{
private:

public:
	// the number of types
	int ntype;
	// the number of particles
	long n;
	// the coordinate limit, info of partile type
	// the time, and time step
	float coorlim[12], *typeinfo, time, dt;
	// particle list pointer
	Particle *List;
};

/********** Cell data **********/
class VtkFile{
private:

public:
	int nx1, nx2, nx3, i, j, k;
	
};

class CellData_Scaler{
private:

public:
	int i;
};

class CellData_Vector{
private:

public:
	int i;
};

	
/********** Read particle list from file **********/
ParticleList *ReadLis(string filename)
{
	ParticleList *pl;
	pl = (ParticleList *)malloc(sizeof(ParticleList));
	ifstream file (filename.c_str(), ios::binary);
	if (file.is_open()) {
		for (int i = 0; i < 12; i++) {
			file.read((char *)(&pl->coorlim[i]), sizeof(float));
			//cout << pl->coorlim[i] << endl;
			// The values in coorlim are the coordinate limits of this grid and domain:
			// x1l, x1u, x2l, x2u, x3l, x3u, x1dl, x1du, x2dl, x2du, x3dl, x3du
			// here l means lower limit, u means upper limit, d means domain
		}
		file.read((char *)(&pl->ntype), sizeof(float));
		pl->typeinfo = (float *)malloc(pl->ntype*sizeof(float));
		for (int i = 0; i < pl->ntype; i++) {
			file.read((char *)(&pl->typeinfo[i]), sizeof(float));
		}
		file.read((char *)(&pl->time), sizeof(float));
		file.read((char *)(&pl->dt), sizeof(float));
		file.read((char *)(&pl->n), sizeof(long));
		pl->List = (Particle *)malloc(pl->n*sizeof(Particle));
		for (long i = 0; i < pl->n; i++) {
			for (int j = 0; j < 3; j++) {
				file.read((char *)(&pl->List[i].x[j]), sizeof(float));
			}
			for (int j = 0; j < 3; j++) {
				file.read((char *)(&pl->List[i].v[j]), sizeof(float));
			}
			file.read((char *)(&pl->List[i].rad), sizeof(float));
			file.read((char *)(&pl->List[i].mass), sizeof(float));
			file.read((char *)(&pl->List[i].pid), sizeof(long));
			file.read((char *)(&pl->List[i].cpuid), sizeof(int));
		}
	} else {
		cout << "Failed to open " << filename << endl;
	}
	file.close();
	return pl;
}





/********** main **********/
int main ()
{
	ParticleList *pl;
	string filename = "comb/Cout.0004.all.lis";
	pl = ReadLis(filename);
	float Hp = 0;
	for (long i = 0; i < pl->n; i++) {
		Hp += pl->List[i].x[2]*pl->List[i].x[2];
	}
	Hp = sqrt(Hp/pl->n);
	cout << "The scaled height of particles is " << Hp << endl;
	

	return 0;
}
