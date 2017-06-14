#ifndef __ENFUNC_H__
#define __ENFUNC_H__

double d3SCdist(double3 ri,double3 rj,double3 ui,double3 uj,double li,double lj);
double lambda(double3 ri,double3 rj,double3 ui,double3 uj,double li);
void sptoca(int i);
void ends(int i);
double PotEnergydev(int i,int j);
double WallEnergy(int i);
double PotEnergy(int i);
double PotEnergy2(int i);
double collider();
double collider2();


#endif
