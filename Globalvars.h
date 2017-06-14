#ifndef __GLOBALVARS_H__
#define __GLOBALVARS_H__

extern double3* r; //array of displacement vectors
extern double3* rc; //array of displacement vectors in cylindrical coords, used in initcond() only
extern double3* r1; //array of displacement vectors of one end of particles
extern double3* r2; //array of displacement vectors of the other end of particles
extern double3* u; //array of orientation vectors
extern double2* theta; //array of theta and phi orientation values
extern double* l; //array of particle lengths
extern double* sigma; //array of particle diameters
extern double3* GU0; //array of cartesian gradient vectors (of U) from previous iteration
extern double3* GU1; //array of cartesian gradient vectors (of U)
extern double2* GU0A; //array of theta 
extern double2* GU1A; //array of theta and phi gradients (of U)
extern double phi; //packing fraction

extern double3* h;
extern double2* hA;

//THESE VALUES ARE IMPORTED FROM A FILE AS THE PARAMETERS OF THE SIMULATION//

extern int npart;
extern double R;
extern double H;
extern double ETA;
extern double DPHI;
extern double PHI;
extern double ALPHA;
extern int CUBE;
extern double LENGTH;
extern double WIDTH;
extern double HEIGHT;

#endif
