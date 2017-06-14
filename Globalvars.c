#include <stdlib.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>

#include "double2.h"
#include "double3.h"

double3* r; //array of displacement vectors
double3* rc; //array of displacement vectors in cylindrical coords, used in initcond() only
double3* r1; //array of displacement vectors of one end of particles
double3* r2; //array of displacement vectors of the other end of particles
double3* u; //array of orientation vectors
double2* theta; //array of theta and phi orientation values
double* l; //array of particle lengths
double* sigma; //array of particle diameters
double3* GU0; //array of cartesian gradient vectors (of U) from previous iteration
double3* GU1; //array of cartesian gradient vectors (of U)
double2* GU0A; //array of theta 
double2* GU1A; //array of theta and phi gradients (of U)
double phi; //packing fraction

double3* h;
double2* hA;

//THESE VALUES ARE IMPORTED FROM A FILE AS THE PARAMETERS OF THE SIMULATION//

int npart;
double R;
double H;
double ETA;
double DPHI;
double PHI;
double ALPHA;
int CUBE;
double LENGTH;
double WIDTH;
double HEIGHT;

