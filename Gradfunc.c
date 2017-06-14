#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "double2.h"
#include "double3.h"
#include "d2func.h"
#include "d3func.h"
#include "Enfunc.h"
#include "Globalvars.h"

void Grad(int i) {
	/*
	 * Arguments: particle id
	 * Returns: the gradient of the particle's potential energy in each cartesian direction, this
	 *  is stored in an array of gradient vectos GU1
	 */

	GU0[i]=GU1[i];
	GU0A[i]=GU1A[i];

	GU1[i].x=0.0;
	GU1[i].y=0.0;
	GU1[i].z=0.0;
	GU1A[i]=d2null();
	double U_0,U_1x,U_1y,U_1z,U_1t,U_1p=0.0;
	
	//Partial Derivative wrt x_i
	double dx=.0001;
	U_0=PotEnergy2(i);
	r[i].x+=dx;
	U_1x=PotEnergy2(i);
	r[i].x-=dx;
	GU1[i].x=(U_1x-U_0)/dx;

	//Partial Derivative wrt y_i
	double dy=.0001;
	r[i].y+=dy;
	U_1y=PotEnergy2(i);
	r[i].y-=dy;
	GU1[i].y=(U_1y-U_0)/dy;

	//Partial Derivative wrt z_i
	double dz=.0001;
	r[i].z+=dz;
	U_1z=PotEnergy2(i);
	r[i].z-=dz;
	GU1[i].z=(U_1z-U_0)/dz;

	//Partial Derivative wrt theta_i
	double dt=.0001;
	theta[i].x+=dt;
	U_1t=PotEnergy2(i);
	theta[i].x-=dt;
	GU1A[i].x=(U_1t-U_0)/dt;

	//Partial Derivative wrt phi_i
	double dp=.0001;
	theta[i].y+=dp;
	U_1p=PotEnergy2(i);
	theta[i].y-=dp;
	GU1A[i].y=(U_1p-U_0)/dp;

	//printf("theta[%i]: %lf %lf\n",i,theta[i].x,theta[i].y);
	//printf("GU1A[%i]: %E %E\n",i,GU1A[i].x,GU1A[i].y);
	//printf("U2A[%i]: %E %E\n",i,U_1t,U_1p);
	//printf("U1[%i]: %E\n",i,U_0);
	//printf("GU1[%i]: %E %E %E\n",i,GU1[i].x,GU1[i].y,GU1[i].z);
}

void Gradient() {
	int i;
	for(i=0;i<npart;i++) {
		Grad(i);
	}
}

double GradSum() {
	/*
	 * Arguments: none
	 * Returns: the sum of the magnitudes of the gradient vector of the potential energy of every
	 *  particle in the cartesian directions, mostly used in calculating whether or not a given 
	 *  packing may be minimized further or not
	 */
	int i;
	double ret=0;
	for(i=0;i<npart;i++) {
		ret+=d3mag(GU1[i]);
	}
	return ret;
}

double GradSumA() {
	/*
	 * Arguments: none
	 * Returns: the sum of the magnitudes of the gradient vector of the potential energy of every
	 *  particle in the theta and phi directions, mostly used in calculating whether or not a given 
	 *  packing may be minimized further or not
	 */
	int i;
	double ret=0;
	for(i=0;i<npart;i++) {
		ret+=hypot(GU1A[i].x,GU1A[i].y);
	}
	return ret;
}

void ConGrad2() {
	double eta=0.0;
	double etaA=0.0;
	double gamma;
	double gammaA;

	eta=-fabs(GradSum());
	if(fabs(eta)>1.0) {
		eta=-1.0;
	}
	else if(fabs(eta)<.001) {
		eta=-.001;
	}
	etaA=-fabs(GradSumA());
	if(fabs(etaA)>1.0) {
		eta=-1.0;
	}
	else if(etaA<.001) {
		etaA=-.001;
	}

	int i;
	for(i=0;i<npart;i++) {
		if(d3dotp(GU0[i],GU0[i])!=0.0) {
			gamma=d3dotp(GU1[i],GU1[i])/d3dotp(GU0[i],GU0[i]);
			if(gamma>1.0) {
				gamma=1.0;
			}
			else if(gamma<-1.0) {
				gamma=-1.0;
			}
		}
		else {
			gamma=1.0;
		}
		
		h[i].x=GU1[i].x+gamma*h[i].x;
		h[i].y=GU1[i].y+gamma*h[i].y;
		h[i].z=GU1[i].z+gamma*h[i].z;
		double3 dr;
		dr=d3multscal(h[i],eta);
		if(fabs(d3mag(dr))>0.5*R) {
			dr=d3multscal(dr,0.5*R/d3mag(dr));
		}
		r[i]=d3add(r[i],dr);

		etaA=-fabs(GradSumA());
		if(fabs(etaA)>1.0) {
			etaA=-1.0;
		}
		else if(fabs(etaA)<.001) {
			etaA=-.001;
		}
		
		if(d2mag(GU0A[i])!=0) {
			gammaA=d2dotp(GU1A[i],GU1A[i])/d2dotp(GU0A[i],GU0A[i]);
			if(gammaA>1.0) {
				gammaA=1.0;
			}
			else if(gammaA<-1.0) {
				gammaA=-1.0;
			}
		}
		else {
			gammaA=1.0;
		}
	
		hA[i]=d2add(GU1A[i],d2multscal(hA[i],gammaA));
		theta[i]=d2add(theta[i],d2multscal(hA[i],etaA));
	
//		printf("r[%i]: %lf %lf %lf\n",i,r[i].x,r[i].y,r[i].z);
//		printf("theta[%i]: %lf %lf\n",i,theta[i].x,theta[i].y);
		//printf("h[%i]: %lf %lf %lf\n",i,h[i].x,h[i].y,h[i].z);
		//printf("eta[%i]: %lf\n",i,eta);
		//printf("GradSumA(): %lf\n",GradSumA());
		//printf("gamma[%i]: %lf\n",i,gamma);
//		printf("GU1A[%i]: %lf %lf\n",i,GU1A[i].x,GU1A[i].y);
		//printf("GU0[%i]: %lf %lf %lf\n",i,GU0[i].x,GU0[i].y,GU0[i].z);
	}
}

