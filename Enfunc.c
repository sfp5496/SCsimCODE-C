#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "double2.h"
#include "double3.h"
#include "d2func.h"
#include "d3func.h"
#include "Globalvars.h"

double d3SCdist(double3 ri, double3 rj, double3 ui, double3 uj, double li, double lj) {
	/*
	 * Arguments: the position and orientation vectors of two spherolines and their corresponding
	 *  lambda values
	 * Returns: the distance between the two spherolines
	 */
	double3 dij;
	dij=d3sub(d3sub(d3add(rj,d3multscal(uj,lj)),ri),d3multscal(ui,li));
	return d3mag(dij);
	
}

double lambda(double3 ri, double3 rj, double3 ui, double3 uj, double li) {
	/*
	 * Arguments: the position and orientation of two spherolines and the length of the 
	 *  spheroline whose lambda we are attempting to find
	 * Returns: the lambda value of the ith spheroline of the ij spheroline pair
	 * What is This: the lambda values are the points on the line of the spherolines, 
	 *  where the two most closely intersect
	 */
	double retn,retd;
	retn=d3dotp(ui,d3sub(rj,ri))-d3dotp(ui,uj)*d3dotp(uj,d3sub(rj,ri));
	retd=1.0-(d3dotp(ui,uj)*d3dotp(ui,uj));
	if((retn/retd)>(li/2.0)) {
		return (li/2.0);
	}
	else if((retn/retd)<(-li/2.0)) {
		return (-li/2.0);
	}
	else {
		return (retn/retd);
	}
}

void sptoca(int i) {
	/*
	 * Arguments: the particle id
	 * Returns: updates the orientation vector of the ith particle based on the theta and phi values
	 */
	u[i].x=sin(theta[i].y)*cos(theta[i].x);
	u[i].y=sin(theta[i].y)*sin(theta[i].x);
	u[i].z=cos(theta[i].y);
}

void ends(int i) {
	/*
	 * Arguments: the particle id
	 * Returns: calculates where the ends of the particle are, useful in calculating potential 
	 *  energy between a particle and a wall
	 */
	r1[i]=d3add(r[i],d3multscal(u[i],0.5*l[i]));
	r2[i]=d3sub(r[i],d3multscal(u[i],0.5*l[i]));
	r1[i].x=hypot(r1[i].x,r1[i].y);
	r2[i].x=hypot(r2[i].x,r2[i].y);
}

double PotEnergydev(int i, int j) {
	/*
	 * Arguments:
	 *    ri: position of ith particle
	 *    rj: position of jth particle
	 *    ui: orientation of ith particle
	 *    uj: orientation of jth particle
	 *    li: length of ith particle
	 *    lj: length of jth particle
	 *    si: diameter of ith particle
	 *    sj: diameter of jth particle
	 *    R: container Radius
	 *    H: container Height
	 * 
	 * Returns: the potential energy between the two particles
	 */
	double U=0.0;
	//Check to see if the particles even have a chance to be touching
	if((d3dist(r[i],r[j])<(l[i]+sigma[i]+l[j]+sigma[j])/2.0)) {
		double lambda_i, lambda_j;
		lambda_i=lambda(r[i],r[j],u[i],u[j],l[i]);
		lambda_j=lambda(r[j],r[i],u[j],u[i],l[j]);
		double d;
		d=d3SCdist(r[i],r[j],u[i],u[j],lambda_i,lambda_j);
		//Check to see if the particles are ACTUALLY in contact
		if(d<(sigma[i]+sigma[j])/2.0) {
			U=0.25*((sigma[i]+sigma[j])/2.0-d)*((sigma[i]+sigma[j])/2.0-d);
		}
	}
	return U;
}

double WallEnergy(int i) {
	/*
	 * Arguments:
	 *    r1: position of one end of the particle
	 *    r2: position of the other end of the particle
	 *    li: length of the particle
	 *    si: diameter of the particle
	 *    R: container radius
	 *    H: container height
	 *
	 * Returns: the potential energy between the particle and the walls
	 */
	
	double U=0.0;

	//check to see if the first end is in contact with the radial wall of the container
	if(r1[i].x>R-(sigma[i]/2.0)) {
		U+=0.5*(r1[i].x-(R-sigma[i]/2.0))*(r1[i].x-(R-sigma[i]/2.0));
	}
	//check to see if the second end is in contact with the radial wall of the container
	if(r2[i].x>R-(sigma[i]/2.0)) {
		U+=0.5*(r2[i].x-(R-sigma[i]/2.0))*(r2[i].x-(R-sigma[i]/2.0));
	}
	//check to see if the first end is in contact with the floor
	if(r1[i].z-sigma[i]/2.0<0.0) {
		U+=0.5*(r1[i].z-sigma[i]/2.0)*(r1[i].z-sigma[i]/2.0);
	}
	//check to see if the first end is in contact with the ceiling
	else if(r1[i].z+sigma[i]/2.0>H) {
		U+=0.5*(r1[i].z+sigma[i]/2.0-H)*(r1[i].z+sigma[i]/2.0-H);
	}
	//check to see if the second end is in contact with the floor
	if(r2[i].z-sigma[i]/2.0<0.0) {
		U+=0.5*(r2[i].z-sigma[i]/2.0)*(r2[i].z-sigma[i]/2.0);
	}
	//check to see if the second end is in contact with the ceiling
	else if(r2[i].z+sigma[i]/2.0>H) {
		U+=0.5*(r2[i].z+sigma[i]/2.0-H)*(r2[i].z+sigma[i]/2.0-H);
	}
	return U;
}

double PotEnergy(int i) {
	/*
	 * Arguments: particle id
	 * Returns: the potential energy of the ith particle with relation to the wall and every other
	 *  particle
	 */
	int k;
	for(k=0;k<npart;k++) {
		sptoca(k);
	}
	double ret=0.0;
	ends(i);
	int j;
	for(j=0;j<npart;j++) {
		//Check to see if the particles are even close to each other
		if((i!=j) && (d3dist(r[i],r[j])<(l[i]+sigma[i]+l[j]+sigma[j])/2.0)) {
			double lambda_i, lambda_j;
			lambda_i=lambda(r[i],r[j],u[i],u[j],l[i]);
			lambda_j=lambda(r[j],r[i],u[j],u[i],l[j]);
			double d;
			d=d3SCdist(r[i],r[j],u[i],u[j],lambda_i,lambda_j);
			//check to see if the particles are in contact
			if(d<(sigma[i]+sigma[j])/2.0) {
				ret+=0.25*((sigma[i]+sigma[j])/2.0-d)*((sigma[i]+sigma[j])/2.0-d);
			}
		}
	}
	
	//check to see if the first end is in contact with the radial wall of the container
	if(r1[i].x>R-(sigma[i]/2.0)) {
		ret+=0.5*(r1[i].x-(R-sigma[i]/2.0))*(r1[i].x-(R-sigma[i]/2.0));
	}
	//check to see if the second end is in contact with the radial wall of the container
	if(r2[i].x>R-(sigma[i]/2.0)) {
		ret+=0.5*(r2[i].x-(R-sigma[i]/2.0))*(r2[i].x-(R-sigma[i]/2.0));
	}
	//check to see if the first end is in contact with the floor
	if(r1[i].z-sigma[i]/2.0<0.0) {
		ret+=0.5*(r1[i].z-sigma[i]/2.0)*(r1[i].z-sigma[i]/2.0);
	}
	//check to see if the first end is in contact with the ceiling
	else if(r1[i].z+sigma[i]/2.0>H) {
		ret+=0.5*(r1[i].z+sigma[i]/2.0-H)*(r1[i].z+sigma[i]/2.0-H);
	}
	//check to see if the second end is in contact with the floor
	if(r2[i].z-sigma[i]/2.0<0.0) {
		ret+=0.5*(r2[i].z-sigma[i]/2.0)*(r2[i].z-sigma[i]/2.0);
	}
	//check to see if the second end is in contact with the ceiling
	else if(r2[i].z+sigma[i]/2.0>H) {
		ret+=0.5*(r2[i].z+sigma[i]/2.0-H)*(r2[i].z+sigma[i]/2.0-H);
	}
	
	return ret;
}

double PotEnergy2(int i) {
	double ret=0.0;
	int j;
	for(j=0;j<npart;j++) {
		sptoca(j);
	}
	ends(i);
	for(j=0;j<npart;j++) {
		if(i==j) {
			ret+=WallEnergy(i);
		}
		else {
			ret+=PotEnergydev(i,j);
		}
	}
	return ret;
}

double collider() {
	/*
	 * Arguments: none
	 * Returns: the total potential energy of the current packing
	 */
	int i;

	double ret=0.0;
	for(i=0;i<npart;i++) {
		ret+=PotEnergy(i);
	}
	return ret;
}

double collider2() {
	int i,j;

	for(i=0;i<npart;i++) {
		sptoca(i);
	}

	double ret=0.0;
	for(i=0;i<npart;i++) {
		for(j=0;j<npart;j++) {
			if(i==j) {
				ends(i);
				ret+=WallEnergy(i);
			}
			else {
				ret+=PotEnergydev(i,j);
			}
		}
	}
	return ret;
}


