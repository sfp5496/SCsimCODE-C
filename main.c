/*
 * Author: Sean Peterson
 * Advisor: Dr. Scott Franklin
 * Program: SCsim.c
 * What it does: randomly generates npart spherolines within a cylindrical container then proceeds
 *  to grow them from volume zero until they fill their container
 *
 * Inputs:
 *  argv[1]: name of the directory where you want to store packing files
 *  argv[2]: parameter file
 *  argv[3]: file which contains an initial packing
 *
 * Parameter Files:
 *  the files should include the following values in order
 *  npart R H ETA DPHI PHI ALPHA CUBE LENGTH WIDTH HEIGHT
 *  
 *  npart: number of particles
 *  R: radius of container, which is a cylinder
 *  H: height of container, which is a cylinder
 *  ETA: used to be a scalar used in the conjugate gradient method but isn't used anymore
 *  DPHI: the increase in our packing fraction each iteration
 *  PHI: the packing fraction that the simulation ends at
 *  ALPHA: the aspect ratio of our particles
 *  CUBE: if set to 1, our container is a cube instead of a cylinder (WIP)
 *  LENGTH: length of the cube (WIP)
 *  WIDTH: width of the cube (WIP)
 *  HEIGHT: height of the cube (WIP)
 *
 * Packing Files:
 *  iter pid r[i].x r[i].y r[i].z u[i].x u[i].y u[i].z
 *  iter: the current iteration number of the program
 *  pid: particle ID, just tells you where its data is in the arrays
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "double2.h"
#include "double3.h"
#include "d2func.h"
#include "d3func.h"
#include "Enfunc.h"
#include "Gradfunc.h"
#include "Misc.h"
#include "Globalvars.h"

FILE* fp;
FILE* inputs;

void start() {
	/*
	 * Arguments: none
	 * Returns: allocates the necessary memory for our arrays and uses initcond() for to generate
	 *  the initial packing and converts it to cartesian coordinates
	 */
	r=malloc(npart*sizeof(double3));
	r1=malloc(npart*sizeof(double3));
	r2=malloc(npart*sizeof(double3));
	u=malloc(npart*sizeof(double3));
	theta=malloc(npart*sizeof(double2));
	l=malloc(npart*sizeof(double));
	sigma=malloc(npart*sizeof(double));
	GU1=malloc(npart*sizeof(double3));
	GU1A=malloc(npart*sizeof(double2));
	rc=malloc(npart*sizeof(double3));
	initcond(rc,theta,l,sigma);
	h=malloc(npart*sizeof(double3));
	hA=malloc(npart*sizeof(double2));
	GU0=malloc(npart*sizeof(double3));
	GU0A=malloc(npart*sizeof(double2));
	int i;
	for(i=0;i<npart;i++) {
		r[i].x=rc[i].x*cos(rc[i].y);
		r[i].y=rc[i].x*sin(rc[i].y);
		r[i].z=rc[i].z;
		sptoca(i);
	}
	phi=packfrac();
}

void doAnIter() {
	/*
	 * Arguments: none
	 * Returns: Runs through the conjugate gradient method for every particle in every dimension 
	 *  until the sum of the gradients is small enough to say the system is in a potential energy 
	 *  minima then increases phi by DPHI and updates l and sigma
	 */
	int i;
	Gradient();
	/*for(i=0;i<npart;i++) {
		ConGrad1(i);
	}*/
	int iter=0;
	int s=0;
	double GS=GradSum();
//	Gradient();
	for(i=0;i<npart;i++) {
		h[i]=d3null();
		hA[i]=d2null();
		GU0[i]=d3null();
		GU0A[i]=d2null();
	}

	//for(i=0;i<1;i++) {
	while((GradSum()>.0000001) && (GradSumA()>.000001)) {
		iter++;
		//for(i=0;i<npart;i++) {
			//ConGrad(i); //perform the conjugate gradient method for the cartesian coords
			//ConGradA(i); //perform the con grad method for theta
			//ConGradA1(i); //perform the con grad method for phi
			//ConGrad2(i);
		ConGrad2();

		//}
		Gradient();
		//This is to ensure we don't get stuck at a certain phi value
		if(iter>20000) {
			fprintf(stderr,"Too many iterations...\n");
			break;
		}
		int k;
		double en=0.0;
		for(k=0;k<npart;k++) {
			en+=PotEnergy(k);
		}
		//This breaks us out of our loop if the sum of the gradients doesn't change by much
		//over many iterations
		if(fabs(GS-GradSum())<0.01) {
			s++;
		}
		else {
			s=0;
		}
		if(s>5000) {
			break;
		}
		GS=GradSum();
	}
	double c=contacts();
	double p=packfrac();
	double U=collider();
	double P=Pressure();
	fprintf(stderr,"Rod Radius: %lf\n",sigma[0]/2.0);
	fprintf(stderr,"Rod Length: %lf\n",l[0]);
	fprintf(stderr,"Packing Fraction: %lf\n",p);
	fprintf(stderr,"Average Number of Contacts: %lf\n",c/npart);
	fprintf(stderr,"Potential Energy: %E\n",U);
	fprintf(stderr,"Average Overlap: %lf\n",sqrt((U/npart)*2.0)/(sigma[0]/2.0));
	fprintf(stderr,"Pressure on Walls: %lf\n",P);
	fprintf(stdout,"%E	%E	%E	%E	%E\n",p,c/npart,U,sqrt((U/npart)*2.0)/(sigma[0]/2.0),P);
	phi+=DPHI; //update phi
	updatephi(); //use the new phi value to get a new sigma and l
}

void DEMend() {
	/*
	 * Arguments: none
	 * Returns: doesn't really return anything, just frees all the memory that we used during the
	 *  simulation
	 */
	free(r);
	free(u);
	free(theta);
	free(rc);
	free(GU1);
	free(GU1A);
	free(r1);
	free(r2);
	free(sigma);
	free(l);
}

int main(int argc, char* argv[]) {
	if(argc>2) { //check for a parameter file
		const char* infile;
		inputs=fopen(argv[2],"r");
		fscanf(inputs,"%d %lf %lf %lf %lf %lf %lf %d %lf %lf %lf",&npart,&R,&H,&ETA,&DPHI,&PHI,&ALPHA,&CUBE,&LENGTH,&WIDTH,&HEIGHT);
	}
	else {
		fprintf(stderr,"No Parameter File given.\n");
		return 1;
	}
	if(argc>1) { //check for a directory to put configuration files
		struct stat st={0};

		if(stat(argv[1],&st)==-1) {
			fprintf(stderr,"Creating Directory %s...\n",argv[1]);
			mkdir(argv[1],0777);
		}
		else {
			fprintf(stderr,"Directory %s Already Exists.\n",argv[1]);
		}
	}
	else {
		fprintf(stderr,"No directory given, if you want to print here use './'.\n");
		return 1;
	}
	char* filename;
	filename=(char*)malloc(20*sizeof(char));
	int j=0;
	start(); //allocate memory and generate an initial packing
	if(argc>3) { //Check for an initial packing file
		inputs=fopen(argv[3],"r");
		int k;
		int iter,pid; //just places to send the variables I don't need
		fprintf(stderr,"Loading initial packing.\n");
		for(k=0;k<npart;k++) {
			fscanf(inputs,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf",&iter,&pid,sigma+k,l+k,&(r[k].x),&(r[k].y),&(r[k].z),&(u[k].x),&(u[k].y),&(u[k].z));
			theta[k].x=atan2(u[k].y,u[k].x);
			theta[k].y=atan(u[k].z);
		}
		phi=packfrac();
		//double X=collider();
		//double X2=collider2();
		//fprintf(stderr,"POTENTIAL ENERGY: %lf %lf\n",X,X2);
	}
	else {
		fprintf(stderr,"No initial packing given...\n");
		fprintf(stderr,"Using randomly generated packing.\n");
	}
	phi=0.0;
	//Gradient();
	//int k;
	//for(k=0;k<1;k++) {
	while(phi<PHI) {
		doAnIter();
		int q;
		//fprintf(stderr,"PRINT CHECK: %i\n",(int)(phi*1000)%10);
		if(((((int)(phi/DPHI))%(int)(.01/DPHI)))==0) {
			sprintf(filename,"%s/%03d.dat",argv[1],(int)(phi*100));
			fp=fopen((const char*)filename,"w");
			for(q=0;q<npart;q++) {
				fprintf(fp,"%i %i %lf %lf %lf %lf %lf %lf %lf %lf\n",j,q,sigma[q],l[q],r[q].x,r[q].y,r[q].z,u[q].x,u[q].y,u[q].z);
			}
			j++;
			fclose(fp);
		}
	}
	DEMend();
	return 0;
}
