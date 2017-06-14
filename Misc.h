#ifndef __MISC_H__
#define __MISC_H__

double unitrand();
void initcond(double3* x,double2* ang,double* len,double* diam);
double packfrac();
void updatephi();
double variable_dphi(double dphi,double U);
double Pressure();
double contacts();

#endif
