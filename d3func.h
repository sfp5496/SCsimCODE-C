#ifndef __D3FUNC_H__
#define __D3FUNC_H__

double3 d3add(double3 a, double3 b);
double3 d3sub(double3 a, double3 b);
double3 d3multscal(double3 a, double b);
double3 d3multscal2(double3 a, double b);
double3 d3divscal(double3 a, double b);
double d3dotp(double3 a, double3 b);
double3 d3crossp(double3 a, double3 b);
double d3mag(double3 a);
double d3dist(double3 a, double3 b);
double3 d3unit(double3 a);
double3 d3null();

#endif
