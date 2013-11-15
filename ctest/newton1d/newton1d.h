#ifndef NEWTON1D_H_
#define NEWTON1D_H_

#include "matrix.h"

#define FDVAL 1e-15
#define EPSILON 2.7182818284590452353602874713526624977572
#define SETITER 100

#ifdef __cplusplus
extern "C" {
#endif

//double (*funcpt)(double *,int) = NULL;

double func1(double *x,int N);

double func2(double *x,int N);

double func3(double *x,int N, double *jac,double *hess);

double func4(double *x,int N);

double solvefunc(double (*funcpt)(double *,int),double *x,int N);

double cholmod(double *A, int N, double *L, double maxinp);

double modelhess2(double *A,int N,double *L);

double modelhess(double *A,int N,double *dx,double *L);

void linsolve_lower(double *L,int N,double *b,double *x);

void newton_func(double (*funcpt)(double *,int,double *,double *),double *xi,int N,double *x);

void jacobian_fd(double (*funcpt)(double *,int),double *x,int N,double *f);

void hessian_fd(double (*funcpt)(double *,int),double *x,int N,double *f);

void newton_fd(double (*funcpt)(double *,int),double *xi,int N,double *x);

int lnsrch(double (*funcpt)(double *,int),double *xi,double *jac,double *p,int N,double maxstep,double * dx,double stol,double *x);

int stopcheck(double fx,int N,double *xc,double *xf,double *jac,double *dx,double fsval,double gtol,double stol,int retval);

int newton_min_func(double (*funcpt)(double *,int),double *xi,int N,double *dx,double fsval,double *xf);

#ifdef __cplusplus
}
#endif

#endif /* NEWTON1D_H_ */