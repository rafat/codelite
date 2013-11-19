#ifndef CONJGRAD_H_
#define CONJGRAD_H_

#include "newtonmin.h"


#ifdef __cplusplus
extern "C" {
#endif

int ichol(double *A, int N);

int stopcheck2(double fx,int N,double *xc,double *xf,double *jac,double *dx,double fsval,double gtol,double stol) ;

int cgr(double *xi,int N,double *A,double *b,double *xf);

int cgpc(double *xi,int N,double *A,double *b,double *xf);

int cgp(double (*funcpt)(double *,int),double *xi,int N,double *A,double *b,double *xf);

int conjgrad_min_lin(double (*funcpt)(double *,int),double *xi,int N,double *xf);

#ifdef __cplusplus
}
#endif

#endif /* CONJGRAD_H_ */