#ifndef SECANT_H_
#define SECANT_H_

#include "conjgrad.h"


#ifdef __cplusplus
extern "C" {
#endif

void bfgs_naive(double *H,int N,double *xi,double *xf,double *jac,double *jacf);

int bfgs_min_naive(double (*funcpt)(double *,int),double *xi,int N,double *dx,double fsval,double *xf);

void bfgs_factored(double *H,int N,double *xi,double *xf,double *jac,double *jacf);

#ifdef __cplusplus
}
#endif

#endif /* SECANT_H_ */