#ifndef NEWTONMIN_H_
#define NEWTONMIN_H_

#include "matrix.h"

#define FDVAL 1e-15
#define EPSILON 2.7182818284590452353602874713526624977572
#define SETITER 1000

#ifdef __cplusplus
extern "C" {
#endif

double signx(double x);

double func4(double *x,int N);

double l2norm(double *vec,int N);

double cholmod(double *A, int N, double *L, double maxinp);

double modelhess(double *A,int N,double *dx,double *L);

void linsolve_lower(double *L,int N,double *b,double *x);

void grad_fd(double (*funcpt)(double *,int),double *x,int N,double *dx,double *f);

void hessian_fd(double (*funcpt)(double *,int),double *x,int N,double *dx,double *f);

void hessian_fd2(double (*funcpt)(double *,int),double *x,int N,double *dx,double *f);

int lnsrch(double (*funcpt)(double *,int),double *xi,double *jac,double *p,int N,double * dx,double maxstep,double stol,double *x);

int lnsrchmod(double (*funcpt)(double *,int),double *xi,double *jac,double *p,int N,double * dx,double maxstep,
	double stol,double *x,double *jacf);
	
int lnsrchcg(double (*funcpt)(double *,int),double *xi,double *jac,double *p,int N,double * dx,double maxstep,
	double stol,double *x,double *jacf);
	
int stopcheck(double fx,int N,double *xc,double *xf,double *jac,double *dx,double fsval,double gtol,double stol,int retval);

int newton_min_func(double (*funcpt)(double *,int),double *xi,int N,double *dx,double fsval,double *xf);

int trsrch(double (*funcpt)(double *,int),double *xi,double *jac,double *sN,int N,double * dx,double maxstep,
		int iter,double *L,double *hess,double stol,double *ioval,double *x);
		
void trstep(double *jac,double *sN,int N,double * dx,double *L,double *hess,double nlen,double *ioval,double *step);

int trupdate(double (*funcpt)(double *,int),double *xi,double *jac,double *step,int N,double * dx,double maxstep,
		int retcode,double *L,double *hess,double stol,int method,double *ioval,double *xprev,double *x);		

void trstep_ddl(double *jac,double *sN,int N,double * dx,double maxstep,double *L,double *hess,double nlen,double *ioval,
		double *ssd,double *v,double *step);		
		
int trsrch_ddl(double (*funcpt)(double *,int),double *xi,double *jac,double *sN,int N,double * dx,double maxstep,
		int iter,double *L,double *hess,double stol,double *ioval,double *x);		

int newton_min_trust(double (*funcpt)(double *,int),double *xi,int N,double *dx,double fsval,double delta,int method,double *xf);		

#ifdef __cplusplus
}
#endif

#endif /* NEWTONMIN_H_ */