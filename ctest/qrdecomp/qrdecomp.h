#ifndef QRDECOMP_H_
#define QRDECOMP_H_

#include "matrix.h"

#define TOL 1e-12

#ifdef __cplusplus
extern "C" {
#endif

void eye(double *mat,int N);

double house(double*x,int N,double *v);

void housemat(double *v, int N,double beta,double *mat);

void qr_house(double *A, int M, int N,double *bvec);

void getQR_house(double *A,int M,int N,double *bvec,double *Q, double *R);

void hessenberg(double *A,int N);

void francisQR(double *A,int N);

void eig22(double *A, int stride,double *eigre,double *eigim);

int francis_iter(double *A, int N, double *H);

void eig(double *A,int N,double *eigre,double *eigim);

#ifdef __cplusplus
}
#endif

#endif /* QRDECOMP_H_ */