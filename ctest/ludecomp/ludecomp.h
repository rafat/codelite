#ifndef LUDECOMP_H_
#define LUDECOMP_H_

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

int pludecomp(double *A,int N,int *ipiv);

int linsolve(double *LU,int N,double *b,double *x,int *ipiv);

#ifdef __cplusplus
}
#endif

#endif /* MATRIX_H_ */