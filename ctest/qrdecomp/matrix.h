/*
 * matrix.h
 *
 *  Created on: Jul 1, 2013
 *      Author: USER
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>

#define CUTOFF 192

#ifdef __cplusplus
extern "C" {
#endif

double array_max_abs(double *array,int N);

//void mmult(double* A, double *B, double *C,int ra,int ca, int rb, int cb);

void dtranspose(double *sig, int rows, int cols,double *col);

void stranspose(double *sig, int rows, int cols,double *col);

void rtranspose(double *m, int rows, int cols,double *n, int r, int c);

void ctranspose(double *sig, int rows, int cols,double *col);

void mtranspose(double *sig, int rows, int cols,double *col);

//int minverse(double *xxt, int p);

void mdisplay(double *A, int row, int col);

void madd(double* A, double* B, double* C,int rows,int cols);

void msub(double* A, double* B, double* C,int rows,int cols);

void scale(double *A, int rows, int cols, double alpha);

void nmult(double* A, double* B, double* C,int m,int n, int p);

void tmult(double* A, double* B, double* C,int m,int n, int p);

void recmult(double* A, double* B, double* C,int m,int n, int p,int sA,int sB, int sC);

void rmult(double* A, double* B, double* C,int m,int n, int p);

int findrec(int *a, int *b, int *c);

void add_zero_pad(double *X, int rows, int cols, int zrow, int zcol,double *Y);

void remove_zero_pad(double *X, int rows, int cols, int zrow, int zcol,double *Y);

void madd_stride(double* A, double* B, double* C,int rows,int cols,int sA,int sB,int sC);

void msub_stride(double* A, double* B, double* C,int rows,int cols,int sA,int sB,int sC);

void rmadd_stride(double* A, double* B, double* C,int rows,int cols,int p,int sA,int sB,int sC);

void rmsub_stride(double* A, double* B, double* C,int rows,int cols,int p,int sA,int sB,int sC);

void srecmult(double* A, double* B, double* C,int m,int n, int p,int sA,int sB,int sC);

void smult(double* A, double* B, double* C,int m,int n, int p);

void mmult(double* A, double* B, double* C,int m,int n, int p);

void ludecomp(double *A,int N,int *ipiv);

void linsolve(double *A,int N,double *b,int *ipiv,double *x);

void minverse(double *A,int M,double *ipiv,double *inv);

#ifdef __cplusplus
}
#endif

#endif /* MATRIX_H_ */
