#include <stdio.h>
#include "cholesky.h"

void hilb(double *A,int N) {
	int i,j,t;
	for (i = 0; i < N;++i) {
		t = i*N;
		for (j = 0; j < N;++j) {
			A[t+j] = 1.0 /(i+j+1.0);
		}
	}
}

void fsub(double *A,int N,int M,double *bb,double *xx) {
	int i,j,t,k,u;
	double sum;
	double *x,*b;
	x = (double*) malloc(sizeof(double) * N);
	b = (double*) malloc(sizeof(double) * N);
	for(k = 0; k < M;++k) {
		u = k * N;
		for(i = 0; i < N;++i) {
			b[i] = bb[u+i];
			x[i] = 0.;
		}
		mdisplay(b,N,1);
		for (i = 0; i < N;++i) {
			t = i*N;
			sum = 0;
			for (j = 0; j < i;++j) {
				sum += A[t+j] * x[j];
				printf("sum %lf \n",sum);
			}
			x[i] = (b[i] - sum) / A[t+i];
		}
		for(i = 0; i < N;++i) {
			xx[u+i] = x[i];
		}
	
	}
	free(x);
	free(b);
}

int main(int argc, char **argv)
{
	int N,i,M,success;
	double *B,*C;
	N = 3;
	//A = (double*) malloc(sizeof(double) * N * N);
	B = (double*) malloc(sizeof(double) * N * N);
	C = (double*) malloc(sizeof(double) * N * N);
	double A[9] = {10,20,30,20,45,80,30,80,171};
	mdisplay(A,N,N);
	chold(A,N);
	mdisplay(A,N,N);
	/*
	hilb(A,N);
	mdisplay(A,N,N);
	success = bcholu(A,N);
	if (success == 0) { 
		mdisplay(A,N,N);
		mtranspose(A,N,N,B);
		mmult(B,A,C,N,N,N);
		mdisplay(C,N,N);
	} else {
		printf("The Input Matrix Is Not Positive Definite");
	}
	 */ 
	/*
	hilb(C,N);
	cholu(C,N);
	mdisplay(C,N,N);
	msub(A,C,B,N,N);
	mdisplay(B,N,N);
	 */ 
	 /*
	N = 3;
	M = 2;
	A = (double*) malloc(sizeof(double) * N * N);
	B = (double*) malloc(sizeof(double) * N * M);
	C = (double*) malloc(sizeof(double) * N * M);

	for (i = 0; i < N*N;++i) {
		A[i] = 0;
	}
	A[0] = 1; A[3] = 2; A[4] = 3; A[6] = 4; A[7] = 5; A[8] = 6;
	B[0] = 7; B[1] = 8; B[2] = 9;B[3] = 10; B[4] = 11; B[5] = 12;
	fsub(A,N,M,B,C);
	mdisplay(C,M,N);
	 */ 
	//free(A);
	free(B);
	free(C);
	return 0;
}
