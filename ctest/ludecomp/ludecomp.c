#include "ludecomp.h"

static int sludecomp(double *A,int N) {
	int k,j,l,c1,c2;
	double ld,mult;
	for(k = 0; k < N-1; ++k) {
		c2 = k*N;
		ld = A[c2 + k];
		if (ld != 0.) {
			for (j = k+1; j < N; ++j) {
				c1 = j*N;
				mult = A[c1+k] /= ld;
				//printf("\n k %d j %d mult %lf \n",k,j,mult);
				for(l = k+1; l < N; ++l) {
					A[c1+l] -= mult * A[c2 + l];
				}
			}
		} else {
			printf("The Matrix is Singular");
			exit(-1);
		}
	}
	return 0;
}

int pludecomp(double *A,int N,int *ipiv) {
	int k,j,l,c1,c2,mind,tempi;
	double ld,mult,mval,temp;
	for(k=0;k < N;++k)
		ipiv[k] = k;
	
	for(k = 0; k < N-1; ++k) {
		//c2 = k*N;
		mval = fabs(A[k*N + k]);
		mind = k;
		for (j=k+1; j < N;++j) {
			if (mval < fabs(A[j*N + k])) {
				mval = A[j*N + k];
				mind = j;
			}
		}
		
		if ( mind != k) {
			c1 = k *N;
			c2 = mind * N;
			tempi = ipiv[mind];
			ipiv[mind] = ipiv[k];
			ipiv[k] = tempi;
			for (j = 0; j < N;j++) {
				temp = A[c1 + j];
				*(A + c1 + j) = *(A + c2 + j);
				*(A + c2 + j) = temp;
			}
		}
		c2 = k*N;
		ld = A[c2 + k];
		if (ld != 0.) {
			for (j = k+1; j < N; ++j) {
				c1 = j*N;
				mult = A[c1+k] /= ld;
				//printf("\n k %d j %d mult %lf \n",k,j,mult);
				for(l = k+1; l < N; ++l) {
					A[c1+l] -= mult * A[c2 + l];
				}
			}
		}
		
	}
	return 0;
	
}

int linsolve(double *A,int N,double *b,double *x,int *ipiv) {
	int i,j,c1,l;
	double *y;
	double sum,temp;
	
	y = (double*) malloc(sizeof(double) *N);
	/*
	 * Two step Solution L * U * x = b
	 * Let U*x = y
	 * Solve L * y = b for y (Forward Substitution
	 * Solve U * x = b for x (Back Substitution)
	 */ 
	for(i = 0; i < N;++i) {
		y[i] = 0.;
		x[i] = 0.;
		if ( A[i*N + i] == 0.) {
			printf("The Matrix system does not have a unique solution");
			exit(1);
		}
		//printf("\n B %d",ipiv[i]);
	}
	
	// Forward Substitution
	
	y[0] = b[ipiv[0]];
	for(i = 1; i < N; ++i) {
		sum = 0.;
		c1 = i*N;
		for(j = 0; j < i; ++j) {
			sum += y[j] * A[c1 + j];
		}
		y[i] = b[ipiv[i]] - sum;
	}
	
	// Back Substitution
	
	x[N - 1] = y[N - 1]/A[N * N - 1];
	
	for (i = N - 2; i >= 0; i--) {
		sum = 0.;
		c1 = i*(N+1);
		l=0;
		for(j = i+1; j < N;j++) {
			l++;
			sum += A[c1 + l] * x[j];
		}
		x[i] = (y[i] - sum) / A[c1];
	}
	
	free(y);
	return 0;
}