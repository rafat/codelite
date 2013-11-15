
#include "cholesky.h"

static int rcholu(double *A,int N, int stride, double *U22) {
	int sc;
	int j,i,u,w;
	double u11;
	
	if (N == 1) {
		if (A[0] > 0) {
			A[0] = sqrt(A[0]);
			return 0;
		} else {
			return -1;
		}
	} else {
		if (A[0] < 0) {
			return -1;
		}
		u11 = sqrt(A[0]);
		A[0] = u11;
		for (j = 1; j < N;++j) {
			A[j] /= u11;
		}
		mmult(A+1,A+1,U22,N-1,1,N-1);
		for (i = 0; i < N-1; ++i) {
			u = stride + 1+ i * stride;
			w = i * (N-1);
			for(j = i; j < N-1;j++) {
				A[j + u] -= U22[j + w];
			}
		}
		
		sc = rcholu(A+stride+1,N-1,stride,U22);
		if (sc == -1) {
			return -1;
		}
		
	}
	
	return sc;
	
}

static int rbcholu(double *A,int N, int stride, double *UB, double *UT) {
	int bs,bb,i,j,Nb,t,k,u,v,w,sc;
	double *b,*x,*U12,*U12T;
	double sum;
	
	bs = (int) BLOCKSIZE;
	bb = bs*bs;
	
	if (N <= BLOCKSIZE) {
		sc = rcholu(A,N,stride,UB);
		if (sc == -1) {
			return -1;
		}
	} else {
		Nb = N - bs;
		x = (double*) malloc(sizeof(double) * bs);
		b = (double*) malloc(sizeof(double) * bs);
		U12T = (double*) malloc(sizeof(double) * Nb * bs);
		U12 = (double*) malloc(sizeof(double) * Nb * bs);
		rcholu(A,bs,stride,UB); // U11
		
		for (i =0; i < bs;++i) {
			t = i *stride;
			u = 0;
			for(j = 0; j < N;++j) {
				UT[u+i] = A[j+t];
				u += bs;
			}
		}
		
		for(k = 0; k < Nb;++k) {
			u = k * bs;
			for(i = 0; i < bs;++i) {
				b[i] = UT[bb+u+i];
				x[i] = 0.;
			}
			for (i = 0; i < bs;++i) {
				t = i*bs;
				sum = 0;
				for (j = 0; j < i;++j) {
					sum += UT[t+j] * x[j];
				}
				x[i] = (b[i] - sum) / UT[t+i];
			}
			v = bs + k;
			for(i = 0; i < bs;++i) {
				A[v] = x[i];
				U12T[u+i] = x[i];
				v += stride;
			}
		}
		
		mtranspose(U12T,Nb,bs,U12);
		mmult(U12T,U12,UT,Nb,bs,Nb);
		free(U12T);
		free(U12);
		free(b);
		free(x);
		for (i = 0; i < Nb; ++i) {
			u = bs * stride + bs + i * stride;
			w = i * Nb;
			for(j = i; j < Nb;j++) {
				A[j + u] -= UT[j + w];
			}
		}
		
		sc = rbcholu(A + bs * stride + bs,Nb,stride,UB,UT);
		if (sc == -1) {
			return -1;
		}
	}
	
	return sc;
}

int cholu(double *A, int N) {
	int stride,i,j,t,sc;
	double *U22;
	U22 = (double*) malloc(sizeof(double) * N * N);
	stride = N; 
	
	sc = rcholu(A,N,stride,U22);
	
	for(i=0; i < N;++i) {
		t = i *N;
		for(j=0;j < i;++j) {
			A[t+j] = 0.;
		}
	}

	free(U22);
	return sc;
	
}

int bcholu(double *A, int N) {
	int stride,i,j,t,b,sc;
	double *UB,*UT;
	b = (int) BLOCKSIZE;
	UT = (double*) malloc(sizeof(double) * N * N);
	UB = (double*) malloc(sizeof(double) * b * b);
	stride = N; 
	
	sc = rbcholu(A,N,stride,UB,UT);
	
	for(i=0; i < N;++i) {
		t = i *N;
		for(j=0;j < i;++j) {
			A[t+j] = 0.;
		}
	}

	free(UB);
	free(UT);
	
	return sc;
	
}

int chol(double *A, int N) {
	int sc;
	sc = cholu(A,N);
	return sc;
}

static void rchold(double *A,int N, int stride, double *U22) {
	int j,i,u,w;
	double d1;
	
	if (N == 1) {
		return;
	} else {
		d1 = A[0];
		for (j = 1; j < N;++j) {
			A[j] /= d1;
		}
		mmult(A+1,A+1,U22,N-1,1,N-1);
		scale(U22,N-1,N-1,d1);
		for (i = 0; i < N-1; ++i) {
			u = stride + 1+ i * stride;
			w = i * (N-1);
			for(j = i; j < N-1;j++) {
				A[j + u] -= U22[j + w];
			}
		}
		
		rchold(A+stride+1,N-1,stride,U22);
	
	}
		
}

void chold(double *A, int N) {
	int stride,i,j,t,sc;
	double *U22;
	U22 = (double*) malloc(sizeof(double) * N * N);
	stride = N; 
	
	rchold(A,N,stride,U22);
	
	for(i=0; i < N;++i) {
		t = i *N;
		for(j=0;j < i;++j) {
			A[t+j] = 0.;
		}
	}

	free(U22);
	
}