#include "qrdecomp.h"

void eye(double *mat,int N) {
	int i,j,t;
	for(i = 0;i < N;++i) {
		for(j =0; j < N;++j) {
			t = i*N;
			if (i == j) {
				mat[t+j] = 1.;
			} else {
				mat[t+j] = 0.;
			}
		}
		
	}
}

static double house_1(double*x,int N,double *v) {
	double beta,mu,temp;
	double *sigma;
	int i;
	
	sigma = (double*) malloc(sizeof(double) * 1);
	
	if (N > 1) {
		mmult(x+1,x+1,sigma,1,N-1,1);
	} else {
		sigma[0] = 0.0;
	}
	
	v[0] =1.;
	for (i = 1; i < N;++i) {
		v[i] = x[i];
	}
	
	if (sigma[0] == 0. && x[0] >= 0.) {
		beta = 0.;
	} else if (sigma[0] == 0. && x[0] < 0.) {
		beta = -2.;
	}else {
		mu = sqrt(sigma[0] + x[0] * x[0]);
		
		if (x[0] <= 0.) {
			v[0] = x[0] - mu;
		} else {
			v[0] = - sigma[0] / (x[0] + mu);
		}
		temp = v[0];
		
		beta = (2.0 * v[0] * v[0]) /(sigma[0] + v[0] * v[0]);
		
		for (i = 0; i < N;++i) {
			v[i] /= temp;
		}
		
	}
	
	free(sigma);
	return beta;
}

static double house_2(double*x,int N,double *v) {
	double sgn,beta,sc;
	double *sigma,*e;
	int i;
	
	sigma = (double*) malloc(sizeof(double) * 1);
	e = (double*) malloc(sizeof(double) * N);
	
	beta = 2.0;
	mmult(x,x,sigma,1,N,1);
	sigma[0] = sqrt(sigma[0]);
	
	e[0] =1.;
	for (i = 1; i < N;++i) {
		e[i] = 0.;
	}
	
	if (x[0] > 0.) {
		sgn = 1.0;
	} else if (x[0] < 0.) {
		sgn = -1.0;
	} else if (x[0] == 0.) {
		sgn = 0.;
	}
	
	sc = sigma[0] * sgn;
	
	//scale(e,N,1,sc);
	
	e[0] *= sc;
	
	for(i = 0; i < N;++i) {
		v[i] = e[i] + x[i];
	}
	
	mmult(v,v,sigma,1,N,1);
	sigma[0] = sqrt(sigma[0]);
	
	for(i = 0; i < N;++i) {
		v[i] = v[i] / sigma[0];
	}
	
	free(sigma);
	free(e);
	return beta;
}

double house(double*x,int N,double *v) {
	double beta;
	beta = house_1(x,N,v);
	return beta;
}


void housemat(double *v, int N,double beta,double *mat) {
	double *temp;
	
	temp = (double*) malloc(sizeof(double) * N * N);
	eye(mat,N);
	mmult(v,v,temp,N,1,N);
	scale(temp,N,N,beta);
	msub(mat,temp,mat,N,N);
	
	free(temp);
}

void qr_house(double *A, int M, int N,double *bvec) {
	int j,i,k,u,t;
	double *x,*v,*AT,*w;
	double beta;
	
	if (M < N) {
			printf("M should be greater than or equal to N");
			exit(1);
	}
	x = (double*) malloc(sizeof(double) * M);
	v = (double*) malloc(sizeof(double) * M);
	AT = (double*) malloc(sizeof(double) * M * N);
	w = (double*) malloc(sizeof(double) * M * M);
	
	for(j = 0; j < N;++j) {
		for(i=j;i < M;++i) {
			x[i-j] = A[i*N+j];
			
		}
		
		beta = house(x,M-j,v);
		bvec[j] = beta;
	
		for (i=j; i < M; i++) {
			t = i * N;
			u = 0;
			for (k=j; k < N; k++) {
				AT[u+i-j] = A[k+t];
				u+=(M-j);
				
			}
			
		}
		
		
		mmult(AT,v,w,N-j,M-j,1);
		scale(w,N-j,1,beta);
		mmult(v,w,AT,M-j,1,N-j);
		for (i=j; i < M; i++) {
			t = i *N;
			for (k=j; k < N; k++) {
				A[t+k] -= AT[(i-j)*(N-j) + k - j];
			}
		}
		if (j < M) {
			
			for(i=j+1;i < M;++i) {
				A[i*N+j] = v[i-j];
			}
		}
		 
	}
	
	free(x);
	free(v);
	free(AT);
	free(w);
	
}

void getQR_house(double *A,int M,int N,double *bvec,double *Q, double *R) {
	int i,j,k,t,u;
	double *x,*v,*AT,*w;
	
	x = (double*) malloc(sizeof(double) * M);
	v = (double*) malloc(sizeof(double) * M);
	AT = (double*) malloc(sizeof(double) * M * N);
	w = (double*) malloc(sizeof(double) * M * M);
	
	for(i = 0; i < N;++i) {
		t = i *N;
		for(j = 0; j < N;++j) {
			if (i > j) {
				R[t+j] = 0.;
			} else {
				R[t+j] = A[t+j];
			}
		}
	}
	
	for(i = 0; i < M;++i) {
		t = i *N;
		for(j = 0; j < N;++j) {
			if (i == j) {
				Q[t+j] = 1.;
			} else {
				Q[t+j] = 0.;
			}
		}
	}
	
	for(j = N-1; j >= 0;--j) {
		v[0] = 1.;
		for(i=j+1;i < M;++i) {
			v[i-j] = A[i*N+j];
			
		}
		
		for (i=j; i < M; i++) {
			t = i * N;
			u = 0;
			for (k=j; k < N; k++) {
				AT[u+i-j] = Q[k+t];
				u+=(M-j);
			}
			
		}
	
		mmult(AT,v,w,N-j,M-j,1);
		scale(w,N-j,1,bvec[j]);
		mmult(v,w,AT,M-j,1,N-j);
		
		for (i=j; i < M; i++) {
			t = i *N;
			for (k=j; k < N; k++) {
				Q[t+k] -= AT[(i-j)*(N-j) + k - j];
			}
		}
	 
	}
	
	free(x);
	free(v);
	free(AT);
	free(w);
}

void hessenberg(double *A,int N) {
	int k,i,j,t,u;
	double *x,*v,*AT,*w;
	double beta;
	x = (double*) malloc(sizeof(double) * N);
	v = (double*) malloc(sizeof(double) * N);
	AT = (double*) malloc(sizeof(double) * N * N);
	w = (double*) malloc(sizeof(double) * N);
	
	for (k = 0; k < N-2;++k) {
		for(i=k + 1;i < N;++i) {
			x[i-k-1] = A[i*N+k];
			//printf("x %lf \n",x[i-k-1]);
			
		}
		
		beta = house(x,N-k-1,v);
		
		
		for (i=k+1; i < N; i++) {
			t = i * N;
			u = 0;
			for (j=k; j < N; j++) {
				AT[u+i-k-1] = A[j+t];
				u+=(N-k-1);				
			}
		}
		//mdisplay(AT,N-k,N-k-1);
		
		
		mmult(AT,v,w,N-k,N-k-1,1);
		scale(w,N-k,1,beta);
		mmult(v,w,AT,N-k-1,1,N-k);
		//mdisplay(AT,N-k-1,N-k);
		
		for (i=k+1; i < N; i++) {
			t = i * N;
			for (j=k; j < N; j++) {
				A[t+j] -= AT[(i-k-1)*(N-k) + j - k];
			}
		}
		//mdisplay(A,N,N);
		
		for (i=0; i < N; i++) {
			t = i * N;
			u = i * (N-k-1);
			for (j=k+1; j < N; j++) {
				AT[u+j-k-1] = A[t+j];
			}
		}
		//mdisplay(AT,N,N-k-1);
		
		mmult(AT,v,w,N,N-k-1,1);
		scale(w,N,1,beta);
		mmult(w,v,AT,N,1,N-k-1);
		//mdisplay(AT,N,N-k-1);
		
		for (i=0; i < N; i++) {
			t = i * N;
			u = i * (N-k-1);
			for (j=k+1; j < N; j++) {
				A[t+j] -= AT[u+j-k-1];
			}
		}
		
	}
	
	free(x);
	free(v);
	free(AT);
	free(w);
}

void francisQR(double *A,int N) {
	int m,n,k,q,r,t,u,i,j;
	double s,t2,beta;
	double *x,*v,*AT,*w;
	int NN;
	/*
	 * Reference - Algorithm 7.5.1 Golub,van Loan Matrix Computations 3rd Edition
	 */ 
	x = (double*) malloc(sizeof(double) * 3);
	v = (double*) malloc(sizeof(double) * 3);
	AT = (double*) malloc(sizeof(double) * 3 * N);
	w = (double*) malloc(sizeof(double) * N);
	n = N-1;
	m = n-1;
	NN = N*N;
	
	s = A[NN-1] + A[NN-N-2];
	t2 = A[NN-1] * A[NN-N-2] - A[NN-2] * A[NN-N-1];
	
	x[0] = A[0]*A[0] + A[1]*A[N] - s*A[0] + t2;
	x[1] = A[N]*(A[0] + A[N+1] - s);
	x[2] = A[N] * A[N+N+1];
	if (N <= 2) {
		return;
	}
	
	for (k = -1; k < N - 3;++k) {
		
		beta = house(x,3,v);
		//mdisplay(x,3,1);
		if (k > 0) {
			q = k;
		} else {
			q = 0;
		}
		
		//printf("q %d \n",q);
		for (i=k+1; i < k+4; i++) {
			t = i * N;
			u = 0;
			for (j=q; j < N; j++) {
				AT[u+i-k-1] = A[j+t];
				u+=3;		
			}
		}
		
		mmult(AT,v,w,N-q,3,1);
		scale(w,N-q,1,beta);
		mmult(v,w,AT,3,1,N-q);
		
		for (i=k+1; i < k+4; i++) {
			t = i * N;
			for (j=q; j < N; j++) {
				A[t+j] -= AT[(i-k-1)*(N-q) + j - q];
			}
		}
		//mdisplay(A,N,N);
		if (k+4 >= n) {
			r = N;
		} else {
			r = k+4+1;
		}
		//printf("r %d \n",r);
		for (i=0; i < r; i++) {
			t = i * N;
			u = i * 3;
			for (j=k+1; j < k+4; j++) {
				AT[u+j-k-1] = A[t+j];
			}
		}
		
		mmult(AT,v,w,r,3,1);
		scale(w,r,1,beta);
		mmult(w,v,AT,r,1,3);
		//mdisplay(AT,N,N-k-1);
		
		for (i=0; i < r; i++) {
			t = i * N;
			u = i * 3;
			for (j=k+1; j < k+4; j++) {
				A[t+j] -= AT[u+j-k-1];
			}
		}
		//mdisplay(A,N,N);
		x[0] = A[N*(k+2) + k+1];
		x[1] = A[N*(k+3) + k+1];
		
		if (k < n-3) {
			x[2] = A[N*(k+4) + k+1];
		} 
		//mdisplay(x,3,1);
		
	}
	//mdisplay(x,2,1);
	beta = house(x,2,v);
	
	for (i=n-1; i < N; i++) {
		t = i * N;
		u = 0;
		for (j=n-2; j < N; j++) {
			AT[u+i-n+1] = A[j+t];
			u+=2;		
		}
	}
	
	mmult(AT,v,w,3,2,1);
	scale(w,3,1,beta);
	mmult(v,w,AT,2,1,3);
	for (i=n-1; i < N; i++) {
		t = i * N;
		for (j=n-2; j < N; j++) {
			A[t+j] -= AT[(i-n+1)*3 + j - n + 2];
		}
	}
	
	for (i=0; i < N; i++) {
		t = i * N;
		u = i * 2;
		for (j=n-1; j < N; j++) {
			AT[u+j-n+1] = A[t+j];
		}
	}
	
	mmult(AT,v,w,N,2,1);
	scale(w,N,1,beta);
	mmult(w,v,AT,N,1,2);
		//mdisplay(AT,N,N-k-1);
		
	for (i=0; i < N; i++) {
		t = i * N;
		u = i * 2;
		for (j=n-1; j < N; j++) {
			A[t+j] -= AT[u+j-n+1];
		}
	}
	
	
	free(x);
	free(v);
	free(AT);
	free(w);
	
}

void eig22(double *A, int stride,double *eigre,double *eigim) {
	int N;
	double a11,a12,a21,a22,c,s,c2,s2,cs,t1,t,t2,at11,at12,at21,at22;
	N = stride;
	
	a11 = A[0];
	a12 = A[1];
	a21 = A[N];
	a22 = A[N+1];
	
	if ( (a12 + a21) == 0) {
		c = 1./sqrt(2.0);
		s = c;
	} else {
		t1 = (a11 - a22) / (a12 + a21);
		t = t1 /(1. + sqrt(1+t1*t1));
		c = 1./sqrt(1 + t*t);
		s = c*t;
	}
	
	c2 = c*c;
	s2 = s*s;
	cs = c*s;

	at11 = c2 * a11 + s2 * a22 - cs * (a12 + a21);
	at12 = c2 * a12 - s2 * a21 + cs * (a11 - a22);
	at21 = c2 * a21 - s2 * a12 + cs * (a11 - a22);
	at22 = c2 * a22 + s2 * a11 + cs * (a12 + a21);
	
	eigre[0] = eigre[1] = at11;
	eigim[0] = sqrt(-at12 * at21);
	eigim[1] = -sqrt(-at12 * at21);
	
	if ( at12*at21 >= 0) {
		if (at12 == 0) {
			c = 0;
			s = 1;
			c2 = 0;
			s2 = 1;
			cs = 0;
		} else {
			t = sqrt(at21/at12);
			t2 = t * t;
			cs = t/(1+t2);
			c2 = (1+t2);
			s2 = t2 /(1+t2);
		}
		eigim[0] = eigim[1] = 0.0;
		eigre[0] = at11 - cs * (at12 + at21);
		eigre[1] = at11 + cs * (at12 + at21);
		
	}
	
}

int francis_iter(double *A, int N, double *H) {
	int success,brkpoint;
	int i,j,it,p,q,t,u;
	double *temp;
	success = 0;
	brkpoint = 30 * N;
	it = 0;
	p = N - 1;
	temp = (double*) malloc(sizeof(double) * N * N);
	for(i = 0; i < N*N;++i) {
		H[i] = A[i];
	}
	
	hessenberg(H,N);
	
	while (p > 1 && it < brkpoint) {
		
		while (p > 1 && (H[N*p + p-1] == 0 || H[N*(p-1) + p-2] == 0)) {
			if (H[N*p + p-1] == 0) {
				p--;
			} else if (H[N*(p-1) + p-2] == 0) {
				p=p-2;
			}
		}
		
		if (p > 0) {
			q = p-1;
			while (q > 0 && fabs(H[N*q + q-1]) != 0) {
				q--;
			}
			//printf("%d %d \n",q,p);
			for (i=q; i <= p; i++) {
				t = i * N;
				u = (i-q) * (p-q+1);
				for (j=q; j <= p; j++) {
					temp[u+j-q] = H[t+j];
				}
			}
			francisQR(temp,p-q+1);
			for (i=q; i <= p; i++) {
				t = i * N;
				u = (i-q) * (p-q+1);
				for (j=q; j <= p; j++) {
					H[t+j] = temp[u+j-q];
				}
			}
			//mdisplay(H,N,N);
			for(i = q; i <= p-1;++i) {
				if ( fabs(H[(i+1)*N+i]) <= TOL * (fabs(H[i*N+i]) + fabs(H[(i+1)*N+i+1]) ) ) {
					H[(i+1)*N+i] = 0.;
				}
			}
			it++;
			//printf("iter %d \n",it);
		}
	}
	
	if (it == brkpoint) {
		success = 0;
	} else {
		success = 1;
	}
	
	free(temp);
	return success;
}

static void eig2t(double *A, int stride) {
	int N;
	double a11,a12,a21,a22,c,s,c2,s2,cs,t1,t,t2,at11,at12,at21,at22;
	N = stride;
	
	a11 = A[0];
	a12 = A[1];
	a21 = A[N];
	a22 = A[N+1];
	
	if ( (a12 + a21) == 0) {
		c = 1./sqrt(2.0);
		s = c;
	} else {
		t1 = (a11 - a22) / (a12 + a21);
		t = t1 /(1. + sqrt(1+t1*t1));
		c = 1./sqrt(1 + t*t);
		s = c*t;
	}
	
	c2 = c*c;
	s2 = s*s;
	cs = c*s;

	at11 = c2 * a11 + s2 * a22 - cs * (a12 + a21);
	at12 = c2 * a12 - s2 * a21 + cs * (a11 - a22);
	at21 = c2 * a21 - s2 * a12 + cs * (a11 - a22);
	at22 = c2 * a22 + s2 * a11 + cs * (a12 + a21);
	A[0] = at11;
	A[1] = at12;
	A[N] = at21;
	A[N+1] = at22;
	/*
	eigre[0] = eigre[1] = at11;
	eigim[0] = sqrt(-at12 * at21);
	eigim[1] = -sqrt(-at12 * at21);
	
	if ( at12*at21 >= 0) {
		if (at12 == 0) {
			c = 0;
			s = 1;
			c2 = 0;
			s2 = 1;
			cs = 0;
		} else {
			t = sqrt(at21/at12);
			t2 = t * t;
			cs = t/(1+t2);
			c2 = (1+t2);
			s2 = t2 /(1+t2);
		}
		eigim[0] = eigim[1] = 0.0;
		eigre[0] = at11 - cs * (at12 + at21);
		eigre[1] = at11 + cs * (at12 + at21);
		
	}
*/	
}

void eig(double *A,int N,double *eigre,double *eigim) {
	int i,t,u,n;
	double *H;
	double t1,t2,cs;
	H = (double*) malloc(sizeof(double) * N * N);
	n = N - 1;
	francis_iter(A,N,H);
	//mdisplay(H,N,N);
	i = 0;
	while (i < n) {
		u = i * N;
		t = (i+1)*N;
		if (H[t+i] != 0.) {
			eig2t(H+u+i,N);
			i = i +2;
		} else {
			i++;
		}
		
	}
	//mdisplay(H,N,N);
	i = 0;
	while (i < n) {
		u = i * N;
		t = (i+1)*N;
		
		if (H[t+i] != 0.) {
			if (H[u+i+1] * H[t+i] < 0.) {
				eigre[i] = H[u+i];
				eigre[i+1] = H[t+i+1];
				eigim[i] = sqrt(-H[u+i+1] * H[t+i]);
				eigim[i+1] = -sqrt(-H[u+i+1] * H[t+i]);
			} else {
				if (H[u+i+1] == 0.) {
					cs = 0.;
				} else {
					t1 = sqrt(H[t+i]/H[u+i+1]);
					t2 = t1 * t1;
					cs = t1/(1+t2);
				}
				eigre[i] = H[u+i] - cs * (H[u+i+1] + H[t+i]);
				eigre[i+1] = H[u+i] + cs * (H[u+i+1] + H[t+i]);
				eigim[i] = 0.;
				eigim[i+1] = 0.;
				
			}
			
			i= i + 2;
			
		} else {
			eigre[i] = H[u+i];
			eigim[i] = 0.;
			i++;
		}
		
	}
	
	if (i == n) {
		eigre[i] = H[N*N - 1];
		eigim[i] = 0.;
	}
	
	free(H);
}