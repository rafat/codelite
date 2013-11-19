#include "conjgrad.h"

static int richolu(double *A,int N, int stride, double *U22) {
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
			if (A[j] != 0) {
				A[j] /= u11;
			}
		}
		mmult(A+1,A+1,U22,N-1,1,N-1);
		for (i = 0; i < N-1; ++i) {
			u = stride + 1 + i * stride;
			w = i * (N-1);
			for(j = i; j < N-1;j++) {
				if (A[j + u] != 0) {
					A[j + u] -= U22[j + w];
				}
			}
		}
		
		sc = richolu(A+stride+1,N-1,stride,U22);
		if (sc == -1) {
			return -1;
		}
		
	}
	
	return sc;
	
}


static int icholu(double *A, int N) {
	int stride,i,j,t,sc;
	double *U22;
	U22 = (double*) malloc(sizeof(double) * N * N);
	stride = N; 
	
	sc = richolu(A,N,stride,U22);
	
	for(i=0; i < N;++i) {
		t = i *N;
		for(j=0;j < i;++j) {
			A[t+j] = 0.;
		}
	}

	free(U22);
	return sc;
	
}

int ichol(double *A, int N) {
	int sc;
	sc = icholu(A,N);
	return sc;
}

int stopcheck2(double fx,int N,double *xc,double *xf,double *jac,double *dx,double fsval,double gtol,double stol) {
	int rcode,i;
	double num,den;
	double stop0;
	double *scheck;
	
	rcode = 0;	
	
	scheck = (double*) malloc(sizeof(double) *N);
	
	if (fabs(fx) > fabs(fsval)) {
			den = fabs(fx);
	} else {
			den = fabs(fsval);
	}
	for(i = 0; i < N;++i) {
		if (fabs(xf[i]) > 1.0 / fabs(dx[i])) {
			num = fabs(xf[i]);
		} else {
			num = 1.0 / fabs(dx[i]);
		}
		scheck[i] = fabs(jac[i]) * num / den;
	}
	
	stop0 = array_max_abs(scheck,N);
	
	if (stop0 <= gtol) {
		rcode = 1;
	} else {
		for(i = 0; i < N;++i) {
			if (fabs(xf[i]) > 1.0 / fabs(dx[i])) {
				den = fabs(xf[i]);
			} else {
				den = 1.0 / fabs(dx[i]);
			}
			num = fabs(xf[i] - xc[i]);
			scheck[i] = num / den;
		}
		stop0 = array_max_abs(scheck,N);
		if (stop0 <= stol) {
			rcode = 2;
		}
	}
	
	free(scheck);
	return rcode;
}

int cgr(double *xi,int N,double *A,double *b,double *xf) {
	int i,rcode;
	int iter,siter;
	double *temp,*rk,*pk,*dx,*jac,*apk;
	double stol,gtol,bn;
	
	temp = (double*) malloc(sizeof(double) *8);
	rk = (double*) malloc(sizeof(double) *N);
	pk = (double*) malloc(sizeof(double) *N);
	apk = (double*) malloc(sizeof(double) *N);
	dx = (double*) malloc(sizeof(double) *N);
	jac = (double*) malloc(sizeof(double) *N);
	
	iter = 0;
	rcode = 0;
	siter = (int) SETITER;
	gtol = pow((double) FDVAL,1.0/3.0);
	stol = gtol * gtol;
	

	mmult(A,xi,rk,N,N,1);
	msub(rk,b,rk,N,1);
	for(i = 0; i < N;++i) {
		pk[i] = -rk[i];
		xf[i] = xi[i];
	}
		
	bn = l2norm(b,N);
	
	while (rcode == 0 && iter < siter) {
		iter++;
		
		mmult(rk,rk,temp,1,N,1);
		mmult(A,pk,apk,N,N,1);
		mmult(pk,apk,temp+1,1,N,1);
		temp[2] = temp[0]/temp[1]; // alpha
		
		for(i = 0; i < N;++i) {
			xf[i] = xi[i] + temp[2] * pk[i];
			apk[i] *= temp[2]; 
		}
		
		madd(rk,apk,rk,1,N);
		mmult(rk,rk,temp+3,1,N,1);
		if (sqrt(temp[3]) <= stol * bn) {
			rcode = 1;
		}
		temp[4] = temp[3] / temp[0];
		scale(pk,N,1,temp[4]);
		mdisplay(xf,1,N);
		msub(pk,rk,pk,N,1);
		
		printf("%d \n",iter);
		
		for(i = 0; i < N;++i) {
			xi[i] = xf[i];
		}
	}
	
	free(temp);
	free(rk);
	free(pk);
	free(dx);
	free(jac);
	free(apk);
	
	return rcode;
}

int cgpc(double *xi,int N,double *A,double *b,double *xf) {
	int i,j,t,rcode;
	int iter,siter;
	double *temp,*rk,*pk,*dx,*jac,*apk,*U,*y,*AB;
	double stol,gtol,bn;
	
	temp = (double*) malloc(sizeof(double) *8);
	rk = (double*) malloc(sizeof(double) *N);
	pk = (double*) malloc(sizeof(double) *N);
	apk = (double*) malloc(sizeof(double) *N);
	dx = (double*) malloc(sizeof(double) *N);
	jac = (double*) malloc(sizeof(double) *N);
	U = (double*) malloc(sizeof(double) *N * N);
	AB = (double*) malloc(sizeof(double) *N * N);
	y = (double*) malloc(sizeof(double) *N);
	
	iter = 0;
	rcode = 0;
	siter = (int) SETITER;
	gtol = pow((double) FDVAL,1.0/3.0);
	stol = gtol * gtol;
	

	mmult(A,xi,rk,N,N,1);
	msub(rk,b,rk,N,1);
	for(i = 0; i < N;++i) {
		t = i *N;
		for(j = 0; j < N;++j) {
			AB[t+j] = A[t+j];
		}
	}
	ichol(AB,N);
	mdisplay(AB,N,N);
	mtranspose(AB,N,N,U);
	linsolve_lower(U,N,rk,y);
	
	for(i = 0; i < N;++i) {
		pk[i] = -y[i];
		xf[i] = xi[i];
	}
		
	bn = l2norm(b,N);
	
	while (rcode == 0 && iter < siter) {
		iter++;
		
		mmult(rk,y,temp,1,N,1);
		mmult(A,pk,apk,N,N,1);
		mmult(pk,apk,temp+1,1,N,1);
		temp[2] = temp[0]/temp[1]; // alpha
		
		for(i = 0; i < N;++i) {
			xf[i] = xi[i] + temp[2] * pk[i];
			apk[i] *= temp[2]; 
		}
		
		madd(rk,apk,rk,1,N);
		linsolve_lower(U,N,rk,y);
		mmult(rk,y,temp+3,1,N,1);
		if (sqrt(temp[3]) <= stol * bn) {
			rcode = 1;
		}
		temp[4] = temp[3] / temp[0];
		scale(pk,N,1,temp[4]);
		mdisplay(xf,1,N);
		msub(pk,y,pk,N,1);
		
		printf("%d \n",iter);
		
		for(i = 0; i < N;++i) {
			xi[i] = xf[i];
		}
	}
	
	free(temp);
	free(rk);
	free(pk);
	free(dx);
	free(jac);
	free(apk);
	free(U);
	free(y);
	free(AB);
	
	return rcode;
}

int cgp(double (*funcpt)(double *,int),double *xi,int N,double *A,double *b,double *xf) {
	int i,rcode;
	int iter,siter;
	double *temp,*rk,*pk,*dx,*jac,*apk;
	double fsval,stol,gtol,fxf;
	
	temp = (double*) malloc(sizeof(double) *8);
	rk = (double*) malloc(sizeof(double) *N);
	pk = (double*) malloc(sizeof(double) *N);
	apk = (double*) malloc(sizeof(double) *N);
	dx = (double*) malloc(sizeof(double) *N);
	jac = (double*) malloc(sizeof(double) *N);
	
	iter = 0;
	rcode = 0;
	siter = (int) SETITER;
	gtol = pow((double) FDVAL,1.0/3.0);
	stol = gtol * gtol;
	
	// Values that may not be needed
	
	fsval = 1.0;
	for(i = 0; i < N;++i) {
		dx[i] = 1.0;
	}
	
	mmult(A,xi,rk,N,N,1);
	msub(rk,b,rk,N,1);
	for(i = 0; i < N;++i) {
		pk[i] = -rk[i];
		xf[i] = xi[i];
	}
		
	
	while (rcode == 0 && iter < siter) {
		iter++;
		
		mmult(rk,rk,temp,1,N,1);
		mmult(A,pk,apk,N,N,1);
		mmult(pk,apk,temp+1,1,N,1);
		temp[2] = temp[0]/temp[1]; // alpha
		
		for(i = 0; i < N;++i) {
			xf[i] = xi[i] + temp[2] * pk[i];
			apk[i] *= temp[2]; 
		}
		
		madd(rk,apk,rk,1,N);
		mmult(rk,rk,temp+3,1,N,1);
		
		temp[4] = temp[3] / temp[0];
		scale(pk,N,1,temp[4]);
		mdisplay(xf,1,N);
		msub(pk,rk,pk,N,1);
		
		fxf = funcpt(xi,N);
		printf("%d %g \n",iter,fxf);
		if (fxf >= DBL_MAX || fxf <= -DBL_MAX) {
			printf("Program Exiting as the function value exceeds the maximum double value");
			exit(1);
		}
		grad_fd(funcpt,xi,N,dx,jac);
		rcode = stopcheck2(fxf,N,xi,xf,jac,dx,fsval,gtol,stol);
		for(i = 0; i < N;++i) {
			xi[i] = xf[i];
		}
	}
	
	free(temp);
	free(rk);
	free(pk);
	free(dx);
	free(jac);
	free(apk);
	
	return rcode;
}

int conjgrad_min_lin(double (*funcpt)(double *,int),double *xi,int N,double *xf) {
	int rcode,i;
	double *A,*b,*dx;
	
	b = (double*) malloc(sizeof(double) *N);
	A = (double*) malloc(sizeof(double) *N * N);
	dx = (double*) malloc(sizeof(double) *N);
	
	for(i = 0; i < N;++i) {
		dx[i] = 1.0;
	}
	
	//hessian_fd(funcpt,xi,N,dx,A);
	//grad_fd(funcpt,xi,N,dx,b);
	A[0] = 4; A[1] = 1; A[2] = 1; A[3] = 3;
	b[0] = 1; b[1] = 2;
	xi[0] = 2; xi[1] = 1;
	mdisplay(A,N,N);
	
	rcode = cgpc(xi,N,A,b,xf);
	
	free(A);
	free(b);
	free(dx);
	return rcode;
}