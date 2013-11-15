#include "newton1d.h"

double func1(double *x,int N) {
	double f;
	f = x[0]*x[0] + x[1]*x[1] - 2;	
	return f;
}

double func2(double *x,int N) {
	double f;
	f = pow((double) EPSILON,x[0] - 1.0) + x[1]*x[1]*x[1] - 2.0;
	return f;
}

double func3(double *x,int N, double *jac,double *hess) {
	double f;
	f = pow((x[0]-2.0),4.0) + pow((x[0]-2.0),2.0) * x[1]*x[1] + (x[1] + 1.0) * (x[1] + 1.0);
	
	//Jacobian
	
	jac[0] = 4.0 * pow((x[0]-2.0),3.0) + 2.0 * (x[0]-2.0) * x[1]*x[1];
	jac[1] = 2.0 * pow((x[0]-2.0),2.0) * x[1] + 2.0 * (x[1] + 1.0);
	
	//Hessian
	
	hess[0] = 12.0 * pow((x[0]-2.0),2.0) + 2.0 * x[1]*x[1];
	hess[3] = 2.0 * pow((x[0]-2.0),2.0) + 2.0;
	hess[1] = 4.0 * (x[0]-2.0) * x[1];
	hess[2] = hess[1];
	
	return f;
}

double func4(double *x,int N) {
	double f;
	f = pow((x[0]-2.0),4.0) + pow((x[0]-2.0),2.0) * x[1]*x[1] + (x[1] + 1.0) * (x[1] + 1.0);
	
	return f;
}

double solvefunc(double (*funcpt)(double *,int),double *x,int N) {
	double result;
	
	result = funcpt(x,N);
	
	return result;
}
/*
static void rcholgmw(double *A,int N, int stride, double *U22) {
	int j,i,u,w;
	double d1;
	
	if (N == 1) {
		A[0] = sqrt(A[0]);
		return;
	} else {
		d1 = sqrt(A[0]);
		A[0] = d1;
		for (j = 1; j < N;++j) {
			A[j] /= d1;
		}
		mmult(A+1,A+1,U22,N-1,1,N-1);
		for (i = 0; i < N-1; ++i) {
			u = stride + 1+ i * stride;
			w = i * (N-1);
			for(j = i; j < N-1;j++) {
				A[j + u] -= U22[j + w];
			}
		}
		
		rcholgmw(A+stride+1,N-1,stride,U22);
	
	}
		
}
*/

double cholmod(double *A, int N, double *L, double maxinp) {
	/*
	 * Algorithm 5.5.2 Dennis, Schnabel, Numerical Methods for Unconstrained Optimization and 
	 * Non-Linear Equations
	 */ 
	int i,j,k,step,step2;
	double beta,d1,ls;
	double *U22;
	double maxadd,maxdiag,minl,minl2,minljj;
	
	U22 = (double*) malloc(sizeof(double) * N * N);
	
	minl = sqrt(sqrt((double) FDVAL)) * maxinp;
	maxadd = 0.0;
	
	for(i = 0; i < N;++i) {
		U22[i] = A[i+i*N];
	}
	
	for(i = 0; i < N * N;++i) {
		L[i] = 0.0;
	}
	
	maxdiag = array_max_abs(U22,N);
	
	if ( maxinp == 0.0) {
		maxinp = sqrt(maxdiag);
	}
	
	minl2 = sqrt((double) FDVAL) * maxinp;
	
	for(j = 0; j < N;++j) {
		step = j * N;
		ls = 0.0;
		for(i = 0; i < j;++i) {
			ls += L[step+i] * L[step+i];
		}
		
		L[step+j] = A[step+j] - ls;
		
		for(i = j+1; i < N;++i) {
			ls = 0.0;
			step2 = i * N;
			for(k = 0;k < j;++k) {
				ls += L[step+k] * L[step2+k];
			}
			L[step2+j] = A[step+i] - ls;
			
			if (fabs(L[step2+j]) > minljj) {
				minljj = fabs(L[step2+j]);
			}
		}
		
		if (minljj/maxinp > minl) {
			minljj = minljj/maxinp;
		} else {
			minljj = minl;
		}
		
		if (L[step+j] > minljj*minljj) {
				L[step+j] = sqrt(L[step+j]);
		} else {
			if (minljj < minl2) {
				minljj = minl2;
			}
			
			if (maxadd < (minljj*minljj - L[step+j])) {
				maxadd = minljj*minljj - L[step+j];
			}
			
			L[step+j] = minljj;
		}
		
		for(i = j+1; i < N;++i) {
			L[i*N+j] /= L[step+j];
		}
		
		
	}
	
	free(U22);
	
	return maxadd;
}

double modelhess2(double *A,int N,double *L) {
	/*
	 * Algorithm 5.5.1 Dennis, Schnabel, Numerical Methods for Unconstrained Optimization and 
	 * Non-Linear Equations
	 */ 
	 double *U22;
	 double sqrteps,maxdiag,mindiag,maxposdiag,u;
	 double maxoffdiag,maxinp,maxadd;
	 double maxev,minev,offrow,sdd;
	 int step,i,j,k;
	 
	 sqrteps = sqrt((double) FDVAL);
	 
	 U22 = (double*) malloc(sizeof(double) * N * N);
	
	
	for(i = 0; i < N;++i) {
		U22[i] = A[i+i*N];
	}
	
	maxdiag = array_max(U22,N);
	mindiag = array_min(U22,N);
	
	maxposdiag = 0.0;
	
	if (maxdiag > maxposdiag) {
		maxposdiag = maxdiag;
	}
	
	for(i = 0; i < N;++i) {
		U22[i] = 0.0;
	}
	u = 0.0;
	if (mindiag <= sqrteps*maxposdiag) {
		u = 2 * (maxposdiag - mindiag) * sqrteps - mindiag;
		maxdiag += u;
	}
	k = 0;
	for (i = 0; i < N;++i) {
		for(j = 0;j < N;++j) {
			if (j > i) {
				U22[k] = A[i*N+j];
				k++;
			}
		}
	}
	
	maxoffdiag = array_max_abs(U22,k);
	
	if ( maxoffdiag*(1+2*sqrteps) > maxdiag) {
		u += (maxoffdiag - maxdiag) + 2*sqrteps*maxoffdiag;
		maxdiag = maxoffdiag*(1+2*sqrteps);
	}
	
	if (maxdiag == 0) {
		u = 1;
		maxdiag = 1;
	}
	
	if (u > 0) {
		for(i=0;i < N;++i) {
			A[i*N+i] += u;
		}
	}
	
	if (maxdiag > maxoffdiag / N) {
		maxinp = sqrt(maxdiag);
	} else {
		maxinp = sqrt(maxoffdiag / N);
	}
	
	maxadd = cholmod(A,N,L,maxinp);
	
	if (maxadd > 0) {
		maxev = minev = A[0];
		for(i = 0; i < N;++i) {
			offrow = 0.0;
			step = i*N;
			
			for(j = 0; j < i;++j) {
				offrow += fabs(A[step+j]);
			}
			
			for(j = i+1; j < N;++j) {
				offrow += fabs(A[step+j]);
			}
			
			if (maxev < A[step+i] + offrow) {
				maxev = A[step+i] + offrow;
			}
			if (minev > A[step+i] - offrow) {
				minev = A[step+i] - offrow;
			}
			
		}
		sdd = (maxev - minev) * sqrteps - minev;
		if (sdd < 0) {
			sdd = 0;
		}
		if (maxadd > sdd) {
			u = sdd;
		} else {
			u = maxadd;
		}
		
		for(i = 0; i < N;++i) {
			A[i*N+i] += u; 
		}
	}
	
	maxadd = cholmod(A,N,L,0.0);
	
	free(U22);
	
	return maxadd;
	 
}

double modelhess(double *A,int N,double *dx,double *L) {
	/*
	 * Algorithm 5.5.1 Dennis, Schnabel, Numerical Methods for Unconstrained Optimization and 
	 * Non-Linear Equations
	 */ 
	 double *U22;
	 double sqrteps,maxdiag,mindiag,maxposdiag,u;
	 double maxoffdiag,maxinp,maxadd;
	 double maxev,minev,offrow,sdd;
	 int step,i,j,k;
	 
	 sqrteps = sqrt((double) FDVAL);
	 
	 U22 = (double*) malloc(sizeof(double) * N * N);
	 
	//scale
	
	for(i = 0;i < N;++i) {
		step = i*N;
		for(j = 0;j < N;++j) {
			A[step+j] /= (dx[i] * dx[j]);
		}
	}
	
	for(i = 0; i < N;++i) {
		U22[i] = A[i+i*N];
	}
	
	maxdiag = array_max(U22,N);
	mindiag = array_min(U22,N);
	
	maxposdiag = 0.0;
	
	if (maxdiag > maxposdiag) {
		maxposdiag = maxdiag;
	}
	
	for(i = 0; i < N;++i) {
		U22[i] = 0.0;
	}
	u = 0.0;
	if (mindiag <= sqrteps*maxposdiag) {
		u = 2 * (maxposdiag - mindiag) * sqrteps - mindiag;
		maxdiag += u;
	}
	k = 0;
	for (i = 0; i < N;++i) {
		for(j = 0;j < N;++j) {
			if (j > i) {
				U22[k] = A[i*N+j];
				k++;
			}
		}
	}
	
	maxoffdiag = array_max_abs(U22,k);
	
	if ( maxoffdiag*(1+2*sqrteps) > maxdiag) {
		u += (maxoffdiag - maxdiag) + 2*sqrteps*maxoffdiag;
		maxdiag = maxoffdiag*(1+2*sqrteps);
	}
	
	if (maxdiag == 0) {
		u = 1;
		maxdiag = 1;
	}
	
	if (u > 0) {
		for(i=0;i < N;++i) {
			A[i*N+i] += u;
		}
	}
	
	if (maxdiag > maxoffdiag / N) {
		maxinp = sqrt(maxdiag);
	} else {
		maxinp = sqrt(maxoffdiag / N);
	}
	
	maxadd = cholmod(A,N,L,maxinp);
	
	if (maxadd > 0) {
		maxev = minev = A[0];
		for(i = 0; i < N;++i) {
			offrow = 0.0;
			step = i*N;
			
			for(j = 0; j < i;++j) {
				offrow += fabs(A[step+j]);
			}
			
			for(j = i+1; j < N;++j) {
				offrow += fabs(A[step+j]);
			}
			
			if (maxev < A[step+i] + offrow) {
				maxev = A[step+i] + offrow;
			}
			if (minev > A[step+i] - offrow) {
				minev = A[step+i] - offrow;
			}
			
		}
		sdd = (maxev - minev) * sqrteps - minev;
		if (sdd < 0) {
			sdd = 0;
		}
		if (maxadd > sdd) {
			u = sdd;
		} else {
			u = maxadd;
		}
		
		for(i = 0; i < N;++i) {
			A[i*N+i] += u; 
		}
	}
	
	maxadd = cholmod(A,N,L,0.0);
	
	//unscale
	
	for(i = 0;i < N;++i) {
		step = i*N;
		for(j = 0;j < N;++j) {
			A[step+j] *= (dx[i] * dx[j]);
		}
	}
	
	for(i = 0;i < N;++i) {
		step = i*N;
		for(j = 0;j < i;++j) {
			L[step+j] *= dx[i];
		}
	}
	
	free(U22);
	
	return maxadd;
	 
}

void linsolve_lower(double *L,int N,double *b,double *x) {
	int i,j,c1,l;
	double *y,*A;
	double sum,temp;
	
	y = (double*) malloc(sizeof(double) *N);
	A = (double*) malloc(sizeof(double) *N *N);
	
	for(i = 0; i < N;++i) {
		y[i] = 0.;
		x[i] = 0.;
		if ( L[i*N + i] == 0.) {
			printf("The Matrix system does not have a unique solution");
			exit(1);
		}
		//printf("\n B %d",ipiv[i]);
	}
	
	// Forward Substitution
	
	y[0] = b[0]/L[0];
	for(i = 1; i < N; ++i) {
		sum = 0.;
		c1 = i*N;
		for(j = 0; j < i; ++j) {
			sum += y[j] * L[c1 + j];
		}
		y[i] = (b[i] - sum)/L[c1+i];
	}
	
	mtranspose(L,N,N,A);
	
	//Back Substitution
	
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
	free(A);
}

void newton_func(double (*funcpt)(double *,int,double *,double *),double *xi,int N,double *x) {
	int it,i;
	double *jac,*hess,*L,*xf;
	double stop,result;
	
	jac = (double*) malloc(sizeof(double) *N);
	xf = (double*) malloc(sizeof(double) *N);
	hess = (double*) malloc(sizeof(double) *N * N);
	L = (double*) malloc(sizeof(double) *N * N);
	
	it = 0;
	stop = 1.0;
	
	while (stop > 1e-05 && it < 100) {
		result = funcpt(xi,N,jac,hess);
		//mdisplay(jac,1,N);
		//mdisplay(hess,N,N);
		modelhess2(hess,N,L);
		//mdisplay(L,N,N);
		scale(jac,1,N,-1.0);
		linsolve_lower(L,N,jac,x);
		madd(x,xi,x,1,N);
		msub(x,xi,xf,1,N);
		stop = array_max_abs(xf,N);
		for(i=0;i < N;++i) {
			xi[i] = x[i];
		}
		it++;
		//printf("result %lf \n",result);
		printf("Values obtained %lf , %lf after %d Iterations \n",x[0],x[1],it);
	}
	
	free(jac);
	free(hess);
	free(xf);
	free(L);
}

void jacobian_fd(double (*funcpt)(double *,int),double *x,int N,double *f) {
	int i,j;
	double step,fd;
	double *xi;
	
	fd = sqrt((double) FDVAL);
	xi = (double*) malloc(sizeof(double) *N);
	
	for(i = 0; i < N; ++i) {
		step = fd * x[i];
		for(j = 0; j < N;++j) {
			xi[j] = x[j];
		}
		xi[i] += step;
		f[i] = (funcpt(xi,N) - funcpt(x,N))/step;
	}
	
	free(xi);
	
}

void hessian_fd(double (*funcpt)(double *,int),double *x,int N,double *f) {
	int i,j;
	double stepi,stepj,fd;
	double *xi,*xj,*xij;
	
	fd = pow((double) FDVAL,1.0/3.0);
	xi = (double*) malloc(sizeof(double) *N);
	xj = (double*) malloc(sizeof(double) *N);
	xij = (double*) malloc(sizeof(double) *N);
	
	for(i = 0; i < N;++i) {
		stepi = fd * x[i];
		for(j = 0; j < N;++j) {
			xi[j] = x[j];
			xj[j]= x[j];
			xij[j] = x[j];
		}
		xi[i] += stepi;
		xij[i] += stepi;
		for(j = 0; j < N;++j) {
			stepj = fd * x[j];
			xj[j] += stepj;
			xij[j] += stepj;
			f[i*N+j] = ((funcpt(xij,N) - funcpt(xi,N)) - (funcpt(xj,N) - funcpt(x,N)))/(stepi * stepj);
			xj[j] -= stepj;
			xij[j] -= stepj;
		}

	}
	
	free(xi);
	free(xj);
	free(xij);
	
}

void newton_fd(double (*funcpt)(double *,int),double *xi,int N,double *x) {
	int it,i;
	double *jac,*hess,*L,*xf;
	double stop,result;
	
	jac = (double*) malloc(sizeof(double) *N);
	xf = (double*) malloc(sizeof(double) *N);
	hess = (double*) malloc(sizeof(double) *N * N);
	L = (double*) malloc(sizeof(double) *N * N);
	
	it = 0;
	stop = 1.0;
	
	while (stop > 1e-05 && it < 100) {
		//result = funcpt(xi,N,jac,hess);
		jacobian_fd(funcpt,xi,N,jac);
		hessian_fd(funcpt,xi,N,hess);
		//mdisplay(jac,1,N);
		//mdisplay(hess,N,N);
		modelhess2(hess,N,L);
		//mdisplay(L,N,N);
		scale(jac,1,N,-1.0);
		linsolve_lower(L,N,jac,x);
		madd(x,xi,x,1,N);
		msub(x,xi,xf,1,N);
		stop = array_max_abs(xf,N);
		for(i=0;i < N;++i) {
			xi[i] = x[i];
		}
		it++;
		//printf("result %lf \n",result);
		printf("Values obtained %lf , %lf after %d Iterations \n",x[0],x[1],it);
	}
	
	free(jac);
	free(hess);
	free(xf);
	free(L);
}

int lnsrch(double (*funcpt)(double *,int),double *xi,double *jac,double *p,int N,double maxstep,double * dx,double stol,double *x) {
	int retval,i,j;
	double alpha,lambda,lambdamin,funcf,funci,lambdaprev,lambdatemp,funcprev;
	double lambda2,lambdaprev2,ll,den,rell;
	double *slopei,*temp1,*temp2,*ab,*rcheck;
	
	slopei = (double*) malloc(sizeof(double) *1);
	temp1 = (double*) malloc(sizeof(double) *4);
	temp2 = (double*) malloc(sizeof(double) *2);
	ab = (double*) malloc(sizeof(double) *2);
	rcheck = (double*) malloc(sizeof(double) *N);
	retval = 100;
	alpha = 1e-04;
	lambda = 1.0;
	
	mmult(jac,p,slopei,1,N,1);
	
	for(i = 0; i < N;++i) {
		if (fabs(xi[i]) > 1.0 /fabs(dx[i])) {
			den = fabs(xi[i]);
		} else {
			den = 1.0 /fabs(dx[i]);
		}
		rcheck[i] = p[i]/den;
	}
	
	rell = array_max_abs(rcheck,N);
	
	lambdamin = stol/rell;
	
	while (retval > 1) {
		scale(p,1,N,lambda);
		madd(xi,p,x,1,N);
		funcf = funcpt(x,N);
		funci = funcpt(xi,N);
		
		if (funcf <= funci + alpha *lambda *slopei[0]) {
			retval = 0;
		} else if (lambda < lambdamin) {
			retval = 1;
			for (i = 0; i < N;++i) {
				x[i] = xi[i]; // Check
			}
		} else {
			if (lambda == 1.0) {
				lambdatemp = - slopei[0] / (funcf - funci - slopei[0]); 
			} else {
				lambda2 = lambda * lambda;
				lambdaprev2 = lambdaprev * lambdaprev;
				ll = lambda - lambdaprev;
				temp1[0] = 1.0 / lambda2; temp1[1] = -1.0 /lambdaprev2;
				temp1[2] = - lambdaprev / lambda2; temp1[3] = lambda /lambdaprev2;
				temp2[0] = funcf - funci - lambda * slopei[0];
				temp2[1] = funcprev - funci - lambdaprev * slopei[0];
				mmult(temp1,temp2,ab,2,2,1);
				scale(ab,1,2,ll);
				if (ab[0] == 0.0) {
					lambdatemp = - slopei[0] / (2.0 * ab[1]);
				} else {
					lambdatemp = (-ab[1] + sqrt( ab[1] * ab[1] - 3.0 * ab[0] *slopei[0]))/ (3.0 * ab[0]);
				}
				
				if (lambdatemp > 0.5 * lambda) {
					lambdatemp = 0.5 * lambda;
				}
			}
			lambdaprev = lambda;
			funcprev = funcf;
			if (lambdatemp <= 0.1 * lambda) {
				lambda = 0.1 * lambda;
			} else {
				lambda = lambdatemp;
			}
		}
	
	}
	
	free(slopei);
	free(temp1);
	free(temp2);
	free(ab);
	free(rcheck);
	return retval;
}

int stopcheck(double fx,int N,double *xc,double *xf,double *jac,double *dx,double fsval,double gtol,double stol,int retval) {
	int rcode,i;
	double num,den;
	double stop0;
	double *scheck;
	
	rcode = 0;	
	if (retval == 1) {
		rcode = 3;
		return rcode;
	}
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

int newton_min_func(double (*funcpt)(double *,int),double *xi,int N,double *dx,double fsval,double *xf) {
	int rcode,iter;
	int i,siter,retval;
	double gtol,stol;
	double fx,num,den,stop0,maxstep,fxf;
	double *jac,*hess,*scheck,*xc,*L,*step;
	
	jac = (double*) malloc(sizeof(double) *N);
	scheck = (double*) malloc(sizeof(double) *N);
	xc = (double*) malloc(sizeof(double) *N);
	step = (double*) malloc(sizeof(double) *N);
	hess = (double*) malloc(sizeof(double) *N * N);
	L = (double*) malloc(sizeof(double) *N * N);
	
	rcode = 0;
	iter = 0;
	siter = (int) SETITER;
	gtol = pow((double) FDVAL,1.0/3.0);
	stol = gtol * gtol;
	fx = funcpt(xi,N);
	jacobian_fd(funcpt,xi,N,jac);
	
	//set values
	maxstep = 1.0;
	
	for(i = 0; i < N;++i) {
		dx[i] = 1.0 / dx[i];
	}
	
	//Check Stop0
	if (fabs(fx) > fabs(fsval)) {
			den = fabs(fx);
	} else {
			den = fabs(fsval);
	}
	for(i = 0; i < N;++i) {
		if (fabs(xi[i]) > 1.0 / fabs(dx[i])) {
			num = fabs(xi[i]);
		} else {
			num = 1.0 / fabs(dx[i]);
		}
		scheck[i] = fabs(jac[i]) * num / den;
	}
	
	stop0 = array_max_abs(scheck,N);
	
	if (stop0 <= gtol * 1e-03) {
		rcode = 1;
		for(i = 0; i < N;++i) {
			xf[i] = xi[i];
		}
		return rcode;
	}
	
	hessian_fd(funcpt,xi,N,hess);
	
	for(i = 0; i < N;++i) {
		xc[i] = xi[i];
	}
	
	while (rcode == 0 && iter < siter) {
		modelhess(hess,N,dx,L);
		scale(jac,1,N,-1.0);
		linsolve_lower(L,N,jac,step);
		scale(jac,1,N,-1.0);
		retval = lnsrch(funcpt,xc,jac,step,N,maxstep,dx,stol,xf); 
		fxf = funcpt(xf,N);
		jacobian_fd(funcpt,xf,N,jac);
		rcode = stopcheck(fxf,N,xc,xf,jac,dx,fsval,gtol,stol,retval);
		hessian_fd(funcpt,xf,N,hess);
		for(i = 0; i < N;++i) {
			xc[i] = xf[i];
		}
	}
	
	free(jac);
	free(hess);
	free(scheck);
	free(xc);
	free(L);
	free(step);
	return rcode;
}