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

static int zoom(double (*funcpt)(double *,int),double *xi,double *jac,double *p,int N,double *dx,double maxstep,
		double lambdalo,double lambdahi,double alpha,double beta,double lambdamin,double *xf) {
	int retval,iter,i;
	double lambda;
	double *slopei,*slopen,*pl,*jacf,*xlo;
	double funci,funcf,funcprev,funclo;
	
	slopei = (double*) malloc(sizeof(double) *1);
	slopen = (double*) malloc(sizeof(double) *1);
	pl = (double*) malloc(sizeof(double) *N);
	jacf = (double*) malloc(sizeof(double) *N);
	xlo = (double*) malloc(sizeof(double) *N);
	
	retval = 5;
	
	mmult(jac,p,slopei,1,N,1);
	
	funci = funcpt(xi,N);
	
	if (funci >= DBL_MAX || funci <= -DBL_MAX) {
		printf("Program Exiting as the function value exceeds the maximum double value");
		exit(1);
	}
	funcprev = funci;
	iter = 0;
	
	
	while (retval > 1 && iter < (int) SETITER) {
		iter++;
		lambda = (lambdalo + lambdahi) / 2.0; // replace bisection
		if (lambda < lambdamin) {
			retval = 1;
			for (i = 0; i < N;++i) {
				xf[i] = xi[i]; // Check
			}
		}
		for(i = 0; i < N;++i) {
			pl[i] = p[i] * lambdalo;
		}
		madd(xi,pl,xlo,1,N);
		funclo = funcpt(xlo,N);
		if (funclo >= DBL_MAX || funclo <= -DBL_MAX) {
			printf("Program Exiting as the function value exceeds the maximum double value");
			exit(1);
		}
		for(i = 0; i < N;++i) {
			pl[i] = p[i] * lambda;
		}
		madd(xi,pl,xf,1,N);
		funcf = funcpt(xf,N);
		if (funcf >= DBL_MAX || funcf <= -DBL_MAX) {
			printf("Program Exiting as the function value exceeds the maximum double value");
			exit(1);
		}
		grad_fd(funcpt,xf,N,dx,jacf);
		mmult(jacf,p,slopen,1,N,1);
		if (funcf > funci + alpha *lambda *slopei[0] || funcf >= funclo) {
			lambdahi = lambda;
		} else {
			//printf("YES %d %g %g \n",iter,fabs(slopen[0]),-beta*slopei[0]);
			if (fabs(slopen[0]) <= -beta*slopei[0]) {
				retval = 0;
			}
			if (slopen[0] * (lambdahi - lambdalo) >= 0) {
				lambdahi = lambdalo;
			}
			lambdalo = lambda;
		}
	} 
	
	free(slopei);
	free(slopen);
	free(pl);
	free(jacf);
	free(xlo);
	
	return retval;
}

int swolfe(double (*funcpt)(double *,int),double *xi,double *jac,double *p,int N,double *dx,double maxstep,double stol,double *xf) {
	double lambda,alpha,beta,nlen,lambdap,rho;
	int i,retval,iter;
	double *slopei,*slopen,*pl,*jacf,*temp1,*temp2,*ab,*rcheck;
	double funci,funcf,funcprev,lambdatemp;
	double lambda2,lambdaprev2,ll,lambdamin,rell,den;
	
	slopei = (double*) malloc(sizeof(double) *1);
	slopen = (double*) malloc(sizeof(double) *1);
	pl = (double*) malloc(sizeof(double) *N);
	jacf = (double*) malloc(sizeof(double) *N);
	temp1 = (double*) malloc(sizeof(double) *4);
	temp2 = (double*) malloc(sizeof(double) *2);
	ab = (double*) malloc(sizeof(double) *2);
	rcheck = (double*) malloc(sizeof(double) *N);
	
	retval = 100;
	lambdap = 0.0;
	alpha = 1.0e-04;
	beta = 0.40;
	nlen = 0.0;
	iter = 0;
	rho = 0.8;
	
	
	for(i = 0; i < N;++i) {
		nlen += dx[i] * p[i] * dx[i] * p[i];
	}
	nlen = sqrt(nlen);
	
	if (nlen > maxstep) {
		scale(p,1,N,maxstep/nlen);
		nlen = maxstep;
	}
	mmult(jac,p,slopei,1,N,1);
	
	funci = funcpt(xi,N);
	lambda = 1.0;
	printf("%d %g \n",iter,funci);
	if (funci >= DBL_MAX || funci <= -DBL_MAX) {
		printf("Program Exiting as the function value exceeds the maximum double value");
		exit(1);
	}
	funcprev = funci;
	
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
	
	while (retval > 1 && iter < (int) SETITER) {
		iter++;
		for(i = 0; i < N;++i) {
			pl[i] = p[i] * lambda;
		}
		
		madd(xi,pl,xf,1,N);
		funcf = funcpt(xf,N);
		printf("%d %g %g \n",iter,funci,funcf);
		//mdisplay(xf,1,N);
		if (funcf >= DBL_MAX || funcf <= -DBL_MAX) {
			printf("Program Exiting as the function value exceeds the maximum double value");
			exit(1);
		}
		if (funcf > funci + alpha *lambda *slopei[0] || (iter > 1 && funcf > funcprev)) {
			//call zoom
			retval = zoom(funcpt,xi,jac,p,N,dx,maxstep,lambdap,lambda,alpha,beta,lambdamin,xf);
		}
		grad_fd(funcpt,xf,N,dx,jacf);
		mmult(jacf,p,slopen,1,N,1);
		if (fabs(slopen[0]) <= -beta * slopei[0]) {
			retval = 0;
		}
		
		if (slopen[0] >= 0) {
			// call zoom
			retval = zoom(funcpt,xi,jac,p,N,dx,maxstep,lambda,lambdap,alpha,beta,lambdamin,xf);
			//printf("YES");
		}
		//lambdap = lambda;
		//lambda = rho * lambda + nlen * (1-rho); // replace with quadratic and cubic interpolation
		//funcprev = funcf;
		
		if (lambda == 1.0) {
			lambdatemp = - slopei[0] / (2.0 * (funcf - funci - slopei[0])); 
		} else if (lambda < lambdamin) {
			if (retval != 5) {
				retval = 1;
			}
			for (i = 0; i < N;++i) {
				xf[i] = xi[i]; // Check
			}
		} else {
			lambda2 = lambda * lambda;
			lambdaprev2 = lambdap * lambdap;
			ll = lambda - lambdap;
			temp1[0] = 1.0 / lambda2; temp1[1] = -1.0 /lambdaprev2;
			temp1[2] = - lambdap / lambda2; temp1[3] = lambda /lambdaprev2;
			temp2[0] = funcf - funci - lambda * slopei[0];
			temp2[1] = funcprev - funci - lambdap * slopei[0];
			mmult(temp1,temp2,ab,2,2,1);
			scale(ab,1,2,1.0/ll);
			//printf("ab %g \n", ab[0]);
			if (ab[0] == 0.0) {
				lambdatemp = - slopei[0] / (2.0 * ab[1]);
			} else if (ab[1] * ab[1] - 3.0 * ab[0] *slopei[0] < 0.0) {
				lambdatemp = 0.5 * lambda;
			} else {
				lambdatemp = (-ab[1] + sqrt( ab[1] * ab[1] - 3.0 * ab[0] *slopei[0]))/ (3.0 * ab[0]);
			}
			if (lambdatemp > 0.5 * lambda) {
				lambdatemp = 0.5 * lambda;
			}
		
			
		}
		lambdap = lambda;
		funcprev = funcf;
		if (lambdatemp <= 0.1 * lambda) {
			lambda = 0.1 * lambda;
		} else {
			lambda = lambdatemp;
		}
		
	}
	
	free(slopei);
	free(slopen);
	free(pl);
	free(jacf);
	free(temp1);
	free(temp2);
	free(ab);
	free(rcheck);
	
	return retval;
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

int cgp(double (*funcpt)(double *,int),double *xi,int N,double *dx,double *xf) {
	int i,rcode,retval,k,restart;
	int iter,siter;
	double *temp,*rk,*pk,*jac,*jacf,*apk;
	double fsval,stol,gtol,fxf;
	double maxstep,dt1,dt2;
	
	temp = (double*) malloc(sizeof(double) *8);
	rk = (double*) malloc(sizeof(double) *N);
	pk = (double*) malloc(sizeof(double) *N);
	apk = (double*) malloc(sizeof(double) *N);
	jac = (double*) malloc(sizeof(double) *N);
	jacf = (double*) malloc(sizeof(double) *N);
	
	iter = 0;
	k = 0;
	rcode = 0;
	restart = N;
	siter = (int) SETITER;
	gtol = pow((double) FDVAL,1.0/3.0);
	stol = gtol * gtol;
	
	// Values that may not be needed
	
	fsval = 1.0;
	
	maxstep = 1000.0; // Needs to be set at a much higher value proportional to l2 norm of dx
	dt1 = dt2 = 0.0;
	for(i = 0; i < N;++i) {
		dt1 += dx[i] * dx[i];
		dt2 += dx[i] * xi[i] * dx[i] * xi[i];
	}

	dt1 = sqrt(dt1);
	dt2 = sqrt(dt2);
	
	if (dt1 > dt2) {
		maxstep *= dt1;
	} else {
		maxstep *= dt2;
	}
	
	
	grad_fd(funcpt,xi,N,dx,jac);
	for(i = 0; i < N;++i) {
		pk[i] = -jac[i];
		xf[i] = xi[i];
	}
		
	if (restart < N) {
		restart = N;
	}	
	
	while (rcode == 0 && iter < siter) {
		iter++;
		k++;
		mmult(jac,jac,temp,1,N,1);
		//retval = lnsrchcg(funcpt,xi,jac,pk,N,dx,maxstep,stol,xf,jacf); 
		retval = swolfe(funcpt,xi,jac,pk,N,dx,maxstep,stol,xf);
		if (retval == 100) {
			printf("The Linesearch Algorithm didn't converge");
			break;
		}
		//mdisplay(xf,1,N);
		grad_fd(funcpt,xf,N,dx,jacf);
		mmult(jacf,jacf,temp+1,1,N,1);
		temp[2] = temp[1]/temp[0]; // beta
		if (k == restart) {
			for(i = 0; i < N;++i) {
				pk[i] = - jacf[i]; 
			}
			k = 0;
		} else {
			for(i = 0; i < N;++i) {
				pk[i] = temp[2] * pk[i] - jacf[i]; 
			}
		}
		
		fxf = funcpt(xi,N);
		printf("ITER %d %g \n",k,fxf);
		if (fxf >= DBL_MAX || fxf <= -DBL_MAX) {
			printf("Program Exiting as the function value exceeds the maximum double value");
			exit(1);
		}
		
		rcode = stopcheck(fxf,N,xi,xf,jac,dx,fsval,gtol,stol,retval);
		for(i = 0; i < N;++i) {
			xi[i] = xf[i];
			jac[i] = jacf[i];
		}
	}
	
	free(temp);
	free(rk);
	free(pk);
	free(jac);
	free(apk);
	free(jacf);
	
	return rcode;
}

int conjgrad_min_lin(double (*funcpt)(double *,int),double *xi,int N,double *dx,double *xf) {
	int rcode,i;
	double *A,*b;
	
	/*
	 * Return Codes
	 * 
	 * Codes 1,2,3 denote possible success.
	 * Codes 0 and 4 denote failure.
	 * 
	 * 1 - df(x)/dx <= gtol achieved so xf may be the local minima.
	 * 2 - Distance between the last two steps is less than stol or |xf - xi| <= stol so xf may be the local minima.
	 * 3 - Global Step failed to locate a lower point than xi so xi may be the local minima.
	 * 4 - Iteration Limit exceeded. Convergence not achieved.
	 * 5.- Only First Strong Wolfe Condition Achieved (Using Conjugate Gradient Method)
	 * 
	 */ 
	
	b = (double*) malloc(sizeof(double) *N);
	A = (double*) malloc(sizeof(double) *N * N);
	
	for(i = 0; i < N;++i) {
		xi[i] *= dx[i];
		dx[i] = 1.0 / dx[i];
	}
	
	
	
	//hessian_fd(funcpt,xi,N,dx,A);
	//grad_fd(funcpt,xi,N,dx,b);
	//A[0] = 4; A[1] = 1; A[2] = 1; A[3] = 3;
	//b[0] = 1; b[1] = 2;
	//xi[0] = 2; xi[1] = 1;
	
	//rcode = cgpc(xi,N,A,b,xf);
	rcode = cgp(funcpt,xi,N,dx,xf);
	free(A);
	free(b);
	
	return rcode;
}