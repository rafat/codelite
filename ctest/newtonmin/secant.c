#include "secant.h"

static void jrotate(double *A,int N,double a,double b,int i) {
	int j,r,r1;
	double c,s,den,y,w;
	r = i * N;
	r1 = r + N;
	if (a == 0.0) {
		c = 0.0;
		s = signx(b);
	} else {
		den = sqrt(a*a+b*b);
		c = a / den;
		s = b / den;
	}
	
	for (j = i; j < N;++j) {
		y = A[r+j];
		w = A[r1+j];
		A[r+j] = c*y - s*w;
		A[r1+j] = s*y + c*w;
	}
}

void qrupdate(double *R,int N,double *u,double *v) {
	int i,j,k;
	
	k = N-1;
	while (u[k] == 0 && k > 0) {
		k--;
	}
	
	for (i = k-1; i >= 0;--i) {
		jrotate(R,N,u[i],-u[i+1],i);
		if (u[i] == 0) {
			u[i] = fabs(u[i+1]);
		} else {
			u[i] = sqrt(u[i]*u[i] + u[i+1]*u[i+1]);
		}
	}
	
	for(i = 0; i < N;++i) {
		R[i] += u[0] * v[i];
	}
	
	for(i = 0; i <= k-1;++i) {
		j = i * N;
		k = j + N;
		jrotate(R,N,R[j+i],-R[k+i],i);
	}
}

void bfgs_naive(double *H,int N,double *xi,double *xf,double *jac,double *jacf) {
	int i,j,supd,ct;
	double sn,yn,fd,yt,jacm;
	double *sk,*yk,*temp,*temp2,*t;
	
	sk = (double*) malloc(sizeof(double) *N);
	yk = (double*) malloc(sizeof(double) *N);
	temp = (double*) malloc(sizeof(double) *1);
	temp2 = (double*) malloc(sizeof(double) *1);
	t = (double*) malloc(sizeof(double) *N);
	
	msub(xf,xi,sk,1,N);
	msub(jacf,jac,yk,1,N);
	
	sn = l2norm(sk,N);
	yn = l2norm(yk,N);
	fd = sqrt((double) FDVAL);
	
	mmult(yk,sk,temp,1,N,1);
	
	if (temp[0] >= fd*sn*yn) {
		supd = 1;
		for(i = 0; i < N;++i) {
			t[i] = 0.0;
			for(j = 0; j < i;++j) {
				ct = j * N;
				t[i] += H[ct+i] * sk[j]
			}
			
			for(j = i; j < N;++j) {
				ct = i * N;
				t[i] += H[ct + j] * sk[j];
			}
			yt = fabs(yk[i] - t[i]);
			if (jac > jacf) {
				jacm = jac;
			} else {
				jacm = jacf;
			}
			if (yt >= fd * jacm) {
				supd = 0;
			}
		}
		if (supd == 0) {
			mmult(sk,t,temp2,1,N,1);
			for(i = 0;i < N; ++i) {
				ct = i * N;
				for(j = 0; j < N;++j) {
					H[ct+j] += (yk[i]*yk[j]/temp[0]);
					H[ct+j] -= (t[i]*t[j]/temp[2]);
				}
			}
		} 
	} 
	
	free(sk);
	free(yk);
	free(temp);
	free(temp2);
	free(t);
}

void inithess_naive(double *H,int N,double fi,double fsval,double *dx) {
	int i,j,ct;
	double temp;
	
	if (fabs(fi) > fsval) {
		temp = fabs(fi);
	} else {
		temp = fsval;
	}
	
	for(i = 0; i < N;++i) {
		ct = i *N;
		for(j = 0; j < N;++j) {
			if (i == j) {
				H[ct+j] = temp * dx[i] * dx[i];
			} else {
				H[ct+j] = 0.0;
			}
			
		}
	}
	
}

int bfgs_min_naive(double (*funcpt)(double *,int),double *xi,int N,double *dx,double fsval,double *xf)  {
	int rcode,iter;
	int i,siter,retval;
	double gtol,stol,dt1,dt2;
	double fx,num,den,stop0,maxstep,fxf;
	double *jac,*hess,*scheck,*xc,*L,*step;
	
	jac = (double*) malloc(sizeof(double) *N);
	scheck = (double*) malloc(sizeof(double) *N);
	xc = (double*) malloc(sizeof(double) *N);
	step = (double*) malloc(sizeof(double) *N);
	hess = (double*) malloc(sizeof(double) *N * N);
	L = (double*) malloc(sizeof(double) *N * N);
	
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
	 * 
	 */ 
	
	rcode = 0;
	iter = 0;
	siter = (int) SETITER;
	gtol = pow((double) FDVAL,1.0/3.0);
	
	stol = gtol * gtol;
	//set values
	for(i = 0; i < N;++i) {
		xi[i] *= dx[i];
		dx[i] = 1.0 / dx[i];
	}
	fx = funcpt(xi,N);
	if (fx >= DBL_MAX || fx <= -DBL_MAX) {
		printf("Program Exiting as the function value exceeds the maximum double value");
		exit(1);
	}
	
	grad_fd(funcpt,xi,N,dx,jac);
	
	
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
	
	printf("dt1 dt2 %g \n", maxstep);

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
	
	hessian_fd(funcpt,xi,N,dx,hess);
	
	for(i = 0; i < N;++i) {
		xc[i] = xi[i];
	}
	
	while (rcode == 0 && iter < siter) {
		iter++;
		
		modelhess(hess,N,dx,L);
		scale(jac,1,N,-1.0);
		//mdisplay(hess,N,N);
		
		linsolve_lower(L,N,jac,step);
		
		scale(jac,1,N,-1.0);
		retval = lnsrch(funcpt,xc,jac,step,N,dx,maxstep,stol,xf); 
		
		fxf = funcpt(xf,N);
		if (fxf >= DBL_MAX || fxf <= -DBL_MAX) {
			printf("Program Exiting as the function value exceeds the maximum double value");
			exit(1);
		}
		printf("%d \n",iter);
		grad_fd(funcpt,xf,N,dx,jac);
		rcode = stopcheck(fxf,N,xc,xf,jac,dx,fsval,gtol,stol,retval);
		hessian_fd(funcpt,xf,N,dx,hess);
		for(i = 0; i < N;++i) {
			xc[i] = xf[i];
		}
	}
	
	if (rcode == 0 && iter >= siter) {
		rcode = 4;
	}
	/*
	for(i = 0; i < N;++i) {
		xf[i] *= dx[i];
		dx[i] = 1.0 / dx[i];
	}
	*/
	free(jac);
	free(hess);
	free(scheck);
	free(xc);
	free(L);
	free(step);
	return rcode;
}