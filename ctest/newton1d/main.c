#include <stdio.h>
#include "newton1d.h"

void newton_1d() {
	double f,fd,x,hc,xc,xp,fp;
	int it;
	
	x = 7;
	it = 0;
	//f = sin(x) - cos(2 * x);
	f = x*x*x - 7*x*x + 11*x - 5;
	fp = f;
	
	//Method 1 : Newton's Method with explicit derivatives
	
	while (fabs(f) >= 1e-06 && it < 100) {
		//f = sin(x) - cos(2 * x);
		f = x*x*x - 7*x*x + 11*x - 5;
		while (fabs(f) > fabs(fp)) {
			x = (x + xp) / 2;
			fp = f;
			xp = x;
			//f = sin(x) - cos(2 * x);
			f = x*x*x - 7*x*x + 11*x - 5;
		}
		//fd = cos(x) + 2 * sin(2*x);
		fd = 3*x*x - 14*x + 11;
		it++;
		fp = f;
		xp = x;
		x = x - f/fd;
		//printf("Value f %g after %d iterations \n",f,it);
		
	}
	
	printf("Method 1 : Value obtained %lf after %d Iterations \n",x,it);
	
	x = 7;
	it = 0;
	//f = sin(x) - cos(2 * x);
	f = x*x*x - 7*x*x + 11*x - 5;
	fp = f;
	
	//Method 2 : Newton's Method with Finite differences
	
	while (fabs(f) >= 1e-06 && it < 100) {
		//f = sin(x) - cos(2 * x);
		f = x*x*x - 7*x*x + 11*x - 5;
		while (fabs(f) > fabs(fp)) {
			x = (x + xp) / 2;
			fp = f;
			xp = x;
			//f = sin(x) - cos(2 * x);
			f = x*x*x - 7*x*x + 11*x - 5;
		}
		//fd = cos(x) + 2 * sin(2*x);
		hc = FDVAL * x;
		xc = x + hc;
		fd = (xc*xc*xc - 7*xc*xc + 11*xc - 5 - f )/hc;
		it++;
		fp = f;
		xp = x;
		x = x - f/fd;
		//printf("Value f %g after %d iterations \n",f,it);
		
	}
	
	printf("Method 2 : Value obtained %lf after %d Iterations \n",x,it);
	
	x = 7;
	it = 0;
	//f = sin(x) - cos(2 * x);
	f = x*x*x - 7*x*x + 11*x - 5;
	fp = f;
	
	//Method 3 : Newton's Secant Method
	
	while (fabs(f) >= 1e-06 && it < 100) {
		//f = sin(x) - cos(2 * x);
		f = x*x*x - 7*x*x + 11*x - 5;
		while (fabs(f) > fabs(fp)) {
			x = (x + xp) / 2;
			fp = f;
			xp = x;
			//f = sin(x) - cos(2 * x);
			f = x*x*x - 7*x*x + 11*x - 5;
		}
		//fd = cos(x) + 2 * sin(2*x);
		if (it == 0) {
			hc = FDVAL * x;
			xc = x + hc;
			fd = (xc*xc*xc - 7*xc*xc + 11*xc - 5 - f )/hc;
		} else {
			fd = (f - fp)/ (x - xp);
		}
		it++;
		xp = x;
		fp = f;
		x = x - f/fd;
		//printf("Value f %g after %d iterations \n",f,it);
		
	}
	
	printf("Method 3 : Value obtained %lf after %d Iterations \n",x,it);
}

void newton_nle_ad() {
	int N,it;
	double *x,*xi,*f,*J;
	double stop;
	int *p;
	
	N = 2;
	it = 0;
	x = (double*) malloc(sizeof(double) * N);
	xi = (double*) malloc(sizeof(double) * N);
	f = (double*) malloc(sizeof(double) * N);
	J = (double*) malloc(sizeof(double) * N * N);
	p = (int*) malloc(sizeof(int) * N);
	
	x[0] = 2;x[1] = 3;
	f[0] = x[0]*x[0] + x[1]*x[1] - 2.0;f[1] = pow((double) EPSILON,x[0] - 1.0) + x[1]*x[1]*x[1] - 2.0;
	
	stop = array_max_abs(f,N);
	
	while (stop >= 1e-06 && it < 100) {
		f[0] = x[0]*x[0] + x[1]*x[1] - 2.0;f[1] = pow((double) EPSILON,x[0] - 1.0) + x[1]*x[1]*x[1] - 2.0;
		stop = array_max_abs(f,N);
		J[0] = 2 * x[0]; J[1] = 2 * x[1];
		J[2] = pow((double) EPSILON,x[0] - 1.0);J[3] = 3*x[1]*x[1];
		ludecomp(J,N,p);
		linsolve(J,N,f,p,xi);
		msub(x,xi,x,N,N);
		it++;
		printf("Method 1 : Values obtained %lf , %lf after %d Iterations \n",x[0],x[1],it);
	
	}
	
	printf("Method 1 : Values obtained %lf , %lf after %d Iterations \n",x[0],x[1],it);
	
	free(x);
	free(xi);
	free(f);
	free(J);
	free(p);
}

void newton_nle_fd() {
	int N,it;
	double *x,*y,*xi,*f,*J;
	double stop,h;
	int *p;
	
	N = 2;
	it = 0;
	x = (double*) malloc(sizeof(double) * N);
	y = (double*) malloc(sizeof(double) * N);
	xi = (double*) malloc(sizeof(double) * N);
	f = (double*) malloc(sizeof(double) * N);
	J = (double*) malloc(sizeof(double) * N * N);
	p = (int*) malloc(sizeof(int) * N);
	
	x[0] = 2;x[1] = 3;
	f[0] = func1(x,N);f[1] = func2(x,N);
	
	stop = array_max_abs(f,N);
	
	while (stop >= 1e-06 && it < 100) {
		f[0] = func1(x,N);f[1] = func2(x,N);
		stop = array_max_abs(f,N);
		//J[0] = 2 * x[0]; J[1] = 2 * x[1];
		//J[2] = pow((double) EPSILON,x[0] - 1.0);J[3] = 3*x[1]*x[1];
		h = FDVAL * x[0];
		y[0] = x[0] + h;
		y[1] = x[1];
		J[0] = (func1(y,N) - func1(x,N) )/ h;
		J[2] = (func2(y,N) - func2(x,N) )/ h;
		
		h = FDVAL * x[1];
		y[0] = x[0];
		y[1] = x[1] + h;
		J[1] = (func1(y,N) - func1(x,N) )/ h;
		J[3] = (func2(y,N) - func2(x,N) )/ h;
		
		ludecomp(J,N,p);
		linsolve(J,N,f,p,xi);
		msub(x,xi,x,N,N);
		it++;
		printf("Method 2 : Values obtained %lf , %lf after %d Iterations \n",x[0],x[1],it);
	
	}
	
	printf("Method 2 : Values obtained %lf , %lf after %d Iterations \n",x[0],x[1],it);
	
	free(x);
	free(y);
	free(xi);
	free(f);
	free(J);
	free(p);
}

double meps() {
	double eps = 1.0;
	while ( (1.0 + eps) != 1.0) {
		eps /= 2.0;
	}
	
	eps *= 2.0;
	return eps;
}

void hilb(double *A,int N) {
	int i,j,t;
	for (i = 0; i < N;++i) {
		t = i*N;
		for (j = 0; j < N;++j) {
			A[t+j] = 1.0 /(i+j+1.0);
		}
	}
}

int main(int argc, char **argv)
{
	/*
	int N;
	double result,eps;
	double t[2] = {1.2,2.3};
	N = 2;
	newton_nle_ad();
	newton_nle_fd();
	//funcpt = &func1;
	result = solvefunc(func2,t,N);
	printf("\n val %lf \n",result);
	eps = meps();
	printf("\n val %g \n",eps);
	
	return 0;
	 */
/*
	int N,i,M,success;
	double *B,*C,*L;
	double madd,minp;
	N = 3;
	A = (double*) malloc(sizeof(double) * N * N);
	//double A[9] = {10,20,30,20,45,80,30,80,171};
	B = (double*) malloc(sizeof(double) * N * N);
	C = (double*) malloc(sizeof(double) * N * N);
	L = (double*) malloc(sizeof(double) * N * N);
	minp = 0.0;
	hilb(A,N);
	mdisplay(A,N,N);
	success = chol(A,N);
	if (success == 0) { 
		mdisplay(A,N,N);
		mtranspose(A,N,N,B);
		mmult(B,A,C,N,N,N);
		mdisplay(C,N,N);
	} else {
		printf("The Input Matrix Is Not Positive Definite");
	}
	
	madd = modelhess(C,N,L);
	mdisplay(L,N,N);
	printf("mdouble %lf \n",madd);
	free(A);
	free(B);
	free(C);
	free(L);
	 */
	int N,i,rcode;
	double fsval;
	N = 2;
	double xi[2] = {11,21}; 
	double x[2] = {0,0}; 
	double dx[2] = {1.0,1.0};
	fsval = 1.0;
	//double f1[2] = {0,0};
	//double f2[4] = {0,0,0,0};
	//newton_fd(func4,xi,N,x);
	rcode = newton_min_func(func4,xi,N,dx,fsval,x);
	mdisplay(x,1,N);
	printf("Termination Code %d \n",rcode);
	//jacobian_fd(func4,xi,N,f1);
	//mdisplay(f1,1,N);
	//hessian_fd(func4,x,N,f2);
	//mdisplay(f2,N,N);

	return 0;
}
