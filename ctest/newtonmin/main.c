#include <stdio.h>
#include "newtonmin.h"
#include "conjgrad.h"


void hilb(double *A,int N)
{
	int i,j,t;
	for (i = 0; i < N; ++i) {
		t = i*N;
		for (j = 0; j < N; ++j) {
			A[t+j] = 1.0 /(i+j+1.0);
		}
	}
}

double func1(double *x,int N)
{
	double f;
	double pi,alpha,alpha2;
	pi = 3.14159;
	alpha = 1;
	alpha2 = alpha*alpha;

	//f = pow((x[0]-2.0),4.0) + pow((x[0]-2.0),2.0) * x[1]*x[1] + (x[1] + 1.0) * (x[1] + 1.0);
	f = 100 * (x[0]*x[0]*alpha2 - x[1]/alpha)* (x[0]*x[0]*alpha2 - x[1]/alpha) + (1 - x[0]*alpha) * (1 - x[0]*alpha);
	//f = (x[0] + 2 * x[1] - 7) * (x[0] + 2 * x[1] - 7) + (2*x[0] + x[1] - 5) * (2*x[0] + x[1] - 5);
	/*
	f = (x[0] + 10 * x[1] ) * (x[0] + 10 * x[1] ) + 5 * (x[2] - x[3]) * (x[2] - x[3])
	+ (x[1] - 2 * x[2]) * (x[1] - 2 * x[2]) * (x[1] - 2 * x[2]) * (x[1] - 2 * x[2])
	+ 10 * (x[0]-x[3]) * (x[0]-x[3]) * (x[0]-x[3]) * (x[0]-x[3]) ;
	 */
	/*
	f = (1.0 /(1.0 + (x[0] - x[1]) * (x[0] - x[1]))) + sin(pi *x[1] * x[2] / 2.0)
	+ exp(-(((x[0]+x[2])/x[1]) - 2) * (((x[0] +x[2])/x[1]) - 2));
	f = -1.0 * f;
	 */

	//f = (10000*x[0]*x[1] - 1) * (10000*x[0]*x[1] - 1) + (exp(-x[0]) + exp(-x[1]) - 1.0001) * (exp(-x[0]) + exp(-x[1]) - 1.0001);
	//f = log(f);

	//f = (x[0] - 1e06)*(x[0] - 1e06) + (x[1] - 2*1e-06)*(x[1] - 2*1e-06) +(x[0]*x[1] - 2)*(x[0]*x[1] - 2);
	//f = 100 * (x[1]-x[0]*x[0]) * (x[1]-x[0]*x[0]) + ( 1.0 - x[0] ) * ( 1.0 - x[0] )  + 90 *(x[3]-x[2]*x[2])*(x[3]-x[2]*x[2])
	//+ ( 1.0 - x[2])*( 1.0 - x[2]) + 10 * (x[1] + x[3] - 2)*(x[1] + x[3] - 2) + 0.1 * (x[1] - x[3]) *(x[1] - x[3]);
	return f;
}

int main(int argc, char **argv)
{

	int N,i,rcode,method;
	double fsval,alpha,delta;
	N = 2;
	method = 1;
	alpha = 1;delta = -1.0;
	double xi[2] = {6.39,-0.221};
	double x[2] = {0,0};
	//double dx[2] = {1.0e5,1.0e-05};
	double dx[2] = {1.0/alpha,1.0*alpha};
	fsval = 1.0;
	//double f1[2] = {0,0};
	//double f2[4] = {0,0,0,0};
	//newton_fd(func4,xi,N,x);
	//rcode = newton_min_func(func1,xi,N,dx,fsval,x);
	//rcode = newton_min_trust(func1,xi,N,dx,fsval,delta,method,x);
	rcode = conjgrad_min_lin(func1,xi,N,x);
	mdisplay(x,1,N);
	printf("Termination Code %d %g \n",rcode,func1(x,N));
	

	return 0;
}
