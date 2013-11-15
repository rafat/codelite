#include <stdio.h>
#include "matrix.h"
#include "matmul.h"
#include "tools/cycle.h"

int main(int argc, char **argv)
{
	int i,M,N;
	N = 3;
	M = N * N;
	//double A[16] = {1,7,6,4,3,9,7,2,2,4,6,9,6,1,2,8};
	double inv[9] = {1,7,6,4,3,9,7,2,2};
	double A[9] = {3,17,10,2,4,-2,6,18,-12};
	int p[3] = {0,0,0};
	double b[4] = {0.0,1.0,0.0,0.0};
	double x[3] = {0,0,0};
	ludecomp(A,N,N,p);
	minverse(A,N,p,inv);
	//linsolve(A,N,b,x,p);
	for(i = 0; i < M; ++i) {
		printf("\n %lf \n",inv[i]);
	}
	/*
	int ra,rb,ca,cb,N,NX,NY,i;
	ticks t0,t1,t2,t3;
	double tn,tr;
	int r2,c2;
	ra = 570;
	ca = 621;
	cb = 355;
	
	rb = ca;
	NX = ra*rb;
	NY = rb*cb;
	double *X,*Y,*Z1,*Z2,*Zd;
	X = (double*) malloc(sizeof(double) * NX);
	Y = (double*) malloc(sizeof(double) * NY);
	
	for (i = 0; i < NX; i++) {
		X[i] = (double) i/ra;
	}
	
	for (i = 0; i < NY; i++) {
		Y[i] = (double) (i+3.0)/rb;
	}
	//mdisplay(X,4,2);
	N = ra * cb;
	Z1 = (double*) malloc(sizeof(double) * N);
	Z2 = (double*) malloc(sizeof(double) * N);
	Zd = (double*) malloc(sizeof(double) * N);
	t0 = getticks();
	//tmult(X,Y,Z1,ra,ca,cb);
	//matmul_matmul(ra, ca, cb, X, ca, Y, cb, Z1, cb);
	rmult(X,Y,Z1,ra,ca,cb);
	t1 = getticks();
	t2 = getticks();
	smult(X,Y,Z2,ra,ca,cb);
	//loop_multiply (ra, Z2, Y,X);
	t3 = getticks();
	//mdisplay(Z,ra,cb);
	tn = elapsed(t1,t0);
	tr = elapsed(t3,t2);
	for (i = 0; i < N; i++) {
		Zd[i] = Z1[i] - Z2[i];
		//printf("\n %lf %lf \n",Z1[i],Z2[i]);
	}
	
	printf("\n Absolute Max Value %g \n",array_max_abs(Zd,N));
	printf("\n TN = %lf , TR = %lf \n",tn,tr);
	//printf("\n Absolute Max Value %g \n",array_max_abs(Z2,N));
	free(X);
	free(Y);
	free(Z1);
	free(Z2);
	free(Zd);
  
/*	
	//i = findrec(&ra,&ca,&cb);
	//printf("\n %d %d %d %d \n",ra,ca,cb,i);
	double *X,*Y,*Z;
	//printf("Pointer %x\n", *Y);
	int r,c;
	X = (double*) malloc(sizeof(double) * ra * ca);
	Z = (double*) malloc(sizeof(double) * ra * ca);
	printf("Pointer %x\n", *X);
	
	for (i = 0; i < ra * ca; i++) {
		X[i] = (double) i / ra;
	}
	printf("Pointer %x\n", *X);
	mdisplay(X,ra,ca);
	r = 3;
	c = 4;
	r2 = ra +r;
	c2 = ca + c;
	Y = (double*) malloc(sizeof(double) * r2 * c2);
	//printf("Pointer %x\n", *Y);
	add_zero_pad(X,ra,ca,r,c,Y);
	mdisplay(Y,ra+r,ca+c);
	
	remove_zero_pad(Y,ra+r,ca+c,r,c,Z);
	
	mdisplay(Z,ra,ca);
	
	madd_stride(Z,Z+(ca/2),Z+(ra/2)*ca,ra/2,ca/2,ca,ca,ca);
	
	mdisplay(Z,ra,ca);
	free(X);
	free(Z);
	//printf("Pointer %x\n", *Y);
	free(Y);
	 */ 
	return 0;
}
