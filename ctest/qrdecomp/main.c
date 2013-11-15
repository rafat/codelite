#include <stdio.h>
#include "qrdecomp.h"

int main(int argc, char **argv)
{
	int i,j,N,row,col,n,k,q;
	double *P,*bvec,*R,*Q;
	N = 4;
	double b;
	b = 1.;
	row = 3;
	col = 3;
	P = (double*) malloc(sizeof(double) * N * N);
	//AA = (double*) malloc(sizeof(double) * row * col);
	bvec = (double*) malloc(sizeof(double) * col);
	R = (double*) malloc(sizeof(double) * col *col);
	Q = (double*) malloc(sizeof(double) * row * col);
	double x[4] = {3,1,5,1};
	double v[4] = {0,0,0,0};
	double AA[9] = {12,-51,4,6,167,-68,-4,24,-41};
	printf("b %lf \n", b);
	b = house(x,N,v);
	printf("beta %lf \n",b);
	
	for(i =0;i < N;++i) {
		printf("v %lf \n",v[i]);
	}
	housemat(v,N,b,P);
	mdisplay(P,N,N);
	
	mmult(P,x,v,N,N,1);
	mdisplay(v,N,1);
	
	
	qr_house(AA,row,col,bvec);
	mdisplay(AA,row,col);
	getQR_house(AA,row,col,bvec,Q,R);
	mdisplay(R,col,col);
	mdisplay(Q,row,col);
	
	mmult(Q,R,AA,row,col,col);
	mdisplay(AA,row,col);
	double A[9] = {1,5,7,3,0,6,4,3,1};
	//double A[16] = {4,1,-1,2,1,4,1,-1,-1,1,4,1,2,-1,1,4};
	double H[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double eigre[3] = {0,0,0,0};
	double eigim[3] = {0,0,0,0};

	/*
	francis_iter(A,4,H);
	mdisplay(H,4,4);
	double eigre[2] = {0,0};
	double eigim[2] = {0,0};
	eig22(H,4,eigre,eigim);
	 * 
	printf("e0 %lf %lf \n",eigre[0],eigim[0]);
	printf("e1 %lf %lf \n",eigre[1],eigim[1]);
	 */
	eig(A,3,eigre,eigim);
	
	for(i=0; i < 3;++i) {
		printf("e%d %g %g \n",i,eigre[i],eigim[i]);
	}
	
	free(P);
	//free(AA);
	free(bvec);
	free(R);
	free(Q);
	
	return 0;
}
