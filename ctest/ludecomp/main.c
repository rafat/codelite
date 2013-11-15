#include <stdio.h>
#include <stdlib.h>
#include "tools/cycle.h"


int main(int argc, char **argv)
{
	int i,M,N;
	N = 4;
	M = N * N;
	double A[16] = {1,7,6,4,3,9,7,2,2,4,6,9,6,1,2,8};
	//double A[9] = {3,17,10,2,4,-2,6,18,-12};
	int p[4] = {0,0,0,0};
	double b[4] = {1.3,22.0,4.5,9.8};
	double x[4] = {0,0,0,0};
	pludecomp(A,N,p);
	linsolve(A,N,b,x,p);
	for(i = 0; i < M; ++i) {
		printf("\n %lf \n",A[i]);
	}
	for(i = 0; i < N; ++i) {
		printf("\n x [%d] %lf \n",i,x[i]);
	}
	return 0;
}
