/*
 * matrix.c
 *
 *  Created on: Jul 1, 2013
 *      Author: Rafat Hussain
 */

#include "matrix.h"

double array_max_abs(double *array,int N) {
	int i;
	double m = 0.0;
	for (i = 0; i < N;++i) {
		if (fabs(array[i]) > m ) {
			m = fabs(array[i]);
		}
	}
	return m;
}

void dtranspose(double *sig, int rows, int cols,double *col) {
    int max,ud,i,k;
    if (rows >= cols) {
    	max = cols;
    } else {
    	max = rows;
    }
    ud = 0;
	for (i= -rows + 1; i < cols; i++) {
		if (i <= 0) {
			ud++;
			if (ud >= max)
				ud = max;
			for (k = 0; k < ud; k++) {
				col[k*rows+k-i] = sig[(k-i)*cols+k];
			}
		} else {
			if (i - cols + rows > 0) {
				ud--;
				if (ud >= max)
					ud = max;
			}
			for (k = 0; k < ud; k++) {
				col[(k+i)*rows+k] = sig[k*cols+k+i];
			}
		}

	}

}

void stranspose(double *sig, int rows, int cols,double *col) {
	int i,j,t,u;
	for (i=0; i < rows; i++) {
		t = i * cols;
		u = 0;
		for (j=0; j < cols; j++) {
			col[u+i] = sig[j+t];
			u+=rows;
		}
	}
	
}

void rtranspose(double *m, int rows, int cols,double *n, int r, int c) {
	register int i,j;
	int rm,cm;
	int rm1,cm1,rm2,cm2;

	if (rows <= 64 && cols <= 64) {
		for (i = 0; i < rows; ++i) {
			for (j = 0; j < cols; ++j) {
				n[i+j*r] = m[j+i*c];
				//cout << *(n+i+j*r) << " ";
			}
		}
		//cout << endl;

	} else if (cols >= rows) {
		rm = rows;
		cm1 = (int) ceil((double) cols/2.0);
		cm2 = cols - cm1;

		rtranspose(m,rm,cm1,n,r,c);
		rtranspose(m+cm1,rm,cm2,n+cm1*r,r,c);
	} else if (rows > cols) {
		rm1 = (int) ceil((double) rows/2.0);
		rm2 = rows - rm1;
		cm = cols;
		rtranspose(m,rm1,cm,n,r,c);
		rtranspose(m+rm1*c,rm2,cm,n+rm1,r,c);
	}

}

void ctranspose(double *sig, int rows, int cols,double *col) {
	int r,c;
	r= rows;
	c = cols;
	if (rows >= 64 || cols >= 64)  {
		rtranspose(sig,rows,cols,col,r,c);
	} else {
		stranspose(sig,rows,cols,col);
	}
}

void mtranspose(double *sig, int rows, int cols,double *col) {
	
	if (rows >= 512 && cols >= 512)  {
		ctranspose(sig,rows,cols,col);
	} else {
		dtranspose(sig,rows,cols,col);
	}
}


void mdisplay(double *A, int row, int col) {
	int i,j;
	printf("\n MATRIX Order : %d X %d \n \n",row,col);
	
	for (i = 0; i < row; i++) {
		printf("R%d: ",i);
		for ( j = 0; j < col;j++) {
			printf("%g ",A[i*col + j]);
		}
		printf(":R%d \n",i);
	}
}

void madd(double* A, double* B, double* C,int rows,int cols) {
	int N,i;
	/*
	 * C = A + B . All matrices have identical dimensions rows X cols
	 */ 
	 
	N = rows * cols;
	
	for (i = 0; i < N; ++i) {
		C[i] = A[i] + B[i];
	}
}

void msub(double* A, double* B, double* C,int rows,int cols) {
	int N,i;
	/*
	 * C = A - B . All matrices have identical dimensions rows X cols
	 */ 
	 
	N = rows * cols;
	
	for (i = 0; i < N; ++i) {
		C[i] = A[i] - B[i];
	}
}

void scale(double *A, int rows, int cols, double alpha) {
	int N,i;
	/*
	 * A = alpha * A
	 * Matrix A is overwritten.
	 */ 
	 
	N = rows * cols;
	
	for (i = 0; i < N;++i) {
		A[i] = alpha * A[i];
	}
}

void nmult(double* A, double* B, double* C,int ra,int ca, int cb) {
	register int i,j,k;
	int u,v,t,rb;
	
	/*
	 * C = A * B , where A is a ra*ca matric while B is a rb*cb
	 * with ca = rb
	 * Matrix C is a ra*cb matrix
	 */ 
	 
	rb = ca;
	for (i = 0; i < ra; ++i) {
		for (j = 0; j < cb; ++j) {
			v = i * rb;
			u = i *cb;
			t = j + u;
			C[t] = 0.;
			for (k = 0; k < rb;++k) {
				C[t] += A[k + v] * B[j + k * cb];
			}
		}
	}

}

void tmult(double* A, double* B, double* C,int ra,int ca, int cb) {
	register int i,j,k;
	int u,v,t,rb;
	double *BT;
	BT = (double*) malloc(sizeof(double) * ca * cb);
	/*
	 * C = A * B , where A is a ra*ca matric while B is a rb*cb
	 * with ca = rb
	 * Matrix C is a ra*cb matrix
	 */ 
	 
	mtranspose(B,ca,cb,BT); 
	rb = ca;
	for (i = 0; i < ra; ++i) {
		for (j = 0; j < cb; ++j) {
			v = i * rb;
			u = i *cb;
			t = j + u;
			C[t] = 0.;
			for (k = 0; k < rb;++k) {
				C[t] += A[k + v] * BT[k + j * rb];
			}
		}
	}
	
	free(BT);

}


void recmult(double* A, double* B, double* C,int m,int n, int p,int sA,int sB, int sC) {
	int m2,n2,p2;
	register int i,j,k;
	int u,v,t;
	if (m + n + p <= CUTOFF) {
		for (i = 0; i < m; ++i) {
			for (j = 0; j < p; ++j) {
				v = i * sB;
				u = i * sC;
				t = j + u;
				for (k = 0; k < n;++k) {
					C[t] += A[k + v] * B[j + k * sC];
				}
			}
		}

		
	} else if (m >= n && m >= p) {
		m2 = (int) ceil((double) m / 2.0);
		recmult(A,B,C,m2,n,p,sA,sB,sC);
		recmult(A + m2*sB,B,C + m2*sC,m-m2,n,p,sA,sB,sC);
		
	} else if (n >= m && n >= p) {
		n2 = (int) ceil((double) n / 2.0);
		recmult(A,B,C,m,n2,p,sA,sB,sC);
		recmult(A+n2,B+n2*sC,C,m,n-n2,p,sA,sB,sC);
		
	} else if (p >= m && p >= n) {
		p2 = (int) ceil((double) p / 2.0);
		recmult(A,B,C,m,n,p2,sA,sB,sC);
		recmult(A,B+p2,C+p2,m,n,p-p2,sA,sB,sC);
	}
}

void rmult(double* A, double* B, double* C,int m,int n, int p) {
	int strA,strB,strC;
	int N;
	register int i;
	strA = m;
	strB = n;
	strC = p;
	N = m * p;
	for(i = 0; i < N; ++i) {
		C[i] = 0.;
	}
	
	recmult(A,B,C,m,n,p,strA,strB,strC);
	
}

int findrec(int *a, int *b, int *c) {
	int rec;
	double da,db,dc,mul;
	da = (double) *a;
	db = (double) *b;
	dc = (double) *c;
	rec = 0;
	mul = 1.;
	
	while (da + db + dc > (double) CUTOFF) {
		rec++;
		mul *= 2;
		da = ceil(da/2.);
		db = ceil(db/2.);
		dc = ceil(dc/2.);
	}
	*a = (int) da * mul;
	*b = (int) db * mul;
	*c = (int) dc * mul;
	
	return rec;
}

void add_zero_pad(double *X, int rows, int cols, int zrow, int zcol,double *Y) {
	int r,c,i,j,u,v;
	r = rows + zrow;
	c = cols + zcol;
	
	for (i = 0; i < rows;++i) {
		u = i*c;
		v = i * cols;
		for (j = 0; j < cols;++j) {
			Y[u + j] = X[v + j];
		}
		for (j = cols; j < c;++j) {
			Y[u + j] = 0.;
		}
	}
	
	for (i = rows; i < r;++i) {
		u = i*c;
		for(j = 0; j < c;++j) {
			Y[u + j] = 0.;
		}
	}
	
}

void remove_zero_pad(double *Y, int rows, int cols, int zrow, int zcol,double *Z) {
	int r,c,i,j,u,v;
	r = rows - zrow;
	c = cols - zcol;
	
	for (i = 0; i < r; ++i) {
		u = i * c;
		v = i * cols;
		for (j = 0; j < c; ++j) {
			Z[j + u] = Y[j + v];
		}
	}
}

void madd_stride(double* A, double* B, double* C,int rows,int cols,int sA,int sB,int sC) {
	int i,j,u,v,w;
	 
	for (i = 0; i < rows; ++i) {
		u = i * sC;
		v = i * sA;
		w = i * sB;
		for(j = 0; j < cols;j++) {
			C[j + u] = A[j + v] + B[j + w];
		}
	}
}

void msub_stride(double* A, double* B, double* C,int rows,int cols,int sA,int sB,int sC) {
	int i,j,u,v,w;
	 
	for (i = 0; i < rows; ++i) {
		u = i * sC;
		v = i * sA;
		w = i * sB;
		for(j = 0; j < cols;j++) {
			C[j + u] = A[j + v] - B[j + w];
		}
	}
}

void rmadd_stride(double* A, double* B, double* C,int rows,int cols,int p,int sA,int sB,int sC) {
	int i,j,u,v,w;
	if (rows + cols + p <= CUTOFF) {
		for (i = 0; i < rows; ++i) {
			u = i * sC;
			v = i * sA;
			w = i * sB;
			for(j = 0; j < cols;j++) {
				C[j + u] = A[j + v] + B[j + w];
			}
		}
		
	 } else {
		 rows/=2;cols/=2;p/=2;
		 rmadd_stride(A,B,C,rows,cols,p,sA,sB,sC);
		 rmadd_stride(A + cols,B + cols,C + cols,rows,cols,p,sA,sB,sC);
		 rmadd_stride(A + rows *sB,B + rows *sC,C + rows *sC,rows,cols,p,sA,sB,sC);
		 rmadd_stride(A + rows *sB + cols,B + rows *sC + cols,C + rows *sC + cols,rows,cols,p,sA,sB,sC);
	 }
}

void rmsub_stride(double* A, double* B, double* C,int rows,int cols,int p,int sA,int sB,int sC) {
	int i,j,u,v,w;
	if (rows + cols + p <= CUTOFF) {
		for (i = 0; i < rows; ++i) {
			u = i * sC;
			v = i * sA;
			w = i * sB;
			for(j = 0; j < cols;j++) {
				C[j + u] = A[j + v] - B[j + w];
			}
		}
		
	 } else {
		 rows/=2;cols/=2;p/=2;
		 rmsub_stride(A,B,C,rows,cols,p,sA,sB,sC);
		 rmsub_stride(A + cols,B + cols,C + cols,rows,cols,p,sA,sB,sC);
		 rmsub_stride(A + rows *sB,B + rows *sC,C + rows *sC,rows,cols,p,sA,sB,sC);
		 rmsub_stride(A + rows *sB + cols,B + rows *sC + cols,C + rows *sC + cols,rows,cols,p,sA,sB,sC);
	 }
}

void srecmult(double* A, double* B, double* C,int m,int n, int p,int sA,int sB, int sC) {
	register int i,j,k;
	int u,v,t;
	double sum;
	double *A1,*B1;
	double *a11,*a12,*a21,*a22;
	double *b11,*b12,*b21,*b22;
	double *c11,*c12,*c21,*c22;
	double *m1,*m2,*m3,*m4,*m5,*m6,*m7;
	int sm1,sm2,sm3,sm4,sm5,sm6,sm7;
	int sA1,sB1;
	if (m + n + p <= CUTOFF) {
		for (i = 0; i < m; ++i) {
			for (j = 0; j < p; ++j) {
				v = i * sA;
				u = i * sC;
				t = j + u;
				sum = 0.;
				for (k = 0; k < n;++k) {
					sum += A[k + v] * B[j + k * sB];
				}
				C[t] = sum;
			}
		}

		
	} else {
		m/=2;n/=2;p/=2;
		// A size mXn, C size mXp
		a11 = A;
		a12 = A + n;
		a21 = A + m * sA;
		a22 = A + n + m * sA;
		
		//B size nXp
		
		b11 = B;
		b12 = B + p;
		b21 = B + n * sB;
		b22 = B + p + n * sB;
		
		//C size mXp
		
		c11 = C;
		c12 = C + p;
		c21 = C + m * sC;
		c22 = C + p + m * sC;
		
		// m matrices have dimension m X p each. See http://en.wikipedia.org/wiki/Strassen_algorithm
		
		m1 = (double*) malloc(sizeof(double) *m * p);
		sm1 = p;
		
		m3 = (double*) malloc(sizeof(double) *m * p);
		sm3 = p;
		
		m4 = (double*) malloc(sizeof(double) *m * p);
		sm4 = p;
		
		m2 = c21;
		sm2 = sC;
		
		m5 = c12;
		sm5 = sC;
		
		m6 = c22;
		sm6 = sC;
		
		
		m7 = c11;
		sm7 = sC;
		
		//m1
		
		sA1 = n;
		sB1 = p;
		
		A1 = (double*) malloc(sizeof(double) * m * n);
		B1 = (double*) malloc(sizeof(double) * n * p);
		
		madd_stride(a11,a22,A1,m,n,sA,sA,sA1);
		
		madd_stride(b11,b22,B1,n,p,sB,sB,sB1);
		
		srecmult(A1,B1,m1,m,n,p,sA1,sB1,sm1);
		
		free(A1);
		free(B1);
		
		
		//m2
		
		A1 = (double*) malloc(sizeof(double) * m * n);
		
		madd_stride(a21,a22,A1,m,n,sA,sA,sA1);
				
		srecmult(A1,b11,m2,m,n,p,sA1,sB,sm2);
		
		free(A1);
		
		
		//m3
		
		B1 = (double*) malloc(sizeof(double) * n * p);
		//rmsub_stride(B + p,B + p + n * sC,B1,n,p,m,sC,sC,sC/2);
		msub_stride(b12,b22,B1,n,p,sB,sB,sB1);
		srecmult(a11,B1,m3,m,n,p,sA,sB1,sm3);
		
		free(B1);
		
		//m4
		
		B1 = (double*) malloc(sizeof(double) * n * p);
		//rmsub_stride(B + p,B + p + n * sC,B1,n,p,m,sC,sC,sC/2);
		msub_stride(b21,b11,B1,n,p,sB,sB,sB1);
		srecmult(a22,B1,m4,m,n,p,sA,sB1,sm4);
		
		free(B1);
		
		//m5
		
		A1 = (double*) malloc(sizeof(double) * m * n);
		
		madd_stride(a11,a12,A1,m,n,sA,sA,sA1);
				
		srecmult(A1,b22,m5,m,n,p,sA1,sB,sm5);
		
		free(A1);
		
		
		//m6
		
		A1 = (double*) malloc(sizeof(double) * m * n);
		B1 = (double*) malloc(sizeof(double) * n * p);
		
		msub_stride(a21,a11,A1,m,n,sA,sA,sA1);
		madd_stride(b11,b12,B1,n,p,sB,sB,sB1);
		srecmult(A1,B1,m6,m,n,p,sA1,sB1,sm6);
	
		free(A1);
		free(B1);
		
		//m7
		
		A1 = (double*) malloc(sizeof(double) * m * n);
		B1 = (double*) malloc(sizeof(double) * n * p);
		
		msub_stride(a12,a22,A1,m,n,sA,sA,sA1);
		madd_stride(b21,b22,B1,n,p,sB,sB,sB1);
		srecmult(A1,B1,m7,m,n,p,sA1,sB1,sm7);
	
		free(A1);
		free(B1);
		
		
		// c11
		
		A1 = (double*) malloc(sizeof(double) * m * p);
		sA1 = p;
		madd_stride(m1,m7,m7,m,p,sm1,sm7,sm7);
		msub_stride(m4,m5,A1,m,p,sm4,sm5,sA1);
		madd_stride(m7,A1,m7,m,p,sm7,sA1,sm7);
		
		free(A1);
		
		
		// c22
		
		A1 = (double*) malloc(sizeof(double) * m * p);
		sA1 = p;
		madd_stride(m1,m6,m6,m,p,sm1,sm6,sm6);
		msub_stride(m3,m2,A1,m,p,sm3,sm2,sA1);
		madd_stride(m6,A1,m6,m,p,sm6,sA1,sm6);
		
		free(A1);
		
		//c12 
		
		madd_stride(m3,m5,m5,m,p,sm3,sm5,sm5);
		
		//c21
		
		madd_stride(m4,m2,m2,m,p,sm4,sm2,sm2);
		
		free(m1);
		free(m3);
		free(m4);
	}
}

void smult(double* A, double* B, double* C,int m,int n, int p) {
	int a,b,c,nrec,i;
	double *X,*Y,*Z,*P,*A1,*B1;
	a = m;
	b = n;
	c = p;
	nrec = findrec(&a,&b,&c);
	X = (double*) malloc(sizeof(double) * a * b);
	Y = (double*) malloc(sizeof(double) * b * c);
	Z = (double*) malloc(sizeof(double) * a * c);
	P = (double*) malloc(sizeof(double) * (a/2) * (c/2));

	
	add_zero_pad(A,m,n,a-m,b-n,X);
	add_zero_pad(B,n,p,b-n,c-p,Y);

	srecmult(X,Y,Z,a,b,c,b,c,c);
	// Memory allocation needs work
	
	remove_zero_pad(Z,a,c,a-m,c-p,C);
	
	// free X,Y,Z
	free(X);
	free(Y);
	free(Z);
	free(P);
	
}

void mmult(double* A, double* B, double* C,int m,int n, int p) {
	if (m+n+p <= CUTOFF/2) {
		nmult(A,B,C,m,n,p);
	} else {
		smult(A,B,C,m,n,p);
	}
}

static int pludecomp(double *A,int N,int *ipiv) {
	int k,j,l,c1,c2,mind,tempi;
	double ld,mult,mval,temp;
	for(k=0;k < N;++k)
		ipiv[k] = k;
	
	for(k = 0; k < N-1; ++k) {
		//c2 = k*N;
		mval = fabs(A[k*N + k]);
		mind = k;
		for (j=k+1; j < N;++j) {
			if (mval < fabs(A[j*N + k])) {
				mval = A[j*N + k];
				mind = j;
			}
		}
		
		if ( mind != k) {
			c1 = k *N;
			c2 = mind * N;
			tempi = ipiv[mind];
			ipiv[mind] = ipiv[k];
			ipiv[k] = tempi;
			for (j = 0; j < N;j++) {
				temp = A[c1 + j];
				*(A + c1 + j) = *(A + c2 + j);
				*(A + c2 + j) = temp;
			}
		}
		c2 = k*N;
		ld = A[c2 + k];
		if (ld != 0.) {
			for (j = k+1; j < N; ++j) {
				c1 = j*N;
				mult = A[c1+k] /= ld;
				//printf("\n k %d j %d mult %lf \n",k,j,mult);
				for(l = k+1; l < N; ++l) {
					A[c1+l] -= mult * A[c2 + l];
				}
			}
		}
		
	}
	return 0;
	
}

void ludecomp(double *A,int N,int *ipiv) {
	pludecomp(A,N,ipiv);
}

void linsolve(double *A,int N,double *b,int *ipiv,double *x) {
	int i,j,c1,l;
	double *y;
	double sum,temp;
	
	y = (double*) malloc(sizeof(double) *N);
	/*
	 * Two step Solution L * U * x = b
	 * Let U*x = y
	 * Solve L * y = b for y (Forward Substitution
	 * Solve U * x = b for x (Back Substitution)
	 */ 
	for(i = 0; i < N;++i) {
		y[i] = 0.;
		x[i] = 0.;
		if ( A[i*N + i] == 0.) {
			printf("The Matrix system does not have a unique solution");
			exit(1);
		}
		//printf("\n B %d",ipiv[i]);
	}
	
	// Forward Substitution
	
	y[0] = b[ipiv[0]];
	for(i = 1; i < N; ++i) {
		sum = 0.;
		c1 = i*N;
		for(j = 0; j < i; ++j) {
			sum += y[j] * A[c1 + j];
		}
		y[i] = b[ipiv[i]] - sum;
	}
	
	// Back Substitution
	
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
}

void minverse(double *A,int N,double *ipiv,double *inv) {
	int i,j,stride;
	double *col,*x;
	
	col = (double*) malloc(sizeof(double) * N);
	x = (double*) malloc(sizeof(double) * N);
	
	for (i = 0; i < N; ++i) {
		col[i] = 0.;
		x[i] = 0.;
	}
	
	for (i = 0; i < N; ++i) {
		col[i] = 1.;
		linsolve(A,N,col,x,ipiv);
		stride = i;
		for(j = 0; j < N;++j) {
			inv[stride] = x[j];
			stride+= N;
		}
		col[i] = 0.;
	}
		
	free(x);
	free(col);
}
