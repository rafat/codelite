
#ifndef CHOLESKY_H_
#define CHOLESKY_H_

#include "matrix.h"

#define BLOCKSIZE 32
#define CODE -1

#ifdef __cplusplus
extern "C" {
#endif

int cholu(double *A, int N);

int bcholu(double *A, int N);

int chol(double *A, int N);

void chold(double *A, int N);

#ifdef __cplusplus
}
#endif

#endif /* CHOLESKY_H_ */
