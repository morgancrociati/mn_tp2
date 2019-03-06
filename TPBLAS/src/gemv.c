#include "mnblas.h"
#include "complexe.h"

/*
a <- matrix 
(i, j) <- elem position 
m <- number of rows 
n <- number of column
*/
#define row_major_order(a, i, j, m) a[i * m + j]
#define column_major_order(a, i, j, n) a[i + j * n]

/*
if trans=CblasNoTrans, then y := alpha*A*x + beta*y;

if trans=CblasTrans, then y := alpha*A'*x + beta*y;

if trans=CblasConjTrans, then y := alpha *conjg(A')*x + beta*y.
*/

/*
use dot between A and X
followed by axpy ?

Maybe just dot and do it in line
*/

void mncblas_sgemv(const MNCBLAS_LAYOUT layout,
				   const MNCBLAS_TRANSPOSE TransA, const int M, const int N,
				   const float alpha, const float *A, const int lda,
				   const float *X, const int incX, const float beta,
				   float *Y, const int incY)
{
	register unsigned int i;
	register unsigned int j;

	if (TransA == MNCblasNoTrans)
	{
		for (i = 0; (i < M); i += incY)
		{
			Y[i] = beta * Y[i];
		}
	}
	else
	{
		for (i = 0; (i < N); i += incY)
		{
			Y[i] = beta * Y[i];
		}
	}

	if (TransA == MNCblasNoTrans)
	{
		for (i = 0; (i < M); i += incY)
		{
			for (j = 0; (j < M); j += incX)
			{
				if (layout == MNCblasRowMajor)
				{
					Y[i] = alpha * X[j] * mnblas_sasum(M*N, A + i*M, 1);
				}
				else
				{
					Y[i] = alpha * X[j] * mnblas_sasum(M*N, A + i, M);
				}
			}
		}
	}
	else
	{
		for (i = 0; (i < N); i += incY)
		{
			for (j = 0; (j < N); j += incX)
			{
				if (layout != MNCblasRowMajor)
				{
					Y[i] = alpha * X[j] * mnblas_sasum(M*N, A + i*M, 1);
				}
				else
				{
					Y[i] = alpha * X[j] * mnblas_sasum(M*N, A + i, M);
				}
			}
		}
	}
}

void mncblas_dgemv(MNCBLAS_LAYOUT layout,
				   MNCBLAS_TRANSPOSE TransA, const int M, const int N,
				   const double alpha, const double *A, const int lda,
				   const double *X, const int incX, const double beta,
				   double *Y, const int incY);

void mncblas_cgemv(MNCBLAS_LAYOUT layout,
				   MNCBLAS_TRANSPOSE TransA, const int M, const int N,
				   const void *alpha, const void *A, const int lda,
				   const void *X, const int incX, const void *beta,
				   void *Y, const int incY);

void mncblas_zgemv(MNCBLAS_LAYOUT layout,
				   MNCBLAS_TRANSPOSE TransA, const int M, const int N,
				   const void *alpha, const void *A, const int lda,
				   const void *X, const int incX, const void *beta,
				   void *Y, const int incY);