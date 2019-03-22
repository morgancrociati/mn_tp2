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
layout == rowMajor
TransA == MNCblasNoTrans
*/

void mncblas_sgemv(const MNCBLAS_LAYOUT layout,
				   const MNCBLAS_TRANSPOSE TransA, const int M, const int N,
				   const float alpha, const float *A, const int lda,
				   const float *X, const int incX, const float beta,
				   float *Y, const int incY)
{
	register unsigned int i;

	//TransA == MNCblasNoTrans)
	for (i = 0; (i < M); i += incY)
	{
		Y[i] = beta * Y[i];

		//layout == MNCblasRowMajor
		Y[i] += alpha * mncblas_sdot(N, A + i * N, 1, X, incX);
	}
}

void mncblas_dgemv(MNCBLAS_LAYOUT layout,
				   MNCBLAS_TRANSPOSE TransA, const int M, const int N,
				   const double alpha, const double *A, const int lda,
				   const double *X, const int incX, const double beta,
				   double *Y, const int incY)
{
	register unsigned int i;

	//TransA == MNCblasNoTrans)
	for (i = 0; (i < M); i += incY)
	{
		Y[i] = beta * Y[i];

		//layout == MNCblasRowMajor
		Y[i] += alpha * mncblas_ddot(N, A + i * N, 1, X, incX);
	}
}

void mncblas_cgemv(MNCBLAS_LAYOUT layout,
				   MNCBLAS_TRANSPOSE TransA, const int M, const int N,
				   const void *alpha, const void *A, const int lda,
				   const void *X, const int incX, const void *beta,
				   void *Y, const int incY)
{
	register unsigned int i;

	register const complexe_float_t *_alpha = alpha;
	register const complexe_float_t *_beta = beta;
	register const complexe_float_t *_X = X;
	register complexe_float_t *_Y = Y;
	register const complexe_float_t *_A = A;

	register complexe_float_t tmp;

	//TransA == MNCblasNoTrans)
	for (i = 0; (i < M); i += incY)
	{
		tmp.imaginary = _Y[i].imaginary;
		tmp.real = _Y[i].real;

		//layout == MNCblasRowMajor
		mncblas_cdotu_sub(N, _A + i * N, 1, _X, incX, _Y + i);

		_Y[i] = add_complexe_float(mult_complexe_float(tmp, *_beta), mult_complexe_float(_Y[i], *_alpha));
	}
}

void mncblas_zgemv(MNCBLAS_LAYOUT layout,
				   MNCBLAS_TRANSPOSE TransA, const int M, const int N,
				   const void *alpha, const void *A, const int lda,
				   const void *X, const int incX, const void *beta,
				   void *Y, const int incY)
{
	register unsigned int i;

	register const complexe_double_t *_alpha = alpha;
	register const complexe_double_t *_beta = beta;
	register const complexe_double_t *_X = X;
	register complexe_double_t *_Y = Y;
	register const complexe_double_t *_A = A;

	register complexe_double_t tmp;

	//TransA == MNCblasNoTrans)
	for (i = 0; (i < M); i += incY)
	{
		tmp.imaginary = _Y[i].imaginary;
		tmp.real = _Y[i].real;

		//layout == MNCblasRowMajor
		mncblas_zdotu_sub(N, _A + i * N, 1, _X, incX, _Y + i);

		_Y[i] = add_complexe_double(mult_complexe_double(tmp, *_beta), mult_complexe_double(_Y[i], *_alpha));
	}
}