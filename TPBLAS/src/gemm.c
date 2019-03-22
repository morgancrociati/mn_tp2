#include "mnblas.h"
#include "complexe.h"

void mncblas_sgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
				   MNCBLAS_TRANSPOSE TransB, const int M, const int N,
				   const int K, const float alpha, const float *A,
				   const int lda, const float *B, const int ldb,
				   const float beta, float *C, const int ldc)
{
	register unsigned int i;

	for (i = 0; (i < M * N); i++)
	{
		C[i] = alpha * mncblas_sdot(K, A + i * K, 1, B + i, N) + beta;
	}
}

void mncblas_dgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
				   MNCBLAS_TRANSPOSE TransB, const int M, const int N,
				   const int K, const double alpha, const double *A,
				   const int lda, const double *B, const int ldb,
				   const double beta, double *C, const int ldc)
{
	register unsigned int i;

	for (i = 0; (i < M * N); i++)
	{
		C[i] = alpha * mncblas_ddot(K, A + i * K, 1, B + i, N) + beta;
	}
}

void mncblas_cgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
				   MNCBLAS_TRANSPOSE TransB, const int M, const int N,
				   const int K, const void *alpha, const void *A,
				   const int lda, const void *B, const int ldb,
				   const void *beta, void *C, const int ldc)
{
	register unsigned int i;
	const register complexe_float_t *_A = A;
	const register complexe_float_t *_B = B;
	register complexe_float_t *_C = C;
	register complexe_float_t tmp;

	const register complexe_float_t *_beta = beta;
	const register complexe_float_t *_alpha = alpha;

	for (i = 0; (i < M * N); i++)
	{
		tmp.imaginary = _C[i].imaginary;
		tmp.real = _C[i].real;

		mncblas_cdotu_sub(K, A + i * K, 1, B + i, N, _C + i);
		_C[i] = mult_complexe_float(*_alpha, add_complexe_float(tmp, *_beta));
	}
}

void mncblas_zgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
				   MNCBLAS_TRANSPOSE TransB, const int M, const int N,
				   const int K, const void *alpha, const void *A,
				   const int lda, const void *B, const int ldb,
				   const void *beta, void *C, const int ldc)
{
	register unsigned int i;
	const register complexe_double_t *_A = A;
	const register complexe_double_t *_B = B;
	register complexe_double_t *_C = C;
	register complexe_double_t tmp;

	const register complexe_double_t *_beta = beta;
	const register complexe_double_t *_alpha = alpha;

	for (i = 0; (i < M * N); i++)
	{
		tmp.imaginary = _C[i].imaginary;
		tmp.real = _C[i].real;

		mncblas_zdotu_sub(K, A + i * K, 1, B + i, N, _C + i);
		_C[i] = mult_complexe_double(*_alpha, add_complexe_double(tmp, *_beta));
	}
}
