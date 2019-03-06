#include "mnblas.h"
#include "complexe.h"

void cblas_sgemm(const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE transa,
				 const CBLAS_TRANSPOSE transb, const MKL_INT m, const MKL_INT n,
				 const MKL_INT k, const float alpha, const float *a, const MKL_INT lda,
				 const float *b, const MKL_INT ldb, const float beta, float *c,
				 const MKL_INT ldc);

void cblas_dgemm(const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE transa,
				 const CBLAS_TRANSPOSE transb, const MKL_INT m, const MKL_INT n,
				 const MKL_INT k, const double alpha, const double *a, const MKL_INT lda,
				 const double *b, const MKL_INT ldb, const double beta, double *c,
				 const MKL_INT ldc);

void cblas_cgemm(const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE transa,
				 const CBLAS_TRANSPOSE transb, const MKL_INT m, const MKL_INT n,
				 const MKL_INT k, const void *alpha, const void *a, const MKL_INT lda,
				 const void *b, const MKL_INT ldb, const void *beta, void *c,
				 const MKL_INT ldc);

void cblas_zgemm(const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE transa,
				 const CBLAS_TRANSPOSE transb, const MKL_INT m, const MKL_INT n,
				 const MKL_INT k, const void *alpha, const void *a, const MKL_INT lda,
				 const void *b, const MKL_INT ldb, const void *beta, void *c,
				 const MKL_INT ldc);
