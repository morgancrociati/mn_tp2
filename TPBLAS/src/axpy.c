#include "mnblas.h"
#include "complexe.h"

void mnblas_saxpy(const int N, const float alpha, const float *X,
				  const int incX, float *Y, const int incY)
{
	register unsigned int i = 0;
	register unsigned int j = 0;

	for (; ((i < N) && (j < N)); i += incX, j += incY)
	{
		Y[j] = Y[j] + alpha * X[i];
	}

	return;
}

void mnblas_daxpy(const int N, const double alpha, const double *X,
				  const int incX, double *Y, const int incY)
{
	register unsigned int i = 0;
	register unsigned int j = 0;

	for (; ((i < N) && (j < N)); i += incX, j += incY)
	{
		Y[j] = Y[j] + alpha * X[i];
	}

	return;
}

void mnblas_caxpy(const int N, const void *alpha, const void *X,
				  const int incX, void *Y, const int incY)
{
	register unsigned int i = 0;
	register unsigned int j = 0;

	register complexe_float_t *_X = X;
	register complexe_float_t *_Y = Y;
	register complexe_float_t *_alpha = alpha;

	for (; ((i < N) && (j < N)); i += incX, j += incY)
	{
		*_Y = add_complexe_float(*_Y, mult_complexe_float(_X[i], *_alpha));
	}

	return;
}

void mnblas_zaxpy(const int N, const void *alpha, const void *X,
				  const int incX, void *Y, const int incY)
{
	register unsigned int i = 0;
	register unsigned int j = 0;

	register complexe_double_t *_X = X;
	register complexe_double_t *_Y = Y;
	register complexe_double_t *_alpha = alpha;

	for (; ((i < N) && (j < N)); i += incX, j += incY)
	{
		*_Y = add_complexe_double(*_Y, mult_complexe_double(_X[i], *_alpha));
	}

	return;
}