#include "mnblas.h"
#include "complexe.h"

void mncblas_scopy(const int N, const float *X, const int incX,
				   float *Y, const int incY)
{
	register unsigned int i = 0;
	register unsigned int j = 0;

	for (; ((i < N) && (j < N)); i += incX, j += incY)
	{
		Y[j] = X[i];
	}

	return;
}

void mncblas_dcopy(const int N, const double *X, const int incX,
				   double *Y, const int incY)
{
	register unsigned int i = 0;
	register unsigned int j = 0;

	for (; ((i < N) && (j < N)); i += incX, j += incY)
	{
		Y[j] = X[i];
	}

	return;
}

void mncblas_ccopy(const int N, const void *X, const int incX,
				   void *Y, const int incY)
{
	register unsigned int i = 0;
	register unsigned int j = 0;

	register complexe_float_t *_X = X;
	register complexe_float_t *_Y = Y;

	for (; ((i < N) && (j < N)); i += incX, j += incY)
	{
		_Y[j] = _X[i];
	}

	return;
}

void mncblas_zcopy(const int N, const void *X, const int incX,
				   void *Y, const int incY)
{
	register unsigned int i = 0;
	register unsigned int j = 0;

	register complexe_double_t *_X = X;
	register complexe_double_t *_Y = Y;

	for (; ((i < N) && (j < N)); i += incX, j += incY)
	{
		_Y[j] = _X[i];
	}

	return;
}
