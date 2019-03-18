#include "mnblas.h"
#include "complexe.h"

void mncblas_sswap(const int N, float *X, const int incX,
				   float *Y, const int incY)
{
	register unsigned int i = 0;
	register unsigned int j = 0;

	register float save;

	for (; ((i < (1 + (N - 1) * abs(incX))) && (j < (1 + (N - 1) * abs(incY)))); i += incX, j += incY)
	{
		save = Y[j];
		Y[j] = X[i];
		X[i] = save;
	}

	return;
}

void mncblas_dswap(const int N, double *X, const int incX,
				   double *Y, const int incY)
{
	register unsigned int i = 0;
	register unsigned int j = 0;

	register double save;

	for (; ((i < (1 + (N - 1) * abs(incX))) && (j < (1 + (N - 1) * abs(incY)))); i += incX, j += incY)
	{
		save = Y[j];
		Y[j] = X[i];
		X[i] = save;
	}

	return;
}

void mncblas_cswap(const int N, void *X, const int incX,
				   void *Y, const int incY)
{
	register unsigned int i = 0;
	register unsigned int j = 0;

	register complexe_float_t *_X = X;
	register complexe_float_t *_Y = Y;

	register complexe_float_t save;

	for (; ((i < (1 + (N - 1) * abs(incX))) && (j < (1 + (N - 1) * abs(incY)))); i += incX, j += incY)
	{
		save = _Y[j];
		_Y[j] = _X[i];
		_X[i] = save;
	}

	return;
}

void mncblas_zswap(const int N, void *X, const int incX,
				   void *Y, const int incY)
{
	register unsigned int i = 0;
	register unsigned int j = 0;

	register complexe_double_t *_X = X;
	register complexe_double_t *_Y = Y;

	register complexe_double_t save;

	for (; ((i < (1 + (N - 1) * abs(incX))) && (j < (1 + (N - 1) * abs(incY)))); i += incX, j += incY)
	{
		save = _Y[j];
		_Y[j] = _X[i];
		_X[i] = save;
	}

	return;
}
