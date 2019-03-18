#include "mnblas.h"
#include "complexe.h"



	float mncblas_sdot(const int N, const float *X, const int incX,
					   const float *Y, const int incY)
{
	register unsigned int i = 0;
	register unsigned int j = 0;

	register float dot = 0.0;

	for (; ((i < (1 + (N - 1) * abs(incX))) && (j < (1 + (N - 1) * abs(incY)))); i += incX, j += incY)
	{
		dot = dot + X[i] * Y[j];
	}

	return dot;
}

double mncblas_ddot(const int N, const double *X, const int incX,
					const double *Y, const int incY)
{
	register unsigned int i = 0;
	register unsigned int j = 0;

	register double dot = 0.0;

	for (; ((i < (1 + (N - 1) * abs(incX))) && (j < (1 + (N - 1) * abs(incY)))); i += incX, j += incY)
	{
		dot = dot + X[i] * Y[j];
	}

	return dot;
}

void mncblas_cdotu_sub(const int N, const void *X, const int incX,
					   const void *Y, const int incY, void *dotu)
{
	register unsigned int i = 0;
	register unsigned int j = 0;

	register complexe_float_t *_X = X;
	register complexe_float_t *_Y = Y;

	register complexe_float_t *_dotu = dotu;
	_dotu->real = 0.0;
	_dotu->imaginary = 0.0;

	for (; ((i < (1 + (N - 1) * abs(incX))) && (j < (1 + (N - 1) * abs(incY)))); i += incX, j += incY)
	{
		*_dotu = add_complexe_float(*_dotu, mult_complexe_float(_X[i], _Y[j]));
	}

	return;
}

void mncblas_cdotc_sub(const int N, const void *X, const int incX,
					   const void *Y, const int incY, void *dotc)
{
	register unsigned int i = 0;
	register unsigned int j = 0;

	register complexe_float_t *_X = X;
	register complexe_float_t *_Y = Y;

	register complexe_float_t *_dotc = dotc;
	_dotc->real = 0.0;
	_dotc->imaginary = 0.0;

	for (; ((i < (1 + (N - 1) * abs(incX))) && (j < (1 + (N - 1) * abs(incY)))); i += incX, j += incY)
	{
		*_dotc = add_complexe_float(*_dotc, mult_complexe_float(conjg_complexe_float(_X[i]), _Y[j]));
	}

	return;
}

void mncblas_zdotu_sub(const int N, const void *X, const int incX,
					   const void *Y, const int incY, void *dotu)
{
	register unsigned int i = 0;
	register unsigned int j = 0;

	register complexe_double_t *_X = X;
	register complexe_double_t *_Y = Y;

	register complexe_double_t *_dotu = dotu;
	_dotu->real = 0.0;
	_dotu->imaginary = 0.0;

	for (; ((i < (1 + (N - 1) * abs(incX))) && (j < (1 + (N - 1) * abs(incY)))); i += incX, j += incY)
	{
		*_dotu = add_complexe_double(*_dotu, mult_complexe_double(_X[i], _Y[j]));
	}

	return;
}

void mncblas_zdotc_sub(const int N, const void *X, const int incX,
					   const void *Y, const int incY, void *dotc)
{
	register unsigned int i = 0;
	register unsigned int j = 0;

	register complexe_double_t *_X = X;
	register complexe_double_t *_Y = Y;

	register complexe_double_t *_dotc = dotc;
	_dotc->real = 0.0;
	_dotc->imaginary = 0.0;

	for (; ((i < (1 + (N - 1) * abs(incX))) && (j < (1 + (N - 1) * abs(incY)))); i += incX, j += incY)
	{
		*_dotc = add_complexe_double(*_dotc, mult_complexe_double(conjg_complexe_double(_X[i]), _Y[j]));
	}

	return;
}
