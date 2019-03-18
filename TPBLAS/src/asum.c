#include "mnblas.h"
#include "complexe.h"

#define abs(a) ((a) < (0) ? (-a) : (a))

float mnblas_sasum(const int N, const float *X, const int incX)
{
	register unsigned int i = 0;

	register float sum = 0.0;

	for (; (i < (1 + (N - 1) * abs(incX))); i += incX)
	{
		sum = sum + abs(X[i]);
	}

	return sum;
}

double mnblas_dasum(const int N, const double *X, const int incX)
{
	register unsigned int i = 0;

	register double sum = 0.0;

	for (; (i < (1 + (N - 1) * abs(incX))); i += incX)
	{
		sum = sum + abs(X[i]);
	}

	return sum;
}

float mnblas_scasum(const int N, const void *X, const int incX)
{
	register unsigned int i = 0;

	register complexe_float_t *_X = X;

	register float sum = 0.0;

	for (; (i < (1 + (N - 1) * abs(incX))); i += incX)
	{
		sum = sum + abs(_X[i].real) + abs(_X[i].imaginary);
	}

	return sum;
}

double mnblas_dzasum(const int N, const void *X, const int incX)
{
	register unsigned int i = 0;

	register complexe_double_t *_X = X;

	register double sum = 0.0;

	for (; (i < (1 + (N - 1) * abs(incX))); i += incX)
	{
		sum = sum + abs(_X[i].real) + abs(_X[i].imaginary);
	}

	return sum;
}