#include "mnblas.h"
#include "complexe.h"

#define abs(a) ((a) < (0) ? (-a) : (a))

CBLAS_INDEX mnblas_isamin(const int N, const float *X, const int incX)
{
	register unsigned int i = 0;

	register CBLAS_INDEX pos = 0;
	register float value = abs(X[0]); //Not mandatory //Might cause a bug if X is NULL

	register float tmp;

	for (; (i < (1 + (N - 1) * abs(incX))); i += incX)
	{
		tmp = abs(X[i]);
		if (tmp < value)
		{
			value = tmp;
			pos = i;
		}
	}

	return pos;
}

CBLAS_INDEX mnblas_idamin(const int N, const double *X, const int incX)
{
	register unsigned int i = 0;

	register CBLAS_INDEX pos = 0;
	register double value = abs(X[0]); //Not mandatory //Might cause a bug if X is NULL

	register double tmp;

	for (; (i < (1 + (N - 1) * abs(incX))); i += incX)
	{
		tmp = abs(X[i]);
		if (tmp < value)
		{
			value = tmp;
			pos = i;
		}
	}

	return pos;
}

CBLAS_INDEX mnblas_icamin(const int N, const void *X, const int incX)
{
	register unsigned int i = 0;

	register complexe_float_t *_X = X;

	register CBLAS_INDEX pos = 0;
	register float value = abs(_X[0].real) + abs(_X[0].imaginary); //Not mandatory //Might cause a bug if X is NULL

	register float tmp;

	for (; (i < (1 + (N - 1) * abs(incX))); i += incX)
	{
		tmp = abs(_X[i].real) + abs(_X[i].imaginary);
		if (tmp < value)
		{
			value = tmp;
			pos = i;
		}
	}

	return pos;
}

CBLAS_INDEX mnblas_izamin(const int N, const void *X, const int incX)
{
	register unsigned int i = 0;

	register complexe_double_t *_X = X;

	register CBLAS_INDEX pos = 0;
	register double value = abs(_X[0].real) + abs(_X[0].imaginary); //Not mandatory //Might cause a bug if X is NULL

	register double tmp;

	for (; (i < (1 + (N - 1) * abs(incX))); i += incX)
	{
		tmp = abs(_X[i].real) + abs(_X[i].imaginary);
		if (tmp < value)
		{
			value = tmp;
			pos = i;
		}
	}

	return pos;
}