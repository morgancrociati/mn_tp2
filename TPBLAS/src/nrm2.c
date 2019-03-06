#include "mnblas.h"
#include "complexe.h"

float mnblas_snrm2(const int N, const float *X, const int incX)
{
	register unsigned int i = 0;

	register float result = 0;

	for (; (i < N); i += incX)
	{
		result = result + X[i]*X[i];
	}

	return sqrt(result);
}

double mnblas_dnrm2(const int N, const double *X, const int incX)
{
	register unsigned int i = 0;

	register double result = 0;

	for (; (i < N); i += incX)
	{
		result = result + X[i]*X[i];
	}

	return sqrt(result);
}

float mnblas_scnrm2(const int N, const void *X, const int incX){
	register unsigned int i = 0;

	register float result;

	register complexe_float_t *_X = X;

	for (; (i < N); i += incX)
	{
		result = result + mult_complexe_float(_X,_X);
	}

	return sqrt(result);
}

double mnblas_dznrm2(const int N, const void *X, const int incX){
	register unsigned int i = 0;

	register double result;

	register complexe_double_t *_X = X;

	for (; (i < N); i += incX)
	{
		result = result + mult_complexe_float(_X,_X);
	}

	return sqrt(result);
}