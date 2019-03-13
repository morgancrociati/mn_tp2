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

	return sqrtf(result);
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

<<<<<<< HEAD
float mnblas_scnrm2(const int N, const void *X, const int incX)
{
	register unsigned int i = 0;

	register complexe_float_t *_X = X;

	register float result = 0;

	for (; (i < N); i += incX)
	{
		result = result + _X[i].imaginary*_X[i].imaginary + _X[i].real*_X[i].real;
	}

	return sqrtf(result);
}

double mnblas_dznrm2(const int N, const void *X, const int incX)
{
	register unsigned int i = 0;

	register complexe_double_t *_X = X;

	register double result = 0;

	for (; (i < N); i += incX)
	{
		result = result + _X[i].imaginary*_X[i].imaginary + _X[i].real*_X[i].real;
=======
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
>>>>>>> b470923ae405a29be627b5ca1eaac99b100dc80f
	}

	return sqrt(result);
}