#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"
#include "blas_test.h"

#define VECSIZE 512
#define MATSIZE 512

#define NB_FOIS 128

int main(int argc, char **argv)
{
	unsigned long long start;
	unsigned long long end;
	unsigned long long i;

	float *vec1 = new_s_full(VECSIZE);
	float *vec2 = new_s_full(VECSIZE);
	float *A1 = new_s_full(VECSIZE * MATSIZE);
	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		mncblas_sgemv(MNCblasRowMajor, MNCblasNoTrans, VECSIZE, MATSIZE, 5.0, A1, 1, vec1, 1, 5.0, vec2, 1);
	end = _rdtsc();
	calcul_flop("sgemv ", NB_FOIS * MATSIZE * (2 * VECSIZE + 3), end - start);

	double *vec3 = new_d_full(VECSIZE);
	double *vec4 = new_d_full(VECSIZE);
	double *A2 = new_d_full(VECSIZE * MATSIZE);
	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		mncblas_dgemv(MNCblasRowMajor, MNCblasNoTrans, VECSIZE, MATSIZE, 5.0, A2, 1, vec3, 1, 5.0, vec4, 1);
	end = _rdtsc();
	calcul_flop("dgemv ", NB_FOIS * MATSIZE * (2 * VECSIZE + 3), end - start);

	complexe_float_t *vec5 = new_c_full(VECSIZE);
	complexe_float_t *vec6 = new_c_full(VECSIZE);
	complexe_float_t *A3 = new_c_full(VECSIZE * MATSIZE);
	complexe_float_t a1 = {2.0, 2.0};
	complexe_float_t b1 = {2.0, 2.0};
	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		mncblas_cgemv(MNCblasRowMajor, MNCblasNoTrans, VECSIZE, MATSIZE, &a1, A3, 1, vec5, 1, &b1, vec6, 1);
	end = _rdtsc();
	calcul_flop("cgemv ", NB_FOIS * MATSIZE * (8 * VECSIZE + 8), end - start);

	complexe_double_t *vec7 = new_z_full(VECSIZE);
	complexe_double_t *vec8 = new_z_full(VECSIZE);
	complexe_double_t *A4 = new_z_full(VECSIZE * MATSIZE);
	complexe_double_t a2 = {2.0, 2.0};
	complexe_double_t b2 = {2.0, 2.0};
	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		mncblas_zgemv(MNCblasRowMajor, MNCblasNoTrans, VECSIZE, MATSIZE, &a2, A4, 1, vec7, 1, &b2, vec8, 1);
	end = _rdtsc();
	calcul_flop("zgemv ", NB_FOIS * MATSIZE * (8 * VECSIZE + 8), end - start);
}
