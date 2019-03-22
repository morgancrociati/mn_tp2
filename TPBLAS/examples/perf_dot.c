#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"
#include "blas_test.h"

#define VECSIZE 4096

#define NB_FOIS 512

int main(int argc, char **argv)
{
	register unsigned long long start;
	register unsigned long long end;
	float res_s;
	double res_d;
	complexe_float_t res_c;
	complexe_double_t res_z;
	register unsigned long i;

	float *vec1 = new_s_full(VECSIZE);
	float *vec2 = new_s_full(VECSIZE);
	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		res_s = mncblas_sdot(VECSIZE, vec1, 1, vec2, 1);
	end = _rdtsc();
	calcul_flop("sdot ", NB_FOIS * 2 * VECSIZE, end - start);

	double *vec3 = new_d_full(VECSIZE);
	double *vec4 = new_d_full(VECSIZE);
	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		res_d = mncblas_ddot(VECSIZE, vec3, 1, vec4, 1);
	end = _rdtsc();
	calcul_flop("ddot ", NB_FOIS * 2 * VECSIZE, end - start);

	complexe_float_t *vec5 = new_c_full(VECSIZE);
	complexe_float_t *vec6 = new_c_full(VECSIZE);
	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		mncblas_cdotu_sub(VECSIZE, vec5, 1, vec6, 1, &res_c);
	end = _rdtsc();
	calcul_flop("cdotu_sub ", NB_FOIS * 8 * VECSIZE, end - start);

	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		mncblas_cdotc_sub(VECSIZE, vec5, 1, vec6, 1, &res_c);
	end = _rdtsc();
	calcul_flop("cdotc_sub ", NB_FOIS * 9 * VECSIZE, end - start);

	complexe_double_t *vec7 = new_z_full(VECSIZE);
	complexe_double_t *vec8 = new_z_full(VECSIZE);
	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		mncblas_zdotu_sub(VECSIZE, vec7, 1, vec8, 1, &res_z);
	end = _rdtsc();
	calcul_flop("zdotu_sub ", NB_FOIS * 8 * VECSIZE, end - start);

	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		mncblas_zdotc_sub(VECSIZE, vec7, 1, vec8, 1, &res_z);
	end = _rdtsc();
	calcul_flop("zdotc_sub ", NB_FOIS * 9 * VECSIZE, end - start);
}
