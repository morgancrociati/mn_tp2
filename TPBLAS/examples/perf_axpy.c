#include <stdio.h>
#include <x86intrin.h>
#include <cblas.h>
#include "mnblas.h"
#include "complexe.h"
#include "blas_test.h"

#define VECSIZE 1024

#define NB_FOIS 512

int main(int argc, char **argv)
{
	register unsigned long long start;
	register unsigned long long end;
	register unsigned long i;

	printf("MNCBLAS:\n");
	float *vec1 = new_s_full(VECSIZE);
	float *vec2 = new_s_full(VECSIZE);
	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		mnblas_saxpy(VECSIZE, 5.0, vec1, 1, vec2, 1);
	end = _rdtsc();
	calcul_flop("saxpy ", NB_FOIS*2 * VECSIZE, end - start);

	double *vec3 = new_d_full(VECSIZE);
	double *vec4 = new_d_full(VECSIZE);
	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		mnblas_daxpy(VECSIZE, 5.0, vec3, 1, vec4, 1);
	end = _rdtsc();
	calcul_flop("daxpy ", NB_FOIS*2 * VECSIZE, end - start);

	complexe_float_t *vec5 = new_c_full(VECSIZE);
	complexe_float_t *vec6 = new_c_full(VECSIZE);
	complexe_float_t a1 = {2.0, 2.0};
	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		mnblas_caxpy(VECSIZE, &a1, vec5, 1, vec6, 1);
	end = _rdtsc();
	calcul_flop("caxpy ", NB_FOIS*8 * VECSIZE, end - start);

	complexe_double_t *vec7 = new_z_full(VECSIZE);
	complexe_double_t *vec8 = new_z_full(VECSIZE);
	complexe_double_t a2 = {2.0, 2.0};
	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		mnblas_zaxpy(VECSIZE, &a2, vec7, 1, vec8, 1);
	end = _rdtsc();
	calcul_flop("zaxpy ", NB_FOIS*8 * VECSIZE, end - start);

//-----------------------------------------------------------------------------
	printf("\nCBLAS:\n");
	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		cblas_saxpy(VECSIZE, 5.0, vec1, 1, vec2, 1);
	end = _rdtsc();
	calcul_flop("saxpy ", NB_FOIS*2 * VECSIZE, end - start);

	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		cblas_daxpy(VECSIZE, 5.0, vec3, 1, vec4, 1);
	end = _rdtsc();
	calcul_flop("daxpy ", NB_FOIS*2 * VECSIZE, end - start);


	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		cblas_caxpy(VECSIZE, &a1, vec5, 1, vec6, 1);
	end = _rdtsc();
	calcul_flop("caxpy ", NB_FOIS*8 * VECSIZE, end - start);


	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		cblas_zaxpy(VECSIZE, &a2, vec7, 1, vec8, 1);
	end = _rdtsc();
	calcul_flop("zaxpy ", NB_FOIS*8 * VECSIZE, end - start);


	return 0;
}
