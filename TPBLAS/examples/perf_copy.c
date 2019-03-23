#include <stdio.h>
#include <x86intrin.h>
#include <cblas.h>
#include "mnblas.h"
#include "complexe.h"
#include "blas_test.h"

#define VECSIZE 4096

#define NB_FOIS 512

int main(int argc, char **argv)
{
	register unsigned long long start;
	register unsigned long long end;
	register unsigned long i;

	printf("MNCBLAS:\n");
	//Perf float
	float *vec1 = new_s_full(VECSIZE);
	float *vec2 = new_s(VECSIZE);
	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		mncblas_scopy(VECSIZE, vec1, 1, vec2, 1);
	end = _rdtsc();
	calcul_octet("scopy ", NB_FOIS * VECSIZE * sizeof(float), end - start);

	//Perf double
	double *vec3 = new_d_full(VECSIZE);
	double *vec4 = new_d(VECSIZE);
	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		mncblas_dcopy(VECSIZE, vec3, 1, vec4, 1);
	end = _rdtsc();
	calcul_octet("dcopy ", NB_FOIS * VECSIZE * sizeof(double), end - start);

	//Perf complexe float
	complexe_float_t *vec5 = new_c_full(VECSIZE);
	complexe_float_t *vec6 = new_c(VECSIZE);
	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		mncblas_ccopy(VECSIZE, vec5, 1, vec6, 1);
	end = _rdtsc();
	calcul_octet("ccopy ", NB_FOIS * 2 * VECSIZE * sizeof(complexe_float_t), end - start);

	//Perf complexe double
	complexe_double_t *vec7 = new_z_full(VECSIZE);
	complexe_double_t *vec8 = new_z(VECSIZE);
	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		mncblas_zcopy(VECSIZE, vec7, 1, vec8, 1);
	end = _rdtsc();
	calcul_octet("zcopy ", NB_FOIS * 2 * VECSIZE * sizeof(complexe_double_t), end - start);

//----------------------------------------------------------------------------------------
	printf("\nCBLAS:\n");
	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		cblas_scopy(VECSIZE, vec1, 1, vec2, 1);
	end = _rdtsc();
	calcul_octet("scopy ", NB_FOIS * VECSIZE * sizeof(float), end - start);


	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		cblas_dcopy(VECSIZE, vec3, 1, vec4, 1);
	end = _rdtsc();
	calcul_octet("dcopy ", NB_FOIS * VECSIZE * sizeof(double), end - start);


	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		cblas_ccopy(VECSIZE, vec5, 1, vec6, 1);
	end = _rdtsc();
	calcul_octet("ccopy ", NB_FOIS * 2 * VECSIZE * sizeof(complexe_float_t), end - start);


	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		cblas_zcopy(VECSIZE, vec7, 1, vec8, 1);
	end = _rdtsc();
	calcul_octet("zcopy ", NB_FOIS * 2 * VECSIZE * sizeof(complexe_double_t), end - start);

	return 0;
}
