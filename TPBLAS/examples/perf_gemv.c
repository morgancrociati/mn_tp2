#include <stdio.h>
#include <x86intrin.h>
#include <cblas.h>

#include "mnblas.h"
#include "complexe.h"
#include "blas_test.h"

#define VECSIZE 64
#define MATSIZE VECSIZE

#define NB_FOIS 512

void mncblas(void){
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
	free(vec1);
	free(vec2);
	free(A1);

	double *vec3 = new_d_full(VECSIZE);
	double *vec4 = new_d_full(VECSIZE);
	double *A2 = new_d_full(VECSIZE * MATSIZE);
	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		mncblas_dgemv(MNCblasRowMajor, MNCblasNoTrans, VECSIZE, MATSIZE, 5.0, A2, 1, vec3, 1, 5.0, vec4, 1);
	end = _rdtsc();
	calcul_flop("dgemv ", NB_FOIS * MATSIZE * (2 * VECSIZE + 3), end - start);
	free(vec3);
	free(vec4);
	free(A2);

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
	free(vec5);
	free(vec6);
	free(A3);

	complexe_double_t *vec7 = new_z_full(VECSIZE);
	complexe_double_t *vec8 = new_z_full(VECSIZE);
	complexe_double_t *A4 = new_z_full(VECSIZE * MATSIZE);
	complexe_double_t a2 = {2.0, 2.0};
	complexe_double_t b2 = {2.0, 2.0};
	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		mncblas_zgemv(MNCblasRowMajor, MNCblasNoTrans, VECSIZE, MATSIZE, &a2, A4, 1, vec7, 1, &b2, vec8, 1);
	end = _rdtsc();
	calcul_flop("cgemv ", NB_FOIS * MATSIZE * (8 * VECSIZE + 8), end - start);
	free(vec7);
	free(vec8);
	free(A4);
}

void cblas(void){
	unsigned long long start;
	unsigned long long end;
	unsigned long long i;

	float *vec1 = new_s_full(VECSIZE);
	float *vec2 = new_s_full(VECSIZE);
	float *A1 = new_s_full(VECSIZE * MATSIZE);
	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		cblas_sgemv(MNCblasRowMajor, MNCblasNoTrans, VECSIZE, MATSIZE, 5.0, A1, MATSIZE, vec1, 1, 5.0, vec2, 1);
	end = _rdtsc();
	calcul_flop("sgemv ", NB_FOIS * MATSIZE * (2 * VECSIZE + 3), end - start);
	free(vec1);
	free(vec2);
	free(A1);

	double *vec3 = new_d_full(VECSIZE);
	double *vec4 = new_d_full(VECSIZE);
	double *A2 = new_d_full(VECSIZE * MATSIZE);
	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		cblas_dgemv(MNCblasRowMajor, MNCblasNoTrans, VECSIZE, MATSIZE, 5.0, A2, MATSIZE, vec3, 1, 5.0, vec4, 1);
	end = _rdtsc();
	calcul_flop("dgemv ", NB_FOIS * MATSIZE * (2 * VECSIZE + 3), end - start);
	free(vec3);
	free(vec4);
	free(A2);

	complexe_float_t *vec5 = new_c_full(VECSIZE);
	complexe_float_t *vec6 = new_c_full(VECSIZE);
	complexe_float_t *A3 = new_c_full(VECSIZE * MATSIZE);
	complexe_float_t a1 = {2.0, 2.0};
	complexe_float_t b1 = {2.0, 2.0};
	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		cblas_cgemv(MNCblasRowMajor, MNCblasNoTrans, VECSIZE, MATSIZE, &a1, A3, MATSIZE, vec5, 1, &b1, vec6, 1);
	end = _rdtsc();
	calcul_flop("cgemv ", NB_FOIS * MATSIZE * (8 * VECSIZE + 8), end - start);
	free(vec5);
	free(vec6);
	free(A3);

	complexe_double_t *vec7 = new_z_full(VECSIZE);
	complexe_double_t *vec8 = new_z_full(VECSIZE);
	complexe_double_t *A4 = new_z_full(VECSIZE * MATSIZE);
	complexe_double_t a2 = {2.0, 2.0};
	complexe_double_t b2 = {2.0, 2.0};
	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		cblas_zgemv(MNCblasRowMajor, MNCblasNoTrans, VECSIZE, MATSIZE, &a2, A4, MATSIZE, vec7, 1, &b2, vec8, 1);
	end = _rdtsc();
	calcul_flop("zgemv ", NB_FOIS * MATSIZE * (8 * VECSIZE + 8), end - start);
}

int main(int argc, char **argv)
{
	printf("MNCBLAS:\n");
	mncblas();
	printf("\n");
	printf("CBLAS:\n");
	cblas();

	return 1;
}
