#include <stdio.h>
#include <assert.h>

#include "mnblas.h"
#include "complexe.h"
#include "blas_test.h"

#define MATSIZE 20
#define M MATSIZE
#define N MATSIZE
#define K MATSIZE

#define NB_FOIS 64

int main(void){
	unsigned long long i;

	float *vec1 = new_s_full(VECSIZE);
	float *vec2 = new_s_full(VECSIZE);
	float *A1 = new_s_full(VECSIZE * MATSIZE);
	for (i = 0; i < NB_FOIS; i++)
		mncblas_sgemv(MNCblasRowMajor, MNCblasNoTrans, VECSIZE, MATSIZE, 5.0, A1, 1, vec1, 1, 5.0, vec2, 1);
	free(vec1);
	free(vec2);
	free(A1);

	double *vec3 = new_d_full(VECSIZE);
	double *vec4 = new_d_full(VECSIZE);
	double *A2 = new_d_full(VECSIZE * MATSIZE);
	for (i = 0; i < NB_FOIS; i++)
		mncblas_dgemv(MNCblasRowMajor, MNCblasNoTrans, VECSIZE, MATSIZE, 5.0, A2, 1, vec3, 1, 5.0, vec4, 1);
	free(vec3);
	free(vec4);
	free(A2);

	complexe_float_t *vec5 = new_c_full(VECSIZE);
	complexe_float_t *vec6 = new_c_full(VECSIZE);
	complexe_float_t *A3 = new_c_full(VECSIZE * MATSIZE);
	complexe_float_t a1 = {2.0, 2.0};
	complexe_float_t b1 = {2.0, 2.0};
	for (i = 0; i < NB_FOIS; i++)
		mncblas_cgemv(MNCblasRowMajor, MNCblasNoTrans, VECSIZE, MATSIZE, &a1, A3, 1, vec5, 1, &b1, vec6, 1);
	free(vec5);
	free(vec6);
	free(A3);

	complexe_double_t *vec7 = new_z_full(VECSIZE);
	complexe_double_t *vec8 = new_z_full(VECSIZE);
	complexe_double_t *A4 = new_z_full(VECSIZE * MATSIZE);
	complexe_double_t a2 = {2.0, 2.0};
	complexe_double_t b2 = {2.0, 2.0};
	for (i = 0; i < NB_FOIS; i++)
		mncblas_zgemv(MNCblasRowMajor, MNCblasNoTrans, VECSIZE, MATSIZE, &a2, A4, 1, vec7, 1, &b2, vec8, 1);
	free(vec7);
	free(vec8);
	free(A4);

	printf("Test gemv: passed\n");

	return 1;
}