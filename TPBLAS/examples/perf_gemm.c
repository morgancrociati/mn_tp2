#include <stdio.h>
#include <x86intrin.h>
#include <cblas.h>

#include "mnblas.h"
#include "complexe.h"
#include "blas_test.h"

#define MATSIZE 20
#define M MATSIZE
#define N MATSIZE
#define K MATSIZE

#define NB_FOIS 256

void mncblas(void){
	unsigned long long start;
	unsigned long long end;
	int i;

	float *C1 = new_s_full(M * N);
	float *B1 = new_s_full(M * K);
	float *A1 = new_s_full(K * N);
	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		mncblas_sgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, M, N, K, 5.0, A1, 1, B1, 1, 5.0, C1, 1);
	end = _rdtsc();
	calcul_flop("sgemm ", NB_FOIS*M*N*(K + 2), end - start);
	free(C1);
	free(B1);
	free(A1);

	double *C2 = new_d_full(M * N);
	double *B2 = new_d_full(M * K);
	double *A2 = new_d_full(K * N);
	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		mncblas_dgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, M, N, K, 5.0, A2, 1, B2, 1, 5.0, C2, 1);
	end = _rdtsc();
	calcul_flop("dgemm ", NB_FOIS*M*N*(K + 2), end - start);
	free(C2);
	free(B2);
	free(A2);

	complexe_float_t *C3 = new_c_full(M * N);
	complexe_float_t *B3 = new_c_full(M * K);
	complexe_float_t *A3 = new_c_full(K * N);
	complexe_float_t a1 = {2.0, 2.0};
	complexe_float_t b1 = {2.0, 2.0};
	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		mncblas_cgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, M, N, K, &a1, A3, 1, B3, 1, &b1, C3, 1);
	end = _rdtsc();
	calcul_flop("cgemm ", NB_FOIS*M*N*(K + 8), end - start);
	free(C3);
	free(B3);
	free(A3);

	complexe_double_t *C4 = new_z_full(M * N);
	complexe_double_t *B4 = new_z_full(M * K);
	complexe_double_t *A4 = new_z_full(K * N);
	complexe_double_t a2 = {2.0, 2.0};
	complexe_double_t b2 = {2.0, 2.0};
	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		mncblas_zgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, M, N, K, &a2, A4, 1, B4, 1, &b2, C4, 1);
	end = _rdtsc();
	calcul_flop("zgemm ", NB_FOIS*M*N*(K + 8), end - start);
	free(C4);
	free(B4);
	free(A4);
}

void cblas(void){
	unsigned long long start;
	unsigned long long end;
	int i;

	float *C1 = new_s_full(M * N);
	float *B1 = new_s_full(M * K);
	float *A1 = new_s_full(K * N);
	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		cblas_sgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, M, N, K, 5.0, A1, MATSIZE, B1, MATSIZE, 5.0, C1, MATSIZE);
	end = _rdtsc();
	calcul_flop("sgemm ", NB_FOIS*M*N*(K + 2), end - start);
	free(C1);
	free(B1);
	free(A1);

	double *C2 = new_d_full(M * N);
	double *B2 = new_d_full(M * K);
	double *A2 = new_d_full(K * N);
	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		cblas_dgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, M, N, K, 5.0, A2, MATSIZE, B2, MATSIZE, 5.0, C2, MATSIZE);
	end = _rdtsc();
	calcul_flop("dgemm ", NB_FOIS*M*N*(K + 2), end - start);
	free(C2);
	free(B2);
	free(A2);

	complexe_float_t *C3 = new_c_full(M * N);
	complexe_float_t *B3 = new_c_full(M * K);
	complexe_float_t *A3 = new_c_full(K * N);
	complexe_float_t a1 = {2.0, 2.0};
	complexe_float_t b1 = {2.0, 2.0};
	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		cblas_cgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, M, N, K, &a1, A3, MATSIZE, B3, MATSIZE, &b1, C3, MATSIZE);
	end = _rdtsc();
	calcul_flop("cgemm ", NB_FOIS*M*N*(K + 8), end - start);
	free(C3);
	free(B3);
	free(A3);

	complexe_double_t *C4 = new_z_full(M * N);
	complexe_double_t *B4 = new_z_full(M * K);
	complexe_double_t *A4 = new_z_full(K * N);
	complexe_double_t a2 = {2.0, 2.0};
	complexe_double_t b2 = {2.0, 2.0};
	start = _rdtsc();
	for (i = 0; i < NB_FOIS; i++)
		cblas_zgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, M, N, K, &a2, A4, MATSIZE, B4, MATSIZE, &b2, C4, MATSIZE);
	end = _rdtsc();
	calcul_flop("zgemm ", NB_FOIS*M*N*(K + 8), end - start);
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
