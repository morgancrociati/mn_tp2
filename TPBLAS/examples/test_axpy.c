#include <stdio.h>
#include <assert.h>

#include "mnblas.h"
#include "complexe.h"
#include "blas_test.h"

int main(void){
	register int i;

	//Test float
	float* s1 = new_s_full(TEST_SIZE);
	float* s1bis = new_s(TEST_SIZE);
	float* s2 = new_s_full(TEST_SIZE);
	float* s2bis = new_s(TEST_SIZE);
	float salpha = 0;

	mncblas_scopy(TEST_SIZE, s1, 1, s1bis, 1);
	mncblas_scopy(TEST_SIZE, s2, 1, s2bis, 1);

	mnblas_saxpy(TEST_SIZE, salpha, s1, 1, s2, 1);
	mnblas_saxpy(TEST_SIZE, salpha, s2, 1, s1, 1);

	for(i = 0; i < TEST_SIZE; i++){
		assert(s1[i] == s1bis[i]);
		assert(s2[i] == s2bis[i]);
	}

	//Test double
	double* d1 = new_d_full(TEST_SIZE);
	double* d1bis = new_d(TEST_SIZE);
	double* d2 = new_d_full(TEST_SIZE);
	double* d2bis = new_d(TEST_SIZE);
	double dalpha = 0;

	mncblas_dcopy(TEST_SIZE, d1, 1, d1bis, 1);
	mncblas_dcopy(TEST_SIZE, d2, 1, d2bis, 1);

	mnblas_daxpy(TEST_SIZE, dalpha, d1, 1, d2, 1);
	mnblas_daxpy(TEST_SIZE, dalpha, d2, 1, d1, 1);

	for(i = 0; i < TEST_SIZE; i++){
		assert(d1[i] == d1bis[i]);
		assert(d2[i] == d2bis[i]);
	}

	//Test complexe float
	complexe_float_t* c1 = new_c_full(TEST_SIZE);
	complexe_float_t* c1bis = new_c(TEST_SIZE);
	complexe_float_t* c2 = new_c_full(TEST_SIZE);
	complexe_float_t* c2bis = new_c(TEST_SIZE);
	complexe_float_t calpha = {0, 0};

	mncblas_ccopy(TEST_SIZE, c1, 1, c1bis, 1);
	mncblas_ccopy(TEST_SIZE, c2, 1, c2bis, 1);

	mnblas_caxpy(TEST_SIZE, &calpha, c1, 1, c2, 1);
	mnblas_caxpy(TEST_SIZE, &calpha, c2, 1, c1, 1);

	for(i = 0; i < TEST_SIZE; i++){
		assert(equal_complexe_float(c1[i], c1bis[i]));
		assert(equal_complexe_float(c2[i], c2bis[i]));
	}

	//Test complexe double
	complexe_double_t* z1 = new_z_full(TEST_SIZE);
	complexe_double_t* z1bis = new_z(TEST_SIZE);
	complexe_double_t* z2 = new_z_full(TEST_SIZE);
	complexe_double_t* z2bis = new_z(TEST_SIZE);
	complexe_double_t zalpha = {0, 0};

	mncblas_zcopy(TEST_SIZE, z1, 1, z1bis, 1);
	mncblas_zcopy(TEST_SIZE, z2, 1, z2bis, 1);

	mnblas_zaxpy(TEST_SIZE, &zalpha, z1, 1, z2, 1);
	mnblas_zaxpy(TEST_SIZE, &zalpha, z2, 1, z1, 1);

	for(i = 0; i < TEST_SIZE; i++){
		assert(equal_complexe_double(z1[i], z1bis[i]));
		assert(equal_complexe_double(z2[i], z2bis[i]));
	}

	printf("Test axpy: passed\n");

	return 1;
}