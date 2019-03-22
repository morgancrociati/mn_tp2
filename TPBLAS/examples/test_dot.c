#include <stdio.h>
#include <assert.h>

#include "mnblas.h"
#include "complexe.h"
#include "blas_test.h"

int main(void){
	//Test float
	float* s1 = new_s_full(TEST_SIZE);
	float* s2 = new_s_full(TEST_SIZE);

	assert(mncblas_sdot(TEST_SIZE, s1, 1, s2, 1) == mncblas_sdot(TEST_SIZE, s2, 1, s1, 1));

	//Test double
	double* d1 = new_d_full(TEST_SIZE);
	double* d2 = new_d_full(TEST_SIZE);

	assert(mncblas_ddot(TEST_SIZE, d1, 1, d2, 1) == mncblas_ddot(TEST_SIZE, d2, 1, d1, 1));

	/*
	//Test complexe float
	complexe_float_t* c1 = new_c_full(TEST_SIZE);
	complexe_float_t* c2 = new_c_full(TEST_SIZE);

	assert(mncblas_cdot(TEST_SIZE, c1, 1, c2, 1) == mncblas_cdot(TEST_SIZE, c2, 1, c1, 1));

	//Test complexe double
	complexe_double_t* z1 = new_z_full(TEST_SIZE);
	complexe_double_t* z2 = new_z_full(TEST_SIZE);

	assert(mncblas_zdot(TEST_SIZE, z1, 1, z2, 1) == mncblas_zdot(TEST_SIZE, z2, 1, z1, 1));
	*/

	printf("Test dot: passed\n");

	return 1;
}