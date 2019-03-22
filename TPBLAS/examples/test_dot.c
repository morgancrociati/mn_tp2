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

	
	//Test complexe float
	complexe_float_t* c1 = new_c_full(TEST_SIZE);
	complexe_float_t* c2 = new_c_full(TEST_SIZE);
	complexe_float_t cres1;
	complexe_float_t cres2;

	mncblas_cdotu_sub(TEST_SIZE, c1, 1, c2, 1, &cres1);
	mncblas_cdotu_sub(TEST_SIZE, c2, 1, c1, 1, &cres2);
	assert(equal_complexe_float(cres1, cres2));

	/*
	mncblas_cdotc_sub(TEST_SIZE, c1, 1, c2, 1, &cres1);
	mncblas_cdotc_sub(TEST_SIZE, c2, 1, c1, 1, &cres2);
	assert(equal_complexe_float(cres1, cres2));
	*/

	//Test complexe double
	complexe_double_t* z1 = new_z_full(TEST_SIZE);
	complexe_double_t* z2 = new_z_full(TEST_SIZE);
	complexe_double_t zres1;
	complexe_double_t zres2;

	mncblas_zdotu_sub(TEST_SIZE, z1, 1, z2, 1, &zres1);
	mncblas_zdotu_sub(TEST_SIZE, z2, 1, z1, 1, &zres2);
	assert(equal_complexe_double(zres1, zres2));

	/*
	mncblas_zdotc_sub(TEST_SIZE, z1, 1, z2, 1, &zres1);
	mncblas_zdotc_sub(TEST_SIZE, z2, 1, z1, 1, &zres2);
	assert(equal_complexe_double(zres1, zres2));
	*/

	printf("Test dot: passed\n");

	return 1;
}
