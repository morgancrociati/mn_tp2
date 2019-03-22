#ifndef _BLAS_TEST_H_
#define _BLAS_TEST_H_

#include <x86intrin.h>
//Pour _rdtsc();

#include "complexe.h"

#define TEST_SIZE 128

static const float duree_cycle = (float) 1 / (float) 2.6 ;
// duree du cycle en nano seconde 10^-9

void calcul_flop (char *message, int nb_operations_flottantes, unsigned long long int cycles);

float* new_s(unsigned int size);

double* new_d(unsigned int size);

complexe_float_t* new_c(unsigned int size);

complexe_double_t* new_z(unsigned int size);

float* new_s_full(unsigned int size);

double* new_d_full(unsigned int size);

complexe_float_t* new_c_full(unsigned int size);

complexe_double_t* new_z_full(unsigned int size);

#endif