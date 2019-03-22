#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "blas_test.h"

#define RANDOM_VALUE rand()/(RAND_MAX/2.0)

char rand_initialized = 0;

void init_rand(void){
	rand_initialized = 1;
	srand(time(NULL));
}

void calcul_flop(char *message, int nb_operations_flottantes, unsigned long long int cycles)
{
	printf("%s %d operations %5.3f GFLOP/s\n", message, nb_operations_flottantes, ((float)nb_operations_flottantes) / (((float)cycles) * duree_cycle));
}

float *new_s(unsigned int size)
{
	return (float *)malloc(sizeof(float) * size);
}

double *new_d(unsigned int size)
{
	return (double *)malloc(sizeof(double) * size);
}

complexe_float_t *new_c(unsigned int size)
{
	return (complexe_float_t *)malloc(sizeof(complexe_float_t) * size);
}

complexe_double_t *new_z(unsigned int size)
{
	return (complexe_double_t *)malloc(sizeof(complexe_double_t) * size);
}

float *new_s_full(unsigned int size)
{
	if(!rand_initialized)
		init_rand();
	register int i;
	register float *result = (float *)malloc(sizeof(float) * size);
	for (i = 0; i < size; i++)
	{
		result[i] = RANDOM_VALUE;
	}
	return result;
}

double *new_d_full(unsigned int size)
{
	if(!rand_initialized)
		init_rand();
	register int i;
	register double *result = (double *)malloc(sizeof(double) * size);
	for (i = 0; i < size; i++)
	{
		result[i] = RANDOM_VALUE;
	}
	return result;
}

complexe_float_t *new_c_full(unsigned int size)
{
	if(!rand_initialized)
		init_rand();
	register int i;
	register complexe_float_t *result = (complexe_float_t *)malloc(sizeof(complexe_float_t) * size);
	for (i = 0; i < size; i++)
	{
		result[i].real = RANDOM_VALUE;
		result[i].imaginary = RANDOM_VALUE;
	}
	return result;
}

complexe_double_t *new_z_full(unsigned int size)
{
	if(!rand_initialized)
		init_rand();
	register int i;
	register complexe_double_t *result = (complexe_double_t *)malloc(sizeof(complexe_double_t) * size);
	for (i = 0; i < size; i++)
	{
		result[i].real = RANDOM_VALUE;
		result[i].imaginary = RANDOM_VALUE;
	}
	return result;
}
