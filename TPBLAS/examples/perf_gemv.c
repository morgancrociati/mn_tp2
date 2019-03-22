#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define VECSIZE 1024
#define M 32
#define N 32

#define NB_FOIS 5

void vector_print_f(float *V)
{
    register unsigned int i;

    for (i = 0; i < VECSIZE; i++)
        printf("%f ", V[i]);
    printf("\n");

    return;
}

void vector_print_d(double *V)
{
    register unsigned int i;

    for (i = 0; i < VECSIZE; i++)
        printf("%f ", V[i]);
    printf("\n");

    return;
}

void vector_print_c(complexe_float_t *V)
{
    register unsigned int i;

    for (i = 0; i < VECSIZE; i++)
        printf("(%f + i%f) ", V[i].real, V[i].imaginary);
    printf("\n");

    return;
}

void vector_print_z(complexe_double_t *V)
{
    register unsigned int i;

    for (i = 0; i < VECSIZE; i++)
        printf("(%f + i%f) ", V[i].real, V[i].imaginary);
    printf("\n");

    return;
}

int main(int argc, char **argv)
{
    unsigned long long start, end;
    int i;

    for (i = 0; i < NB_FOIS; i++)
    {
        float *vec1 = new_s_full(VECSIZE);
        float *vec2 = new_s_full(VECSIZE);
        float *A = new_s_full(M * N);

        start = _rdtsc();
        mncblas_sgemv(MNCblasRowMajor, MNCblasNoTrans, M, N, 5.0, A, 1, vec1, 1, 5.0, vec2, 1);
        end = _rdtsc();

        printf("mncblas_sgemv %d : ", i);
        vector_print_f(vec2);
        printf("nombre de cycles: %Ld \n", end - start);
        calcul_flop("sgemv ", 2 * VECSIZE, end - start);
    }

    for (i = 0; i < NB_FOIS; i++)
    {
        double *vec3 = new_d_full(VECSIZE);
        double *vec4 = new_d_full(VECSIZE);
        double *A = new_d_full(M * N);

        start = _rdtsc();
        mncblas_dgemv(MNCblasRowMajor, MNCblasNoTrans, M, N, 5.0, A, 1, vec3, 1, 5.0, vec4, 1);
        end = _rdtsc();

        printf("mncblas_dgemv %d : ", i);
        vector_print_f(vec4);
        printf("nombre de cycles: %Ld \n", end - start);
        calcul_flop("dgemv ", 2 * VECSIZE, end - start);
    }

    for (i = 0; i < NB_FOIS; i++)
    {
        complexe_float_t *vec5 = new_c_full(VECSIZE);
        complexe_float_t *vec6 = new_c_full(VECSIZE);
        complexe_float_t *A = new_c_full(M * N);

        complexe_float_t *a, *b;
        a->imaginary = 2.0;
        a->real = 2.0;
        b->imaginary = 2.0;
        b->real = 2.0;

        start = _rdtsc();
        mncblas_cgemv(MNCblasRowMajor, MNCblasNoTrans, M, N, a, A, 1, vec5, 1, b, vec6, 1);
        end = _rdtsc();

        printf("mncblas_cgemv %d : ", i);
        vector_print_f(vec6);
        printf("nombre de cycles: %Ld \n", end - start);
        calcul_flop("cgemv ", 2 * VECSIZE, end - start);
    }

    for (i = 0; i < NB_FOIS; i++)
    {
        complexe_double_t *vec7 = new_z_full(VECSIZE);
        complexe_double_t *vec8 = new_z_full(VECSIZE);
        complexe_double_t *A = new_z_full(M * N);

        complexe_double_t *a, *b;
        a->imaginary = 2.0;
        a->real = 2.0;
        b->imaginary = 2.0;
        b->real = 2.0;

        start = _rdtsc();
        mncblas_zgemv(MNCblasRowMajor, MNCblasNoTrans, M, N, a, A, 1, vec7, 1, b, vec8, 1);
        end = _rdtsc();

        printf("mncblas_zgemv %d : ", i);
        vector_print_f(vec8);
        printf("nombre de cycles: %Ld \n", end - start);
        calcul_flop("zgemv ", 2 * VECSIZE, end - start);
    }
}
