#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define VECSIZE 1024
#define M 32
#define N 32
#define K 32

#define NB_FOIS 5

void vector_print_f(float *V)
{
    register unsigned int i;
    register unsigned int j;

    for (i = 0; i < M; i++){
        for(j = 0; j < N; j < 0){
            printf("%f ", V[i]);
        }
        printf("\n");
    }

    printf("\n");

    return;
}

void vector_print_d(double *V)
{
    register unsigned int i;
    register unsigned int j;

    for (i = 0; i < M; i++){
        for(j = 0; j < N; j < 0){
            printf("%f ", V[i]);
        }
        printf("\n");
    }

    printf("\n");

    return;
}

void vector_print_c(complexe_float_t *V)
{
    register unsigned int i;
    register unsigned int j;

    for (i = 0; i < M; i++){
        for(j = 0; j < N; j < 0){
            printf("(%f + i%f) ", V[i].real, V[i].imaginary);
        }
        printf("\n");
    }

    printf("\n");

    return;
}

void vector_print_z(complexe_double_t *V)
{
    register unsigned int i;
    register unsigned int j;

    for (i = 0; i < M; i++){
        for(j = 0; j < N; j < 0){
            printf("(%f + i%f) ", V[i].real, V[i].imaginary);
        }
        printf("\n");
    }

    return;
}

int main(int argc, char **argv)
{
    unsigned long long start, end;
    int i;

    for (i = 0; i < NB_FOIS; i++)
    {
        float *C = new_s_full(M * N);
        float *B = new_s_full(M * K);
        float *A = new_s_full(K * N);

        start = _rdtsc();
        mncblas_sgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, M, N, K, 5.0, A, 1, B, 1, 5.0, C, 1);
        end = _rdtsc();

        printf("mncblas_sgemm %d : ", i);
        vector_print_f(C);
        printf("nombre de cycles: %Ld \n", end - start);
        calcul_flop("sgemm ", 2 * VECSIZE, end - start);
    }

    for (i = 0; i < NB_FOIS; i++)
    {
        double *C = new_d_full(M * N);
        double *B = new_d_full(M * K);
        double *A = new_d_full(K * N);

        start = _rdtsc();
        mncblas_dgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, M, N, K, 5.0, A, 1, B, 1, 5.0, C, 1);
        end = _rdtsc();

        printf("mncblas_dgemm %d : ", i);
        vector_print_f(C);
        printf("nombre de cycles: %Ld \n", end - start);
        calcul_flop("dgemm ", 2 * VECSIZE, end - start);
    }

    for (i = 0; i < NB_FOIS; i++)
    {
        complexe_float_t *C = new_c_full(M * N);
        complexe_float_t *B = new_c_full(M * K);
        complexe_float_t *A = new_c_full(K * N);

        complexe_float_t *a, *b;
        a->imaginary = 2.0;
        a->real = 2.0;
        b->imaginary = 2.0;
        b->real = 2.0;

        start = _rdtsc();
        mncblas_cgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, M, N, K, a, A, 1, B, 1, b, C, 1);
        end = _rdtsc();

        printf("mncblas_cgemm %d : ", i);
        vector_print_f(C);
        printf("nombre de cycles: %Ld \n", end - start);
        calcul_flop("cgemm ", 2 * VECSIZE, end - start);
    }

    for (i = 0; i < NB_FOIS; i++)
    {
        complexe_double_t *C = new_z_full(M * N);
        complexe_double_t *B = new_z_full(M * K);
        complexe_double_t *A = new_z_full(K * N);

        complexe_double_t *a, *b;
        a->imaginary = 2.0;
        a->real = 2.0;
        b->imaginary = 2.0;
        b->real = 2.0;

        start = _rdtsc();
        mncblas_zgemm(MNCblasRowMajor, MNCblasNoTrans, MNCblasNoTrans, M, N, K, a, A, 1, B, 1, b, C, 1);
        end = _rdtsc();

        printf("mncblas_zgemm %d : ", i);
        vector_print_f(C);
        printf("nombre de cycles: %Ld \n", end - start);
        calcul_flop("zgemm ", 2 * VECSIZE, end - start);
    }
}
