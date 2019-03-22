#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define VECSIZE 1024

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
        printf("%f + i%f ", V[i].real, V[i].imaginary);
    printf("\n");

    return;
}

void vector_print_z(complexe_double_t *V)
{
    register unsigned int i;

    for (i = 0; i < VECSIZE; i++)
        printf("%f + i%f ", V[i].real, V[i].imaginary);
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

        start = _rdtsc();
        mncblas_scopy(VECSIZE, vec1, 1, vec2, 1);
        end = _rdtsc();

        printf("mncblas_scopy %d : ", i);
        vector_print_f(vec2);
        printf("nombre de cycles: %Ld \n", end - start);
        calcul_flop("scopy ", 2 * VECSIZE, end - start);
    }
    for (i = 0; i < NB_FOIS; i++)
    {
        double *vec3 = new_d_full(VECSIZE);
        double *vec4 = new_d_full(VECSIZE);

        start = _rdtsc();
        mncblas_dcopy(VECSIZE, vec3, 1, vec4, 1);
        end = _rdtsc();

        printf("mncblas_dcopy %d : ", i);
        vector_print_f(vec4);
        printf("nombre de cycles: %Ld \n", end - start);
        calcul_flop("dcopy ", 2 * VECSIZE, end - start);
    }
    for (i = 0; i < NB_FOIS; i++)
    {
        complexe_float_t *vec5 = new_c_full(VECSIZE);
        complexe_float_t *vec6 = new_c_full(VECSIZE);

        start = _rdtsc();
        mncblas_ccopy(VECSIZE, vec5, 1, vec6, 1);
        end = _rdtsc();

        printf("mncblas_ccopy %d : ", i);
        vector_print_f(vec6);
        printf("nombre de cycles: %Ld \n", end - start);
        calcul_flop("ccopy ", 2 * VECSIZE, end - start);
    }
    for (i = 0; i < NB_FOIS; i++)
    {
        complexe_double_t *vec7 = new_z_full(VECSIZE);
        complexe_double_t *vec8 = new_z_full(VECSIZE);

        start = _rdtsc();
        mncblas_zcopy(VECSIZE, vec7, 1, vec8, 1);
        end = _rdtsc();

        printf("mncblas_zcopy %d : ", i);
        vector_print_f(vec8);
        printf("nombre de cycles: %Ld \n", end - start);
        calcul_flop("zcopy ", 2 * VECSIZE, end - start);
    }
}
