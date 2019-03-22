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
    float res_s;
    double res_d;
    complexe_float_t res_c;
    complexe_double_t res_z;
    int i;

    for (i = 0; i < NB_FOIS; i++)
    {
        float *vec1 = new_s_full(VECSIZE);
        float *vec2 = new_s_full(VECSIZE);

        start = _rdtsc();
        res_s = mncblas_sdot(VECSIZE, vec1, 1, vec2, 1);
        end = _rdtsc();

        printf("mncblas_sdot %d : res = %3.2f nombre de cycles: %Ld \n", i, res_s, end - start);
        calcul_flop("sdot ", 2 * VECSIZE, end - start);
    }

    for (i = 0; i < NB_FOIS; i++)
    {
        double *vec3 = new_d_full(VECSIZE);
        double *vec4 = new_d_full(VECSIZE);

        start = _rdtsc();
        res_d = mncblas_ddot(VECSIZE, vec3, 1, vec4, 1);
        end = _rdtsc();

        printf("mncblas_ddot %d : res = %3.2f nombre de cycles: %Ld \n", i, res_d, end - start);
        calcul_flop("ddot ", 2 * VECSIZE, end - start);
    }

    for (i = 0; i < NB_FOIS; i++)
    {
        complexe_float_t *vec5 = new_c_full(VECSIZE);
        complexe_float_t *vec6 = new_c_full(VECSIZE);

        start = _rdtsc();
        mncblas_cdotc_sub(VECSIZE, vec5, 1, vec6, 1, &res_c);
        end = _rdtsc();

        printf("mncblas_cdotc_sub %d : res = %3.2f nombre de cycles: %Ld \n", i, res_c, end - start);
        calcul_flop("cdotc_sub ", 2 * VECSIZE, end - start);
    }

    for (i = 0; i < NB_FOIS; i++)
    {
        complexe_float_t *vec5 = new_c_full(VECSIZE);
        complexe_float_t *vec6 = new_c_full(VECSIZE);

        start = _rdtsc();
        mncblas_cdotu_sub(VECSIZE, vec5, 1, vec6, 1, &res_c);
        end = _rdtsc();

        printf("mncblas_cdotu_sub %d : res = %3.2f nombre de cycles: %Ld \n", i, res_c, end - start);
        calcul_flop("cdotu_sub ", 2 * VECSIZE, end - start);
    }

    for (i = 0; i < NB_FOIS; i++)
    {
        complexe_double_t *vec7 = new_z_full(VECSIZE);
        complexe_double_t *vec8 = new_z_full(VECSIZE);

        start = _rdtsc();
        mncblas_zdotu_sub(VECSIZE, vec7, 1, vec8, 1, &res_z);
        end = _rdtsc();

        printf("mncblas_zdotu_sub %d : res = %3.2f nombre de cycles: %Ld \n", i, res_z, end - start);
        calcul_flop("zdotu_sub ", 2 * VECSIZE, end - start);
    }

    for (i = 0; i < NB_FOIS; i++)
    {
        complexe_double_t *vec7 = new_z_full(VECSIZE);
        complexe_double_t *vec8 = new_z_full(VECSIZE);

        start = _rdtsc();
        mncblas_zdotc_sub(VECSIZE, vec7, 1, vec8, 1, &res_z);
        end = _rdtsc();

        printf("mncblas_zdotc_sub %d : res = %3.2f nombre de cycles: %Ld \n", i, res_z, end - start);
        calcul_flop("zdotc_sub ", 2 * VECSIZE, end - start);
    }
}
