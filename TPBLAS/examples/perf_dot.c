#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define VECSIZE    1024

#define NB_FOIS    5

typedef float vfloat [VECSIZE] ;
typedef double vdouble [VECSIZE];
typedef complexe_float_t vcfloat [VECSIZE];
typedef complexe_double_t vcdouble [VECSIZE];

vfloat vec1, vec2 ;
vdouble vec3, vec4;
vcfloat vec5, vec6;
vcdouble vec7, vec8;

void vector_init_f (vfloat V, float x)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    V [i] = x ;

  return ;
}

void vector_init_d (vdouble V, double x)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    V [i] = x ;

  return ;
}

void vector_init_c (vcfloat V, complexe_float_t x)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    V [i] = x ;

  return ;
}

void vector_init_z (vcdouble V, complexe_double_t x)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    V [i] = x ;

  return ;
}

void vector_print_f (vfloat V)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    printf ("%f ", V[i]) ;
  printf ("\n") ;
  
  return ;
}

void vector_print_d (vdouble V)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    printf ("%f ", V[i]) ;
  printf ("\n") ;
  
  return ;
}

void vector_print_c (vcfloat V)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    printf ("%f + i%f ", V[i].real, V[i].imaginary) ;
  printf ("\n") ;
  
  return ;
}

void vector_print_z (vcdouble V)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    printf ("%f + i%f ", V[i].real, V[i].imaginary) ;
  printf ("\n") ;
  
  return ;
}

int main (int argc, char **argv)
{
 unsigned long long start, end ;
 float res_f ;
 double res_d;
 complexe_float_t res_c;
 complexe_double_t res_z;
 int i ;

 for (i = 0 ; i < NB_FOIS; i++)
   {
     float* vec1 = new_s_full(VECSIZE);
     float* vec2 = new_s_full(VECSIZE);

     start = _rdtsc () ;
        res_f = mncblas_sdot (VECSIZE, vec1, 1, vec2, 1) ;
     end = _rdtsc () ;
     
     printf ("mncblas_sdot %d : res = %3.2f nombre de cycles: %Ld \n", i, res_f, end-start) ;
     calcul_flop ("sdot ", 2 * VECSIZE, end-start) ;
   }
    for (i = 0 ; i < NB_FOIS; i++)
   {
     double* vec3 = new_d_full(VECSIZE);
     double* vec4 = new_d_full(VECSIZE);

     start = _rdtsc () ;
        res_f = mncblas_ddot (VECSIZE, vec3, 1, vec4, 1) ;
     end = _rdtsc () ;
     
     printf ("mncblas_ddot %d : res = %3.2f nombre de cycles: %Ld \n", i, res_f, end-start) ;
     calcul_flop ("ddot ", 2 * VECSIZE, end-start) ;
   }
    for (i = 0 ; i < NB_FOIS; i++)
   {
     complexe_float_t * vec5 = new_c_full(VECSIZE);
     complexe_float_t * vec6 = new_c_full(VECSIZE);

     start = _rdtsc () ;
        res_f = mncblas_cdot (VECSIZE, vec5, 1, vec6, 1) ;
     end = _rdtsc () ;
     
     printf ("mncblas_cdot %d : res = %3.2f nombre de cycles: %Ld \n", i, res_f, end-start) ;
     calcul_flop ("cdot ", 2 * VECSIZE, end-start) ;
   }
    for (i = 0 ; i < NB_FOIS; i++)
   {
     complexe_double_t *vec7 = new_z_full(VECSIZE);
     complexe_double_t *vec8 = new_z_full(VECSIZE);

     start = _rdtsc () ;
        res_f = mncblas_zdot (VECSIZE, vec7, 1, vec8, 1) ;
     end = _rdtsc () ;
     
     printf ("mncblas_zdot %d : res = %3.2f nombre de cycles: %Ld \n", i, res_f, end-start) ;
     calcul_flop ("zdot ", 2 * VECSIZE, end-start) ;
   }

}

