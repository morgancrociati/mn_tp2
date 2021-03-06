#ifndef _COMPLEXE_H_
#define _COMPLEXE_H_

typedef struct c
{
	float real;
	float imaginary;
} complexe_float_t;

typedef struct z
{
	double real;
	double imaginary;
} complexe_double_t;

complexe_float_t add_complexe_float(const complexe_float_t c1, const complexe_float_t c2);

complexe_double_t add_complexe_double(const complexe_double_t c1, const complexe_double_t c2);

complexe_float_t mult_complexe_float(const complexe_float_t c1, const complexe_float_t c2);

complexe_double_t mult_complexe_double(const complexe_double_t c1, const complexe_double_t c2);

complexe_float_t conjg_complexe_float(const complexe_float_t c);

complexe_double_t conjg_complexe_double(const complexe_double_t c);

int equal_complexe_float(const complexe_float_t c1, const complexe_float_t c2);

int equal_complexe_double(const complexe_double_t z1, const complexe_double_t z2);

#endif