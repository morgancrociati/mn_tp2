#include "complexe.h"

complexe_float_t add_complexe_float(const complexe_float_t c1, const complexe_float_t c2)
{
	complexe_float_t r;

	r.real = c1.real + c2.real;
	r.imaginary = c1.imaginary + c2.imaginary;

	return r;
}

complexe_double_t add_complexe_double(const complexe_double_t c1, const complexe_double_t c2)
{
	complexe_double_t r;

	r.real = c1.real + c2.real;
	r.imaginary = c1.imaginary + c2.imaginary;

	return r;
}

complexe_float_t mult_complexe_float(const complexe_float_t c1, const complexe_float_t c2)
{
	complexe_float_t r;

	r.real = c1.real * c2.real - c1.imaginary * c2.imaginary;
	r.imaginary = c1.real * c2.imaginary + c1.imaginary * c2.real;

	return r;
}

complexe_double_t mult_complexe_double(const complexe_double_t c1, const complexe_double_t c2)
{
	complexe_double_t r;

	r.real = c1.real * c2.real - c1.imaginary * c2.imaginary;
	r.imaginary = c1.real * c2.imaginary + c1.imaginary * c2.real;

	return r;
}

complexe_float_t conjg_complexe_float(const complexe_float_t c)
{
	complexe_float_t r;

	r.real = c.real;
	r.imaginary = -c.imaginary;

	return r;
}

complexe_double_t conjg_complexe_double(const complexe_double_t c)
{
	complexe_double_t r;

	r.real = c.real;
	r.imaginary = -c.imaginary;

	return r;
}

int equal_complexe_float(const complexe_float_t c1, const complexe_float_t c2)
{
	return (c1.real == c2.real) && (c1.imaginary == c2.imaginary);
}

int equal_complexe_double(const complexe_double_t z1, const complexe_double_t z2)
{
	return (z1.real == z2.real) && (z1.imaginary == z2.imaginary);
}
