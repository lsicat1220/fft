#ifndef COMPLEX_H
#define COMPLEX_H

#include <math.h>
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832
#endif

typedef struct {
	double real;
	double imag;
} Complex;

double magnitude(Complex num);

Complex complexExp(double num);

Complex complexAdd(Complex num1, Complex num2);

Complex complexMult(Complex num1, Complex num2);

Complex* doubleToComplexArr(double* arr, int N);

#endif
