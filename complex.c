#include "complex.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

double magnitude(Complex num) {
	return sqrt((num.real)*(num.real) + (num.imag)*(num.imag));
} 

Complex complexExp(double num) {
	Complex output = {cos(num), sin(num)};
	return output;
}

Complex complexAdd(Complex num1, Complex num2) {
	Complex output = {num1.real + num2.real, num1.imag + num2.imag};
	return output;
}

Complex complexMult(Complex num1, Complex num2) {
	double output_real = num1.real * num2.real - num1.imag * num2.imag;
	double output_imag = num1.real * num2.imag + num1.imag * num2.real;
	Complex output = {output_real, output_imag};
	return output;
}

Complex* doubleToComplexArr(double* arr, int N) {
	Complex* output = malloc(sizeof(Complex) * N);
	memset(output, 0, sizeof(Complex) * N);
	for (int i = 0; i < N; i++) {
		output[i].real = arr[i];	
	}
	return output;
}
