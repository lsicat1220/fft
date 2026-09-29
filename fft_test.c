#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "complex.h"
#include "fft.h"

#define record_time(function, output) \
	do { \
		clock_t start = clock(); \
		function; \
		clock_t end = clock(); \
		double exe_time = (double) (end - start) / CLOCKS_PER_SEC; \
		printf("Execution time: %lf\n", exe_time); \
		output = exe_time; \
	} while (0) \

double f(double x) {
	return sin(2 * M_PI * 15 * x) + sin(2 * M_PI * 2.3 * x);
}

double* linspace(int time, int N) {
	double* output = malloc(sizeof(double) * N);
	double stride = (double)time / N;
	double val = 0;
	for (int i = 0; i < N; i++) {
		output[i] = val;
		val += stride;
	}
	return output;
}

int main() {
	printf("Hello world!\n");
	int time = 4;
	int f_s = 512;
	int N = time * f_s;
	
	double* x = linspace(time, N);
	double* y = malloc(sizeof(double) * N);
	for (int i = 0; i < N; i++) {
		y[i] = f(x[i]);
	}
	Complex* new_output = malloc(sizeof(Complex) * N);
	Complex* dft_output = malloc(sizeof(Complex) * N);
	printf("DFT output:\n");

	double old_time;
	record_time({
		Complex* input = doubleToComplexArr(y, N);
		dft_output = dft(input, N, new_output);
	}, old_time);

	Complex* twiddles = malloc(sizeof(Complex) * N);
	
	printf("New output:\n");

	double new_time;
	record_time({
		new_output = fft(y, N, new_output, twiddles);
	}, new_time);

	printf("Time difference: %lf, %lf percent\n", old_time - new_time, ((new_time - old_time) / old_time) * 100);
	int equal = 1;
	for (int i = 0; i < N; i++) {
		if (magnitude(dft_output[i]) - magnitude(new_output[i]) > 0.000001) {
			equal = 0; 
			break;
		}
	}
	printf("Equal = %d\n", equal);
	return 0;
}
