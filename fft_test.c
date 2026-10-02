#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "complex.h"
#include "fft.h"
#include "assert.h"
#define N 256

#define record_time(function, output) \
	do { \
		clock_t start = clock(); \
		function; \
		clock_t end = clock(); \
		double exe_time = (double) (end - start) / CLOCKS_PER_SEC; \
		printf("Execution time: %lf\n", exe_time); \
		output = exe_time; \
	} while (0) \

double* linspace(int time, int n) {
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
	
	double y[N];
	int iterations = whatPowerOf2(N);

	Complex fft_output[N];
	Complex dft_output[N];
	Complex twiddles[N];

	printf("Beginning verification...\n");
	printf("Test 1: Impulse/DC Signal\n\n");
	
	y[0] = 1;
	for (int i = 1; i < N; i++) {
		y[i] = 0;	
	}
	int unexpected_count = 0;
	int discrepancy_count = 0;
	for (int i = 0; i < iterations; i++) {
		int n = 1 << i;
		printf("--- n = %d ---\n", n);
		dft(y, n, dft_output);
		fft(y, n, fft_output, twiddles);	
		for (int i = 0; i < n; i++) {
			if (fft_output[i].real != 1 || fft_output[i].imag != 0) {
				unexpected_count++;
			}
			double delta_real = fabs(fft_output[i].real - dft_output[i].real);
			double delta_imag = fabs(fft_output[i].imag - dft_output[i].imag);
			if (delta_real > 1e-7 + 1e-7 * dft_output[i].real) {
				printf("Significant real delta at i = %d: %lf\n", i, delta_real);
				discrepancy_count++;
			}
			if (delta_imag > 1e-7 + 1e-7 * dft_output[i].imag) {
				printf("Significant imaginary delta at i = %d: %lf\n", i, delta_imag);
				discrepancy_count++;
			}
		}
	}
	if (unexpected_count + discrepancy_count != 0) {
		return -1;
	}
	return 0;
}
