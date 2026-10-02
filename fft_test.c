#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "complex.h"
#include "fft.h"
#include "assert.h"
#define MAX_N 256

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
	double* output = malloc(sizeof(double) * MAX_N);
	double stride = (double)time / MAX_N;
	double val = 0;
	for (int i = 0; i < MAX_N; i++) {
		output[i] = val;
		val += stride;
	}
	return output;
}

int compareComplexArr(Complex* reference, Complex* arr, int n) {
	int discrepancy_count = 0;
	for (int i = 0; i < n; i++) {
		double delta_real = fabs(arr[i].real - reference[i].real);
		double delta_imag = fabs(arr[i].imag - reference[i].imag);
		if (!isfinite(delta_real) || !isfinite(delta_imag)) {
			discrepancy_count++;
		}
		if (delta_real > 1e-7 + 1e-7 * fabs(reference[i].real)) {
			printf("Significant real delta at i = %d: %lf\n", i, delta_real);
			discrepancy_count++;
		}
		if (delta_imag > 1e-7 + 1e-7 * fabs(reference[i].imag)) {
			printf("Significant imaginary delta at i = %d: %lf\n", i, delta_imag);
			discrepancy_count++;
		}
	}
	return discrepancy_count;
}

int main() {
	printf("Hello world!\n");
	
	double y[MAX_N];
	int iterations = whatPowerOf2(MAX_N);
	int failures = 0;

	Complex fft_output[MAX_N];
	Complex dft_output[MAX_N];
	Complex twiddles[MAX_N];
	

	printf("Beginning verification...\n");
	printf("Test 1: Impulse\n\n");
	
	y[0] = 1;
	for (int i = 1; i < MAX_N; i++) {
		y[i] = 0;	
	}
	for (int i = 0; i <= iterations; i++) {
		int unexpected_count = 0;
		int n = 1 << i;
		printf("--- n = %d ---\n", n);
		dft(y, n, dft_output);
		fft(y, n, fft_output, twiddles);	
		for (int i = 0; i < n; i++) {
			if (fabs(fft_output[i].real - 1) > 1e-7 || fabs(fft_output[i].imag) > 1e-7) {
				printf("At index %d: Expected 1, got %lf\n", i, fft_output[i].real);
				unexpected_count++;
			}
		}
		int discrepancy_count = compareComplexArr(dft_output, fft_output, n);
		if (unexpected_count + discrepancy_count != 0) {
			printf("Unexpected count: %d\n", unexpected_count);
			printf("Discrepancy count: %d\n", discrepancy_count);
			failures++;
		}
	}

	printf("Test 2: DC sigmal\n\n");
	for (int i = 0; i < MAX_N; i++) {
		y[i] = 1;
	}
	for (int i = 0; i <= iterations; i++) {
		int unexpected_count = 0;
		int n = 1 << i;
		printf("--- n = %d ---\n", n);
		dft(y, n, dft_output);
		fft(y, n, fft_output, twiddles);
		if (fabs(fft_output[0].real - n) > 1e-7 || fabs(fft_output[0].imag) > 1e-7) {
			unexpected_count++;
			printf("At index 0: Expected %d, got %lf\n", n, fft_output[0].real);
		}
		for (int j = 1; j < n; j++) {
			if (fabs(fft_output[j].real) > 1e-7 || fabs(fft_output[j].imag) > 1e-7) {
				unexpected_count++;
				printf("At index %d: Expected 0, got %lf\n", j, fft_output[j].real);
			}
		}
		int discrepancy_count = compareComplexArr(dft_output, fft_output, n);
		if (unexpected_count + discrepancy_count != 0) {
			printf("Unexpected count: %d\n", unexpected_count);
			printf("Discrepancy count: %d\n", discrepancy_count);
			failures++;
		}
	}

	printf("Test 3: No signal\n\n");
	for (int i = 0; i < MAX_N; i++) {
		y[i] = 0;
	} 
	for (int i = 0; i <= iterations; i++) {
		int unexpected_count = 0;
		int n = 1 << i;
		printf("--- n = %d ---\n", n);
		dft(y, n, dft_output);
		fft(y, n, fft_output, twiddles);
		for (int i = 0; i < n; i++) {
			if (fabs(fft_output[i].real) > 1e-7 || fabs(fft_output[i].imag) > 1e-7) {
				unexpected_count++;
				printf("At index %d: Expected 0, got %lf\n", i, fft_output[i].real);
			}
		}
		int discrepancy_count = compareComplexArr(dft_output, fft_output, n);
		if (unexpected_count + discrepancy_count != 0) {
			printf("Unexpected count: %d\n", unexpected_count);
			printf("Discrepancy count: %d\n", discrepancy_count);
			failures++;
		}
	}

	printf("Test 4: Alternating -1 and 1\n\n");
	for (int i = 0; i < MAX_N; i++) {
		if (i % 2) {
			y[i] = -1;
		} else {
			y[i] = 1;
		}
	}
	for (int i = 0; i <= iterations; i++) {
		int unexpected_count = 0;
		int n = 1 << i;
		printf("--- n = %d ---\n", n);
		dft(y, n, dft_output);
		fft(y, n, fft_output, twiddles);
		for (int i = 0; i < n; i++) {
			if (i == n / 2) {
				if (fabs(fft_output[i].real - n) > 1e-7 || fabs(fft_output[i].imag) > 1e-7) {
					unexpected_count++;
					printf("At index %d: Expected %d, got %lf\n", i, n, fft_output[i].real);
				}
			} else {
				if (fabs(fft_output[i].real) > 1e-7 || fabs(fft_output[i].imag) > 1e-7) {
					unexpected_count++;
					printf("At index %d: Expected 0, got %lf\n", i, fft_output[i].real);
				}
			}
		}
		int discrepancy_count = compareComplexArr(dft_output, fft_output, n);
		if (unexpected_count + discrepancy_count != 0) {
			printf("Unexpected count: %d\n", unexpected_count);
			printf("Discrepancy count: %d\n", discrepancy_count);
			failures++;
		}
	}
	printf("Test 5: Cosine wave\n\n");
	const int k = 3;
	for (int i = 2; i <= iterations; i++) {
		int unexpected_count = 0;
		int n = 1 << i;
		printf("--- n = %d ---\n", n);
		for (int j = 0; j < n; j++) {
			y[j] = cos(2 * M_PI * ((double) k / (double) n) * j);
		}
		dft(y, n, dft_output);
		fft(y, n, fft_output, twiddles);
		for (int i = 0; i < n; i++) {
			if (i == k || i == n - k) {
				if (fabs(fft_output[i].real - (double)n/2) > 1e-7 || fabs(fft_output[i].imag) > 1e-7) {
					unexpected_count++;
					printf("At index %d: Expected %d, got %lf\n", i, n/2, fft_output[i].real);
				}
			} else {
				if (fabs(fft_output[i].real) > 1e-7 || fabs(fft_output[i].imag) > 1e-7) {
					unexpected_count++;
					printf("At index %d: Expected 0, got %lf\n", i, fft_output[i].real);
				}
			}
		}
		int discrepancy_count = compareComplexArr(dft_output, fft_output, n);
		if (unexpected_count + discrepancy_count != 0) {
			printf("Unexpected count: %d\n", unexpected_count);
			printf("Discrepancy count: %d\n", discrepancy_count);
			failures++;
		}
	}
	if (failures) {
		return EXIT_FAILURE;
	}
	return 0;
}
