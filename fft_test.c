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
	int prev_failures = 0;

	Complex fft_output[MAX_N];
	Complex dft_output[MAX_N];
	Complex twiddles[MAX_N];
	

	printf("Beginning verification...\n");
	printf("[Test 1: Impulse]\n\n");
	
	y[0] = 1;
	for (int i = 1; i < MAX_N; i++) {
		y[i] = 0;	
	}
	for (int i = 0; i <= iterations; i++) {
		int unexpected_count = 0;
		int n = 1 << i;
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
			printf("Failure at n = %d\n", n);
			printf("Unexpected count: %d\n", unexpected_count);
			printf("Discrepancy count: %d\n", discrepancy_count);
			failures++;
		}
	}
	if (failures == prev_failures) {
		printf("Passed\n\n");
	} 	
	printf("[Test 2: DC sigmal]\n\n");
	prev_failures = failures;
	for (int i = 0; i < MAX_N; i++) {
		y[i] = 1;
	}
	for (int i = 0; i <= iterations; i++) {
		int unexpected_count = 0;
		int n = 1 << i;
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
			printf("Failure at n = %d\n", n);
			printf("Unexpected count: %d\n", unexpected_count);
			printf("Discrepancy count: %d\n", discrepancy_count);
			failures++;
		}
	}
	if (failures == prev_failures) {
		printf("Passed\n\n");
	} 	

	printf("[Test 3: No signal]\n\n");
	
	prev_failures = failures;
	for (int i = 0; i < MAX_N; i++) {
		y[i] = 0;
	} 
	for (int i = 0; i <= iterations; i++) {
		int unexpected_count = 0;
		int n = 1 << i;
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
			printf("Failure at n = %d\n", n);
			printf("Unexpected count: %d\n", unexpected_count);
			printf("Discrepancy count: %d\n", discrepancy_count);
			failures++;
		}
	}
	if (failures == prev_failures) {
		printf("Passed\n\n");
	} 	


	printf("[Test 4: Alternating -1 and 1]\n\n");
	prev_failures = failures;
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
			printf("Failure at n = %d\n", n);
			printf("Unexpected count: %d\n", unexpected_count);
			printf("Discrepancy count: %d\n", discrepancy_count);
			failures++;
		}
	}
	if (failures == prev_failures) {
		printf("Passed\n\n");
	} 	
	printf("[Test 5: Cosine wave]\n\n");
	prev_failures = failures;
	const int k = 3;
	for (int i = 2; i <= iterations; i++) {
		int unexpected_count = 0;
		int n = 1 << i;
		for (int j = 0; j < n; j++) {
			y[j] = cos(2 * M_PI * ((double) k / (double) n) * j);
		}
		dft(y, n, dft_output);
		fft(y, n, fft_output, twiddles);
		for (int j = 0; j < n; j++) {
			if (j == k || j == n - k) {
				if (fabs(fft_output[j].real - (double)n/2) > 1e-7 || fabs(fft_output[j].imag) > 1e-7) {
					unexpected_count++;
					printf("At index %d: Expected %d, got %lf\n", j, n/2, fft_output[j].real);
				}
			} else {
				if (fabs(fft_output[j].real) > 1e-7 || fabs(fft_output[j].imag) > 1e-7) {
					unexpected_count++;
					printf("At index %d: Expected 0, got %lf\n", j, fft_output[j].real);
				}
			}
		}
		int discrepancy_count = compareComplexArr(dft_output, fft_output, n);
		if (unexpected_count + discrepancy_count != 0) {
			printf("Failure at n = %d\n", n);
			printf("Unexpected count: %d\n", unexpected_count);
			printf("Discrepancy count: %d\n", discrepancy_count);
			failures++;
		}
	}
	if (failures == prev_failures) {
		printf("Passed\n\n");
	} 	
	printf("[Test 6: Random data & Reverse FFT]\n\n");
	prev_failures = failures;
	const int seed = 676767;
	srand(seed);
	Complex reverse[MAX_N];
	Complex y_complex[MAX_N];
	for (int set = 0; set < 10; set++) {
		prev_failures = failures;
		printf("Set %d\n", set);
		for (int i = 0; i < MAX_N; i++) {
			y[i] = ((double)rand() / RAND_MAX) * 2 - 1;
			//bounded -1 to 1
		}	
		for (int i = 0; i <= iterations; i++) {
			int n = 1<<i;
			dft(y, n, dft_output);
			fft(y, n, fft_output, twiddles);
			int discrepancy_count = compareComplexArr(dft_output, fft_output, n);
			if (discrepancy_count) {
				printf("Unexpected count: %d\n", discrepancy_count);
				failures++;
			}
			for (int j = 0; j < n; j++) {
				fft_output[j].imag *= -1;
			}	
			fft_complex(fft_output, n, reverse, twiddles);
			for (int j = 0; j < n; j++) {
				reverse[j].real /= n;
				reverse[j].imag /= n;
				reverse[j].imag *= -1;
				y_complex[j].real = y[j];
				y_complex[j].imag = 0;
			}	
			discrepancy_count = compareComplexArr(y_complex, reverse, n);
			if (discrepancy_count) {
				printf("Failure at n = %d\n", n);
				printf("Unexpected count for reverse FFT: %d\n", discrepancy_count);
				failures++;
			}	
		}	
		if (failures == prev_failures) {
			printf("Passed\n\n");
		} 	
	}	
	prev_failures = failures;
	printf("[Test 7: Invalid data]\n\n");
	if (fft(y, 3, fft_output, twiddles) != NULL) {
		printf("Incorrectly outputted non-null for a non-power of 2 n\n");
		failures++;
	}
	if (fft(y, -1, fft_output, twiddles) != NULL) {
		printf("Incorrectly outputted non-null for n = -1\n");
		failures++;
	}
	if (fft(y, 0, fft_output, twiddles) != NULL) {
		printf("Incorrectly outputted non-null for n = 0\n");
		failures++;
	}
	if (failures == prev_failures) {
		printf("Passed\n\n");
	}
	if (failures) {
		return EXIT_FAILURE;
	}
	return 0;
}
