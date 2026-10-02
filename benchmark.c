#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "complex.h"
#include "fft.h"

#define MAX_N 256
#define NUM_CALLS 50 
int main() {
	printf("Benchmark test\n\n");
	const int seed = 676767;
	srand(seed);
	struct timespec start, end;
	double dft_elapsed_time = 0;
	double fft_elapsed_time = 0;
	Complex fft_output[MAX_N];
	Complex dft_output[MAX_N];
	Complex twiddles[MAX_N / 2 - 1];
	double y[MAX_N];
	for (int i = 0; i < MAX_N; i++) {
		y[i] = ((double) rand() / RAND_MAX) * 2 - 1;
	}	
	dft(y, MAX_N, dft_output);
	fft(y, MAX_N, fft_output, twiddles);
	static volatile double sink;
	
	for (int power = 4; power <= 8; power++) {
		int n = 1<<power;
		double checksum = 0;
		printf("[n = %d]\n", n);
		clock_gettime(CLOCK_MONOTONIC, &start);
		for (int i = 0; i < NUM_CALLS; i++) {
			dft(y, n, dft_output);
			for (int j = 0; j < n; j++) {
				checksum += dft_output[j].real + dft_output[j].imag;
			}
			sink = checksum;
		}
		clock_gettime(CLOCK_MONOTONIC, &end);
		dft_elapsed_time = (end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) / 1e9;
		dft_elapsed_time /= NUM_CALLS;
		printf("DFT elapsed time per call: %lf\n", dft_elapsed_time);

		clock_gettime(CLOCK_MONOTONIC, &start);
		for (int i = 0; i < NUM_CALLS; i++) {
			fft(y, n, fft_output, twiddles);
			for (int j = 0; j < n; j++) {
				checksum += fft_output[j].real + fft_output[j].imag;
			}
			sink = checksum;
		}
		clock_gettime(CLOCK_MONOTONIC, &end);
		fft_elapsed_time = (end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec) / 1e9;
		fft_elapsed_time /= NUM_CALLS;
		printf("FFT elapsed time per call: %lf\n", fft_elapsed_time);
		printf("Speed up: %lf\n\n", dft_elapsed_time / fft_elapsed_time);
	}
}
