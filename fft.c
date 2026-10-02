#include "complex.h"
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

Complex* dft(double* arr, int N, Complex* output) {
	for (int k = 0; k < N; k++) {
		// printf("k = %d\n", k);
		output[k].real = 0;
		output[k].imag = 0;
		for (int n = 0; n < N; n++) {
			//printf("\tn = %d\n", n);
			Complex term = {arr[n], 0};
			Complex omega = complexExp(-2 * M_PI * ((double)k / (double)N) * n);
			Complex product = complexMult(term, omega);
			output[k] = complexAdd(output[k], product);
			// printf("\t\toutput[k] = %lf + j%lf\n", output[k].real, output[k].imag);
		}
	}
	return output;
}

Complex* naive_fft(Complex* arr, int N, int stride) {
	Complex* output = malloc(sizeof(Complex) * N);

	if ((N & (N - 1)) != 0)	{
		fputs("ERROR: Must be a power of 2\n", stderr);
		return NULL;
	} else if (N == 1) {
		output[0] = arr[0];
		return output;
	}

	Complex* evens = naive_fft(arr, N/2, stride * 2);
	Complex* odds = naive_fft(arr + stride, N/2, stride * 2);
	
	for(int k = 0; k < N/2; k++) {
		Complex omega = complexExp(-2 * M_PI * ((double)k / N));
		Complex product = complexMult(odds[k], omega);
		output[k] = complexAdd(evens[k], product);
		product.real *= -1;
		product.imag *= -1;
		output[k + N/2] = complexAdd(evens[k], product);
	}

	free(evens);
	free(odds);
	return output;
}

int whatPowerOf2 (int N) {
	int output = 0;
	if ((N & (N - 1)) != 0)	{
		fputs("ERROR: Must be a power of 2\n", stderr);
		return -1;
	}
	for (; N > 1; N>>=1) {
		output++;
	}
	return output;
}

Complex* recursive_fft(Complex* arr, int N, int stride) {
	Complex* output = malloc(sizeof(Complex) * N);
	if (N == 1) {
		output[0] = arr[0];
		return output;
	}

	int power = whatPowerOf2(N);
	if (power == -1) {
		return NULL;
	}
	static Complex* factors = NULL;	
	if (factors == NULL) {
		factors = malloc(sizeof(Complex) * (power + 1));
		memset(factors, 0, sizeof(Complex) * (power + 1));
	} 
	if (factors[power].real + factors[power].imag == 0) {
		factors[power] = complexExp(-2 * M_PI / (double)N);	
	}

	Complex* evens = recursive_fft(arr, N/2, stride * 2);
	Complex* odds = recursive_fft(arr + stride, N/2, stride * 2);
	Complex omega = {1.0, 0};	
	for(int k = 0; k < N/2; k++) {
		Complex rhs;
		if (k == 0) {
			rhs = odds[k];
		} else {
			omega = complexMult(omega, factors[power]);	
			rhs = complexMult(odds[k], omega);
		}
		output[k] = complexAdd(evens[k], rhs);
		rhs.real *= -1;
		rhs.imag *= -1;
		output[k + N/2] = complexAdd(evens[k], rhs);
	}

	free(evens);
	free(odds);
	return output;
}

int reverseBits(int power, unsigned int input) {
	int output = 0;
	for (int i = 0; i < power; i++) {
		output = output<<1 | (input & 1);
		input >>= 1;
	}
	return output;
}

Complex* fft_notwiddle(double* arr, int N, Complex* output) {
	int power = whatPowerOf2(N);
	if (power == -1) {
		return NULL;
	}
	int stored_num = (N / 2) - 1;
	Complex* twiddles = malloc(sizeof(Complex) * ((N/2) - 1));
	for (int k = 1; k <= stored_num; k++) {
		twiddles[k - 1] = complexExp(-2 * M_PI * (double)k / N);
	}
	for (int i = 0; i < N; i++) {
		output[i].real = arr[reverseBits(power, i)];
		output[i].imag = 0;
	}	
	for (int stride = 1; stride < N; stride *= 2) {
		for (int k = 0; k < stride; k++) {
			for (int i = k; i < N; i += stride * 2) {
				Complex rhs;
				if (k != 0) {
					int index = (k * (N / (stride * 2))) - 1;
					rhs = complexMult(output[i + stride], twiddles[index]);
				} else {
					rhs = output[i + stride];
				}
				Complex sum1 = complexAdd(output[i], rhs);
				rhs.real *= -1;
				rhs.imag *= -1;
				Complex sum2 = complexAdd(output[i], rhs);
				output[i] = sum1;
				output[i + stride] = sum2;
			}	
		}
	}
	free(twiddles);
	return output;
}

Complex* fft(double* arr, int N, Complex* output, Complex* twiddles) {
	int power = whatPowerOf2(N);
	if (power == -1) {
		return NULL;
	}
	int stored_num = (N / 2) - 1;
	for (int k = 1; k <= stored_num; k++) {
		twiddles[k - 1] = complexExp(-2 * M_PI * (double)k / N);
	}
	for (int i = 0; i < N; i++) {
		output[i].real = arr[reverseBits(power, i)];
		output[i].imag = 0;
	}	
	for (int stride = 1; stride < N; stride *= 2) {
		for (int k = 0; k < stride; k++) {
			for (int i = k; i < N; i += stride * 2) {
				Complex rhs;
				if (k != 0) {
					int index = (k * (N / (stride * 2))) - 1;
					rhs = complexMult(output[i + stride], twiddles[index]);
				} else {
					rhs = output[i + stride];
				}
				Complex sum1 = complexAdd(output[i], rhs);
				rhs.real *= -1;
				rhs.imag *= -1;
				Complex sum2 = complexAdd(output[i], rhs);
				output[i] = sum1;
				output[i + stride] = sum2;
			}	
		}
	}
	return output;
}



