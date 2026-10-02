#ifndef FFT_H
#define FFT_H

#include "complex.h"

Complex* recursive_fft(Complex* arr, int N, int stride);

Complex* naive_fft(Complex* arr, int N, int stride);

Complex* dft(double* arr, int N, Complex* output);
 
Complex* fft_nocache(double* arr, int N, Complex* output);

Complex* fft(double* arr, int N, Complex* output, Complex* twiddles);

Complex* fft_complex(Complex* arr, int N, Complex* output, Complex* twiddles);

int whatPowerOf2 (int N);

#endif
