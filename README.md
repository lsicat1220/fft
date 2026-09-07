# The FFT

This repository contains the FFT code from the (Work in progress at the time of writing) audio visualizer/instrument tuner project.

In the `fft.c` file, you can find the various iterations of the FFT algorithm that I created based on my research into the fast Fourier transform.

- `dft` is self explanatory. It's used as a baseline for performance as well as a reference benchmark for numerical accuracy.
- `naive_fft` is the first successful iteration of the algorithm, using the basic recursive solution. 
Repeated twiddle calculations, dynamically allocated memory for the output.
- `recursive_fft` is my first exploration into optimizing the algorithm, and uses a static variable to store twiddle factors.
Clearly, this also runs into memory leakage.
- `fft_notwiddle` is similar to my final iteration of the FFT. Here, the twiddles are calculated before the FFT is performed. 
Recursion is also no longer used, and instead we use bit reversal to hold values in one array, which is allocated before calling the function.
- `fft`, the final iteration, simply asks the caller to provide its own array for the twiddles.

