fft_test: fft.c complex.c fft_test.c
	gcc -Wall -lm fft.c complex.c fft_test.c -o ./bin/fft_test
