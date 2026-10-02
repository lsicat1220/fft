fft_test: fft.c complex.c fft_test.c
	gcc -Wall -lm fft.c complex.c fft_test.c -o ./bin/fft_test

benchmark: fft.c complex.c benchmark.c
	gcc -lm -O3 fft.c complex.c benchmark.c -o ./bin/benchmark
