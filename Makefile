stencil: stencil.c
	gcc -Ofast -pg -std=c99 -fopenmp-simd -Wall $^ -o $@

