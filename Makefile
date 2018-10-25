stencil: stencil.c
	icc -xHost -Ofast -pg -std=c99 -Wall $^ -o $@

