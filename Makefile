main: main.c medstr.c
	gcc -O3 -m64 -march=corei7 -msse4.2 -o main main.c medstr.c	
