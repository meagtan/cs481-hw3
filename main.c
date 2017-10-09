#include <stdio.h>
#include "medstr.h"
#include <ctype.h>

#define INPUT  "input.txt"
#define OUTPUT "output.txt"

#define MAXSIZ 10
#define MAXCHAR 100

int main(void)
{

	// read input file, retrieve k and store sequences in dynamic array
	FILE *fin = fopen(INPUT, "r");
	if (!fin) {
		fprintf(stderr, "Error: input file %s does not exist.\n", INPUT);
		return -1;
	}


	// char **dna = malloc(MAXSIZ * sizeof(char *));
	size_t t = 0, maxsiz = MAXSIZ;

        BITSEQ **dnab = malloc(MAXSIZ * sizeof(BITSEQ *));
        size_t n;
	int k;

	fscanf(fin, "%d", &k);
	if (k >= SEQSIZ) {
		fprintf(stderr, "Error: only patterns of length < %d are allowed. We don't have 600 years to wait for the program to finish.\n", SEQSIZ);
		return -1;
	}
	fgetc(fin); // newline

	// won't be using this in final version, writebits should take stream
	while (!feof(fin)) {
		// if dynamic array full
		if (t == maxsiz) {
			maxsiz <<= 1;
			dnab = realloc(dnab, maxsiz * sizeof(BITSEQ *)); // dnab
		}

		// dna[t] = malloc(MAXCHAR * sizeof(char));
		// fgets(dna[t], MAXCHAR, fin); // may add newlines to strings

		dnab[t] = writebitsf(fin, &n);
		if (dnab[t] != NULL) ++t;

		// print sequence
		// printf("%d\t", t);
                // for (int i = 0; i < n; ++i)
                //         printf("%c", nucs[(dnab[t-1][i / SEQSIZ] >> (2 * SEQSIZ - 2 * (i % SEQSIZ) - 2)) & 3]); // CHAR(dnab[t-1])
                // printf("\n");

	}

	fclose(fin);

	// printf("File read\n");

	/*
	// create bit array
	BITSEQ **dnab = malloc(t * sizeof(BITSEQ *));
	size_t n;
	FILE *fin = fopen(INPUT, "r");
	// for (char c = fgetc(fin); c != '\n'; c = fgetc(fin));

	for (size_t i = 0; i < t; ++i) {
		dnab[i] = writebitsf(fin, &n); // later change to read directly from stream
	}

	fclose(fin);

	for (size_t j = 0; j < t; ++j) {
		for (int i = 0; i < n; ++i)
                	printf("%c", nucs[(dnab[j][i / SEQSIZ] >> (2 * SEQSIZ - 2 * (i % SEQSIZ) - 2)) & 3]);
		printf("\n");
	}
	*/

	// call median string function
	BITSEQ pat =  medstr_bit(dnab, t, n, k);  // medstr_char(dna, t, k);

	// size_t *pos = malloc(t * sizeof(size_t));
	// printf("Distance: %d\n", /* totdist_char(dna, t, pat, k, pos), */ totdist_bit(dnab, t, n, pat, k, pos));
	// free(pos);

	// write output value
	FILE *fout = fopen(OUTPUT, "w");
	if (!fout) {
		fprintf(stderr, "Error creating file %s.\n", OUTPUT);
		return -1;
	}

	// printf("Pattern: ");
	for (size_t i = 0; i < k; ++i) {
		fputc(CHAR(pat, i, k), fout);
		// putchar(CHAR(pat, i, k));
	}
	// printf("\n");
	fclose(fout);


	// free arrays
	for (size_t i = 0; i < t; ++i) {
		// free(dna[i]);
		free(dnab[i]); // comment out this line if not using dnab
	}
	// free(dna);
	free(dnab);

	return 0;
}
