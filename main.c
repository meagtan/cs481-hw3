#include <stdio.h>
#include "medstr.h"
#include <ctype.h>

#define INPUT  "input.txt"
#define OUTPUT "output.txt"

#define MAXSIZ 10
#define MAXCHAR 100

// extern int *dists;

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
	fgetc(fin); // skip newline

	// read each line from stream to dnab
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
	}

	fclose(fin);

	// call median string function
	BITSEQ pat = medstr_bit(dnab, t, n, k);  // */ medstr_char(dna, t, k);

	// report distance to stdout
	// size_t *pos = malloc(t * sizeof(size_t));
	// dists = malloc(t * sizeof(int));
	// printf("Distance: %d\n", /* totdist_char(dna, t, pat, k, pos), */ totdist_bit(dnab, t, n, pat, k, pos));
	// free(pos);
	// free(dists);

	// write output value
	FILE *fout = fopen(OUTPUT, "w");
	if (!fout) {
		fprintf(stderr, "Error creating file %s.\n", OUTPUT);
		return -1;
	}

	printbits(fout, pat, k);
	fputc('\n', fout);
	fclose(fout);


	// free arrays
	for (size_t i = 0; i < t; ++i) {
		// free(dna[i]);
		free(dnab[i]);
	}
	// free(dna);
	free(dnab);

	return 0;
}
