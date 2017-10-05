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

	int k;
	fscanf(fin, "%d", &k);
	if (k >= SEQSIZ) {
		fprintf(stderr, "Error: only patterns of length < %d are allowed. We don't have 600 years to wait for the program to finish.\n", SEQSIZ);
		return -1;
	}
	fgetc(fin); // newline

	char **dna = malloc(MAXSIZ * sizeof(char *));
	size_t t = 0, maxsiz = MAXSIZ;

        // BITSEQ **dnab = malloc(MAXSIZ * sizeof(BITSEQ *));
        // size_t n;

	// won't be using this in final version, writebits should take stream
	while (!feof(fin)) {
		// if dynamic array full
		if (t == maxsiz) {
			maxsiz <<= 1;
			dna = realloc(dna, maxsiz * sizeof(BITSEQ *)); // dnab
		}

		dna[t] = malloc(MAXCHAR * sizeof(char));
		fgets(dna[t], MAXCHAR, fin); // may add newlines to strings
		if (isalpha(dna[t][0])) ++t;
		else free(dna[t]);

		/*
		dnab[t] = writebitsf(fin, &n);
		if (dnab[t] != NULL) ++t;
		*/

		/* print sequence
		printf("%d\t", t);
                for (int i = 0; i < n; ++i)
                        printf("%c", nucs[(dnab[t-1][i / SEQSIZ] >> (2 * SEQSIZ - 2 * (i % SEQSIZ) - 2)) & 3]); // CHAR(dnab[t-1])
                printf("\n");
		*/
	}

	fclose(fin);

	// printf("File read\n");


	// create bit array
	BITSEQ **dnab = malloc(t * sizeof(BITSEQ *));
	size_t n;
	for (size_t i = 0; i < t; ++i) {
		dnab[i] = writebits(dna[i], &n); // later change to read directly from stream
	}

/*
	for (size_t j = 0; j < t; ++j) {
		for (int i = 0; i < n; ++i)
                	printf("%c", nucs[(dnab[j][i / SEQSIZ] >> (2 * SEQSIZ - 2 * (i % SEQSIZ) - 2)) & 3]);
		printf("\n");
	}
*/

	// call median string function
	BITSEQ pat = medstr_char(dna, t, k);

	// printf("%d\n", pat);

	// write output value
	FILE *fout = fopen(OUTPUT, "w");
	if (!fout) {
		fprintf(stderr, "Error creating file %s.\n", OUTPUT);
		return -1;
	}

	for (size_t i = 0; i < k; ++i) {
		fputc(CHAR(pat, i, k), fout);
		putchar(CHAR(pat, i, k));
	}
	// printf("\n%llx\n", pat);
	fclose(fout);


	// free arrays
	for (size_t i = 0; i < t; ++i) {
		free(dna[i]);
		free(dnab[i]); // comment out this line if not using dnab
	}
	free(dna);
	free(dnab);

	return 0;
}
