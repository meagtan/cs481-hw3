/*
 * CS 481
 * Ata Deniz Aydin
 * 21502637
 *
 * Main function, reads input sequences and computes optimal naive and affine-gap alignments.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "sw.h"

#define INPUT "sequences.fasta"
#define OUT1  "naiveGap.txt"
#define OUT2  "affineGap.txt"

#define MAXLEN 1024

int main()
{
	// two sequences and their lengths
	char *seq1 = malloc(MAXLEN), // sizeof(char) == 1
	     *seq2 = malloc(MAXLEN);
	int n = 0, m = 0;

	// read input sequences
	FILE *f = fopen(INPUT, "r");
	if (!f) {
		fprintf(stderr, "Error: could not read file %s.\n", INPUT);
		return 1;
	}

	fscanf(f, ">Sequence1\n");
	// read until new sequence starts
	do {
		seq1[n] = fgetc(f);         // write current character to string
		if (isalpha(seq1[n])) ++n;  // if it is whitespace, skip it and override it with next character
	} while (n < MAXLEN && seq1[n] != '>');
	seq1[n] = '\0';

	fscanf(f, "Sequence2\n"); // > already read
	// read second sequence until eof
	do {
		seq2[m] = fgetc(f);
		if (isalpha(seq2[m])) ++m;
	} while (m < MAXLEN && !feof(f));
	seq2[m] = '\0';

	fclose(f);

	// for testing
	// printf("%d\t%s\n%d\t%s\n", n, seq1, m, seq2);

	// calculate alignments, write them to file
	f = fopen(OUT1, "w");
	if (!f) {
		fprintf(stderr, "Error: could not open file %s.\n", OUT1);
		return 1;
	}

	naivegap(f, seq1, seq2, n, m);
	fclose(f);

	/*
	f = fopen(OUT2, "w");
	if (!f) {
		fprintf(stderr, "Error: could not open file %s.\n", OUT2);
		return 1;
	}

	affinegap(f, seq1, seq2, n, m);
	fclose(f);
	*/

	free(seq1);
	free(seq2);
}
