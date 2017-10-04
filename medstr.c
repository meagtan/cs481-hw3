#include "medstr.h"
#include <stdio.h> // testing

// initial maximum size for dynamic arrays
#define MAXSIZ 10

const char *nucs = "ACTG";

BITSEQ *writebits(char *seq, size_t *size)
{
	size_t last = 0, maxsiz = MAXSIZ;
	BITSEQ *res = calloc(sizeof(BITSEQ), maxsiz);
	int i; // number of nucleotides in current BITSEQ, after one iteration of for loop

	for (i = 0; *seq; ++seq, ++i) {
		// if current BITSEQ full, create another one
		if (i == SEQSIZ) {
			i = 0;
			// extend array
			if (last + 1 == maxsiz) {
				maxsiz <<= 1;
				res = realloc(res, sizeof(BITSEQ)*maxsiz);
			}
			res[++last] = BITS(*seq);
		} else {
			// push current char onto topmost element of res
			res[last] = (res[last] << 2) | BITS(*seq);
			// printf("%llx\n", res[0]);
		}
	}

	// if last element of res not complete, left justify it
	res[last] <<= 2 * (SEQSIZ - i); // verify this

	*size = last * SEQSIZ + i;
	return res;
}

// branch-and-bound algorithm
BITSEQ medstr_char(char *dna[], size_t t, size_t k)
{
	BITSEQ str = 0;
	size_t i = 0; // depth in search tree

	// current distance, best distance, best string, etc.

	return str;
}
