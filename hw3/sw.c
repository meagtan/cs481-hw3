/*
 * CS 481
 * Ata Deniz Aydin
 * 21502637
 *
 * Implementation of Smith-Waterman algorithm as given in sw.h.
 */

#include <stdlib.h>

#include "sw.h"

// substitution matrix
const int subst[4][4] = {{ 4, -3, -2, -1},  // A
                         {-3,  4, -1, -2},  // C
                         {-2, -1,  4, -3},  // T
                         {-1, -2, -3,  4}}; // G
//                         A   C   T   G

const int gap = -4, gapop = -16, gapex = -4;

#define SETMAX(arr, rval, previ, prevj) do {    \
	if (arr[i+1][j+1] < (rval)) {   \
		arr[i+1][j+1] = (rval); \
		bestis[i+1][j+1] = previ; \
		bestjs[i+1][j+1] = prevj; \
	}} while (0)

void naivegap(FILE *out, char *seq1, char *seq2, int n, int m)
{
	// scoring and traceback matrix
	// the parent of (i,j) is (bestis[i][j], bestjs[i][j])
	int **scores, **bestis, **bestjs, i, j;
	int besti = 0, bestj = 0; // position of highest score

	// allocate matrices
	scores = calloc(n+1, sizeof(int *));
	bestis = calloc(n+1, sizeof(int *));
	bestjs = calloc(n+1, sizeof(int *));
	for (i = 0; i < n+1; ++i) {
		scores[i] = calloc(m+1, sizeof(int));
		bestis[i] = calloc(m+1, sizeof(int));
		bestjs[i] = calloc(m+1, sizeof(int));
	}

	// fill scoring matrix row by row
	for (i = 0; i < n; ++i) {
		for (j = 0; j < m; ++j) {
			// find best parent for scores[i+1][j+1], initially zero
			SETMAX(scores, scores[i][j] + S(seq1[i], seq2[j]), i, j);
			SETMAX(scores, scores[i+1][j] + gap, i+1, j);
			SETMAX(scores, scores[i][j+1] + gap, i, j+1);

			// update maximum score
			if (scores[i+1][j+1] > scores[besti][bestj]) {
				besti = i+1;
				bestj = j+1;
			}
		}
	}

	// traceback
	traceback(out, seq1, seq2, besti, bestj, scores, bestis, bestjs);

	// free matrices
	for (i = 0; i < n+1; ++i) {
		free(scores[i]);
		free(bestis[i]);
		free(bestjs[i]);
	}
	free(scores);
	free(bestis);
	free(bestjs);
}

void affinegap(FILE *out, char *seq1, char *seq2, int n, int m)
{
	// scoring and traceback matrices
	// the parent of (i,j) is (bestis[i][j], bestjs[i][j])
	int **v, **e, **f, **g, **bestis, **bestjs, i, j;
	int besti = 0, bestj = 0; // position of highest score

	// allocate matrices
	v = calloc(n+1, sizeof(int *));
	e = calloc(n+1, sizeof(int *));
	f = calloc(n+1, sizeof(int *));
	g = calloc(n+1, sizeof(int *));
	bestis = calloc(n+1, sizeof(int *));
	bestjs = calloc(n+1, sizeof(int *));
	for (i = 0; i < n+1; ++i) {
		v[i] = calloc(m+1, sizeof(int));
		e[i] = calloc(m+1, sizeof(int));
		f[i] = calloc(m+1, sizeof(int));
		g[i] = calloc(m+1, sizeof(int));
		bestis[i] = calloc(m+1, sizeof(int));
		bestjs[i] = calloc(m+1, sizeof(int));
	}

	// fill scoring matrices row by row
	// no need to initialize first row of e and f, will be 0 anyway
	for (i = 0; i < n; ++i) {
		for (j = 0; j < m; ++j) {
			g[i+1][j+1] = v[i][j] + S(seq1[i], seq2[j]);

			// maximize e[i+1][j+1]
			e[i+1][j+1] = e[i+1][j] + gapex;
			SETMAX(e, g[i+1][j] + gapop + gapex, 0, 0);
			SETMAX(e, f[i+1][j] + gapop + gapex, 0, 0);

			// maximize f[i+1][j+1]
			f[i+1][j+1] = f[i][j+1] + gapex;
			SETMAX(f, g[i][j+1] + gapop + gapex, 0, 0);
			SETMAX(f, e[i][j+1] + gapop + gapex, 0, 0);

			// maximize v[i+1][j+1], set bestis, bestjs
			SETMAX(v, e[i+1][j+1], i+1, j); // j+1 matches a space
			SETMAX(v, f[i+1][j+1], i, j+1); // i+1 matches a space
			SETMAX(v, g[i+1][j+1], i, j);   // i+1, j+1 match

			// update maximum score
			if (v[besti][bestj] < v[i+1][j+1]) {
				besti = i+1;
				bestj = j+1;
			}
		}
	}

	// traceback
	traceback(out, seq1, seq2, besti, bestj, v, bestis, bestjs);

	// free matrices
	for (i = 0; i < n+1; ++i) {
		free(v[i]);
		free(e[i]);
		free(f[i]);
		free(g[i]);
		free(bestis[i]);
		free(bestjs[i]);
	}
	free(v);
	free(e);
	free(f);
	free(g);
	free(bestis);
	free(bestjs);
}

void traceback(FILE *out, char *seq1, char *seq2, int besti, int bestj, int **scores, int **bestis, int **bestjs)
{
	int i, j;

	// write alignment backwards to arrays
	char *align1, *align2;
	int temp, len = 0; // length of alignment
	align1 = calloc(besti+bestj, sizeof(char)); // alignment can be of length at most n+m, corresponding to all indels
	align2 = calloc(besti+bestj, sizeof(char));

	// start from (besti, bestj), trace back until score 0
	for (i = besti, j = bestj; scores[i][j]; ++len) {
		align1[len] = (bestis[i][j] == i) ? '-' : seq1[i-1]; // if i does not change from previous, j matches a space
		align2[len] = (bestjs[i][j] == j) ? '-' : seq2[j-1]; // if j does not change from previous, i matches a space

		// (i, j) = (bestis[i][j], bestjs[i][j])
		temp = bestis[i][j];
		j    = bestjs[i][j];
		i    = temp;
	}

	// write arrays backwards to file
	fprintf(out, "%d\n", scores[besti][bestj]);
	for (i = len-1; i >= 0; --i)
		fputc(align1[i], out);
	fputc('\n', out);
	for (j = len-1; j >= 0; --j)
		fputc(align2[j], out);
	// fputc('\n', out);

	// free arrays
	free(align1);
	free(align2);
}
