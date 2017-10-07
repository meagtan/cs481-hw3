#include "medstr.h"
#include <limits.h>    // may not need this
#include <ctype.h>
#include <x86intrin.h> // builtin popcnt, why use a for loop for something the CPU handles

// initial maximum size for dynamic arrays
#define MAXSIZ 10

#define BYPASS(str) do { ++str; while (!(str & 3)) {str >>= 2; --i;}} while (0)
#define NEXT(str)   do { if (i == k) BYPASS(str); else {str <<= 2; ++i;} } while (0)

#define EVENDIGITS 0x5555555555555555ll

const char *nucs = "ACTG";

BITSEQ *writebits(char *seq, size_t *size)
{
	size_t last = 0, maxsiz = MAXSIZ;
	BITSEQ *res = calloc(sizeof(BITSEQ), maxsiz);
	int i; // number of nucleotides in current BITSEQ, after one iteration of for loop

	for (i = 0; isalpha(*seq); ++seq, ++i) {
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

BITSEQ *writebitsf(FILE *f, size_t *size)
{
        size_t last = 0, maxsiz = MAXSIZ;
        BITSEQ *res = calloc(sizeof(BITSEQ), maxsiz);
        int i; // number of nucleotides in current BITSEQ, after one iteration of for loop
	char c;

	// skip non-alphabetic characters
	for (c = fgetc(f); c != EOF && !isalpha(c); c = fgetc(f));

	if (c == EOF) {
		free(res);
		return NULL;
	}

        for (i = 0; !feof(f) && isalpha(c); c = fgetc(f), ++i) {
                // if current BITSEQ full, create another one
                if (i == SEQSIZ) {
                        i = 0;
                        // extend array
                        if (last + 1 == maxsiz) {
                                maxsiz <<= 1;
                                res = realloc(res, sizeof(BITSEQ)*maxsiz);
                        }
                        res[++last] = BITS(c);
                } else {
                        // push current char onto topmost element of res
                        res[last] = (res[last] << 2) | BITS(c);
                        // printf("%llx\n", res[0]);
                }
        }

        // if last element of res not complete, left justify it
        res[last] <<= 2 * (SEQSIZ - i); // verify this

        *size = last * SEQSIZ + i;
        return res;
}

// branch-and-bound algorithm
/* Possible optimization strategy exploiting some optimal substructure:
 * Store the best positions of the last pattern examined, and try to bound the distance of the current pattern to each string.
 * - If str is reached via a bypass, omitting any backtracking, its last character was changed to another character.
 *   If, for any sequence seq, the best position j of last str on seq satisfies seq[j+i-1] == str[i-1] for current str,
 *   then we can say the distance of str to seq is d-1 with best position j, where d is the best distance of str to last seq. Otherwise run usual distance.
 * - If str is reached from a parent, its last character was added to the parent. If the best distance of seq to parent was d at position j, and
 *   seq[j+i-1] == str[i-1] for current str, we can extend the best position of the parent to the child, as any other occurrence of the child has to be
 *   at least d away from the parent, independently of the matching of the last character. Otherwise calculate distance as usual again.
 * Try this after you test _bits functions.
 */
BITSEQ medstr_char(char *dna[], size_t t, size_t k)
{
	BITSEQ str = 0, best = 0;
	size_t i = 1; // depth in search tree, < k
	int bestdist = t * k, dist;
	size_t *bestpos, *pos; // may use bestpos and pos later, not necessary now
	pos = calloc(t, sizeof *pos);

	// to move to child, shift str by 2 and increment i
	// to go to next leaf or bypass, increment str
	// if str & 3 == 3, backtrack by shifting str to the right and decrementing 1
	while (i) {
		/*
		printf("%d\t%llx\t", i, str);
		printbits(stdout, str, i);
		printf("%llx", str);
		*/

		dist = totdist_char(dna, t, str, i, pos);
		/*
		if (str == 0xa07)
			printf("\t%d\t%d", mindist_char(dna[0], str, k, pos), dist_char(dna[0] + 24, str, k));
		printf("\t%d\n", dist);
		*/
		if (i < k) {
			if (dist > bestdist)
				BYPASS(str);
			else
				NEXT(str);
		} else {
			if (dist < bestdist) {
				bestdist = dist;
				best = str;
				// printf("%llx\t%d\n", best, dist);
				// bestpos = pos;
			}
			NEXT(str);
		}// */
		// NEXT(str);
	}

	return best;
}

BITSEQ medstr_bit(BITSEQ *dna[], size_t t, size_t n, size_t k)
{
	BITSEQ str = 0, best = 0;
	size_t i = 1; // depth in search tree, < k
	int bestdist = t * k, dist;
	size_t *bestpos, *pos; // may use bestpos and pos later, not necessary now
	pos = calloc(t, sizeof *pos);

	// to move to child, shift str by 2 and increment i
 	// to go to next leaf or bypass, increment str
	// if str & 3 == 3, backtrack by shifting str to the right and decrementing 1
	while (i) {
		dist = totdist_bit(dna, t, n, str, i, pos);

		if (i < k) {
			if (dist > bestdist)
				BYPASS(str);
			else
				NEXT(str);
		} else {
			if (dist < bestdist) {
				bestdist = dist;
				best = str;
				// bestpos = pos; // do I need to keep this?
                        }
			NEXT(str);
		}
	}
	free(pos); // remove when not using pos

	return best;
}

// totdist

int totdist_char(char *dna[], size_t t, BITSEQ str, size_t k, size_t *pos)
{
	int res = 0;
	// iterate through each sequence in dna
	for (size_t i = 0; i < t; ++i)
		res += mindist_char(dna[i], str, k, pos+i);
	return res;
}

int totdist_bit(BITSEQ *dna[], size_t t, size_t n, BITSEQ str, size_t k, size_t *pos)
{
	int res = 0;
	// printf("totdist_bit(dna, %zu, %zu, %llx, %zu, pos)\n", t, n, str, k);
	for (size_t i = 0; i < t; ++i)
		res += mindist_bit(dna[i], n, str, k, pos+i);
	return res;
}

// mindist

// utilize pos here, also take the previous distance and whether it's a parent or sibling
int mindist_char(char *seq, BITSEQ str, size_t k, size_t *pos)
{
	int res = k, dist;
	for (size_t i = 0; seq[i+k-1]; ++i) { // stop if substring of length k contains null terminator
		dist = dist_char(seq+i, str, k);
		if (dist < res) {
			res = dist;
			// printf("*pos = %zu\n", i);
			*pos = i;
		}
	}
	return res;
}

int mindist_bit(BITSEQ *seq, size_t n, BITSEQ str, size_t k, size_t *pos)
{
	int res = k, dist;
	// i represents ending position, not included in string
	for (size_t i = k; i <= n; ++i) {
	        // printf("dist_bit(seq, %zu, %llx, %zu)\n", i, str, k);
		dist = dist_bit(seq, i, str, k);
		if (dist < res) {
			res = dist;
			*pos = i - k + 1;
		}
	}
	// printf("\n");
	return res;
}

// dist

// assumes seq has at least k characters
int dist_char(char *seq, BITSEQ str, size_t k)
{
	int res = 0;
	for (size_t i = 0; i < k && seq[i]; ++i)
		res += (seq[i] != CHAR(str, i, k));
	return res;
}

// test this
int dist_bit(BITSEQ *seq, size_t i, BITSEQ str, size_t k)
{
        // printf("dist_bit(seq, %zu, %llx, %zu)\n", i, str, k);
	seq += i / SEQSIZ;
	i %= SEQSIZ;

	// printf("dist_bit(seq, %zu, %llx, %zu)\n", i, str, k);
	// whether pattern contained in one BITSEQ
	if (i >= k) {
		// right justify and compare
		BITSEQ pat = (*seq >> 2*(SEQSIZ - i)) & (1ll << 2*k) - 1;
		// printbits(stdout, pat, k); printbits(stdout,str,k);
		// putchar('\n');
		// printf("%llx\t%llx\t", str, pat);
		str ^= pat;
	} else {
		// left shift end of *(seq-1) by i and add to right justified portion in *seq
		BITSEQ patl = (*(seq-1) << 2*i) & (1ll<<2*k)-1,
		       patr = (*seq >> 2*(SEQSIZ - i));
		// printbits(stdout,patl|patr,k); printbits(stdout,patr,i); printbits(stdout,str,k);
		// putchar('\n');
		// printf("%llx\t%llx\t%llx\t", str, patl, patr);
		str ^= (patl | patr) & (1ll << 2*k)-1;
	}

	// need popcount in base 4
	str = (str & EVENDIGITS) | ((str & (EVENDIGITS << 1)) >> 1);
	int res = __builtin_popcountll(str); // count nonzero bits in str
	// printf("%d\n", res);
	return res;
}

// print

void printbits(FILE *f, BITSEQ str, size_t k)
{
	for (int i = 0; i < k; ++i)
		fputc(nucs[(str >> (2 * (k - i - 1))) & 3], f); // fputc(CHAR(str, i), f);
	fputc('\t', f);
}
