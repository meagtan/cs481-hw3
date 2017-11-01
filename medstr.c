#include "medstr.h"
#include <ctype.h>
#include <x86intrin.h> // builtin popcnt, why use a for loop for something the CPU handles

// initial maximum size for dynamic arrays
#define MAXSIZ 10

#define BYPASS(str) do { ++str; /* prevpat = NEIGHBOR; */ while (!(str & 3)) {str >>= 2; --i; /* prevpat = FIRST; */}} while (0)
#define NEXT(str)   do { if (i == k) BYPASS(str); else {str <<= 2; ++i; /* prevpat = PARENT; */} } while (0)

#define EVENDIGITS 0x5555555555555555ll

const char *nucs = "ACTG";

/* Possible optimization strategy exploiting some optimal substructure:
 * Store the best positions of the last pattern examined, and try to bound the distance of the current pattern to each string.
 * - If str is reached via a bypass, omitting any backtracking, its last character was changed to another character.
 *   If, for any sequence seq, the best position j of last str on seq satisfies seq[j+i-1] == str[i-1] for current str,
 *   then we can say the distance of str to seq is d-1 with best position j, where d is the best distance of str to last seq. Otherwise run usual distance.
 * - If str is reached from a parent, its last character was added to the parent. If the best distance of seq to parent was d at position j, and
 *   seq[j+i-1] == str[i-1] for current str, we can extend the best position of the parent to the child, as any other occurrence of the child has to be
 *   at least d away from the parent, independently of the matching of the last character. Otherwise calculate distance as usual again.
 * Or simply, if the first k-1 characters of the current pattern are the same as those of the previous pattern, and the kth character of the current pattern
 *  occurs on the best position of the previous pattern, then the best position of the current pattern will also be the best position of the previous pattern,
 *  with distance decremented by 1 if the previous pattern was also k characters.
 * After testing: strangely and disappointingly, this actually slows down computation by about 20%. Related code commented out, including the global variables below.
 */
// for convenience, will use global variables instead of changing function arguments
// true if the previous pattern shares its first k-1 characters with the current pattern
// enum {FIRST, PARENT, NEIGHBOR} prevpat;
// previous distance values calculated
// int *dists = NULL;


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
		}
	}

	// if last element of res not complete, left justify it
	res[last] <<= 2 * (SEQSIZ - i);

	*size = last * SEQSIZ + i;
	return res;
}

// reads one row
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
                }
        }

        // if last element of res not complete, left justify it
        res[last] <<= 2 * (SEQSIZ - i);

        *size = last * SEQSIZ + i;
        return res;
}

// branch-and-bound algorithm
BITSEQ medstr_char(char *dna[], size_t t, size_t k)
{
	BITSEQ str = 0, best = 0;
	size_t i = 1; // depth in search tree, < k
	int bestdist = t * k, dist;
	size_t *bestpos, *pos; // may use bestpos and pos later, not necessary now
	pos = calloc(t, sizeof *pos);
	// prevpat = 0;
	// dists = calloc(t, sizeof(int));

	// to move to child, shift str by 2 and increment i
	// to go to next leaf or bypass, increment str
	// if str & 3 == 3, backtrack by shifting str to the right and decrementing 1
	while (i) {
		dist = totdist_char(dna, t, str, i, pos);
		if (i < k) {
			if (dist > bestdist)
				BYPASS(str);
			else
				NEXT(str);
		} else {
			if (dist < bestdist) {
				bestdist = dist;
				best = str;
				// bestpos = pos;
			}
			NEXT(str);
		}
	}

	// free(dists);
	// dists = NULL;
	return best;
}

BITSEQ medstr_bit(BITSEQ *dna[], size_t t, size_t n, size_t k)
{
	BITSEQ str = 0, best = 0;
	size_t i = 1; // depth in search tree, < k
	int bestdist = t * k, dist;
	size_t *bestpos, *pos; // may use bestpos and pos later, not necessary now
	pos = calloc(t, sizeof *pos);
	// prevpat = 0;
	// dists = calloc(t, sizeof(int));

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
				// bestpos = pos;
                        }
			NEXT(str);
		}
	}
	free(pos); // remove when not using pos

	// free(dists);
	// dists = NULL;
	// prevpat = 0;
	return best;
}

// totdist

int totdist_char(char *dna[], size_t t, BITSEQ str, size_t k, size_t *pos)
{
	int res = 0;
	// iterate through each sequence in dna
	for (size_t i = 0; i < t; ++i) {
		// if (!prevpat || pos[i]+k-1 >= strlen(dna[i]) || dna[i][pos[i]+k-1] != (str & 3))
		// 	dists[i] = mindist_char(dna[i], str, k, pos+i);
		// else if (prevpat == NEIGHBOR)
		// 	--dists[i];
		// res += dists[i];
		res += mindist_char(dna[i], str, k, pos+i);
	}
	return res;
}

int totdist_bit(BITSEQ *dna[], size_t t, size_t n, BITSEQ str, size_t k, size_t *pos)
{
	int res = 0;
	for (size_t i = 0; i < t; ++i) {
		// if (!prevpat || pos[i]+k-1 >= n || CHAR(dna[i][(pos[i]+k-1)/SEQSIZ], (pos[i]+k-1)%SEQSIZ, SEQSIZ) != CHAR(str, k-1, k))
		// 	dists[i] = mindist_bit(dna[i], n, str, k, pos+i);
		// else if (prevpat == NEIGHBOR)
		// 	--dists[i]; // last character was mismatch, now match
		// res += dists[i];
		res += mindist_bit(dna[i], n, str, k, pos+i);
	}
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
		dist = dist_bit(seq, i, str, k);
		if (dist < res) {
			res = dist;
			*pos = i - k;
		}
	}
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

int dist_bit(BITSEQ *seq, size_t i, BITSEQ str, size_t k)
{
	seq += i / SEQSIZ;
	i %= SEQSIZ;

	// whether pattern contained in one BITSEQ
	if (i >= k) {
		// right justify and compare
		BITSEQ pat = (*seq >> 2*(SEQSIZ - i)) & (1ll << 2*k) - 1;
		str ^= pat;
	} else {
		// left shift end of *(seq-1) by i and add to right justified portion in *seq
		BITSEQ patl = (*(seq-1) << 2*i) & (1ll<<2*k)-1,
		       patr = (*seq >> 2*(SEQSIZ - i));
		str ^= (patl | patr) & (1ll << 2*k)-1;
	}

	// need popcount in base 4
	str = (str & EVENDIGITS) | ((str & (EVENDIGITS << 1)) >> 1);
	return __builtin_popcountll(str); // count nonzero bits in str
}

// print, for testing (could be used in main as well)

void printbits(FILE *f, BITSEQ str, size_t k)
{
	for (int i = 0; i < k; ++i)
		fputc(CHAR(str, i, k), f);
	// fputc('\t', f);
}
