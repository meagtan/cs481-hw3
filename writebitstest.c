#include "medstr.h"
#include <stdio.h>
#include <string.h>

char *itoa(int, char *, int);

int main(void)
{
	char *seq = "TAAGTCTATACCATCGTAGTCTAATTAACGTTATGGTAGGAT";
	size_t siz;
	char str[100];

	BITSEQ *res = writebits(seq, &siz);

	for (int i = 0; i < 4; ++i)
		printf("%d\t%c\t%d\n", i, nucs[i], BITS(nucs[i]));

	printf("%lu\t%lu\n", strlen(seq), siz);

	printf("%s\n", seq);
	for (int i = 0; i < siz; ++i)
		// for (int j = 0; j < SEQSIZ; ++j)
		printf("%c", nucs[(res[i / SEQSIZ] >> (2 * SEQSIZ - 2 * (i % SEQSIZ) - 2)) & 3]);
	// for (int i = 0; i < siz; ++i) 
	// 	printf("%s%s", itoa(res[i] >> 16, str, 4), itoa(res[i] & ((1 << 16) - 1), str, 4));
	printf("\n");
}
