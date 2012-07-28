#include "../real.h"
#include <stdlib.h>

void free_integer_matrix(int **m, int lr, int ur, int lc)
{
	/*  Frees an integer matrix of range [lr..ur][lc..uc].  */

	int i;

	for (i=ur; i>=lr; i--) free((char*) (m[i]+lc));
	free((char*) (m+lr));
}

