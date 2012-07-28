#include "../real.h"
#include <stdlib.h>

void free_real_matrix(real_t **m, int lr, int ur, int lc)
{
	/*  Frees a real matrix of range [lr..ur][lc..uc].  */

	int i;

	for (i=ur; i>=lr; i--) free((char*) (m[i]+lc));
	free((char*) (m+lr));
}

