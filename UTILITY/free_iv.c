#include "../real.h"
#include <stdlib.h>

void free_integer_vector(int *v, int l)
{
	/*  Frees an integer vector of range [l..u].  */

	free((char*) (v+l));
}

