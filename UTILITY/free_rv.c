#include "../real.h"
#include <stdlib.h>

void free_real_vector(real_t *v, int l)
{
	/*  Frees a real vector of range [l..u].  */

	free((char*) (v+l));
}

