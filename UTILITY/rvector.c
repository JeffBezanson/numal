#include "../real.h"
#include <stdlib.h>

real_t *allocate_real_vector(int l, int u)
{
	/*  Allocates a real vector of range [l..u].  */

	void system_error(char *);
	real_t *p;

	p=(real_t *)malloc((unsigned) (u-l+1)*sizeof(real_t));
	if (!p)
            system_error("Failure in allocate_real_vector().");
	return p-l;
}
